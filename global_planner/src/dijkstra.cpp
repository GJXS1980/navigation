/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2008, 2013, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: Eitan Marder-Eppstein
 *         David V. Lu!!
 *********************************************************************/
#include<global_planner/dijkstra.h>
#include <algorithm>
namespace global_planner {

DijkstraExpansion::DijkstraExpansion(PotentialCalculator* p_calc, int nx, int ny) :
        Expander(p_calc, nx, ny), pending_(NULL), precise_(false) {
    // 优先缓冲区
    buffer1_ = new int[PRIORITYBUFSIZE];
    buffer2_ = new int[PRIORITYBUFSIZE];
    buffer3_ = new int[PRIORITYBUFSIZE];

    priorityIncrement_ = 2 * neutral_cost_;
}

DijkstraExpansion::~DijkstraExpansion() {
  delete[] buffer1_;
  delete[] buffer2_;
  delete[] buffer3_;
  if (pending_)
      delete[] pending_;
}

//
// 设置/处置地图尺寸
//
void DijkstraExpansion::setSize(int xs, int ys) {
    Expander::setSize(xs, ys);
    if (pending_)
        delete[] pending_;

    pending_ = new bool[ns_];
    memset(pending_, 0, ns_ * sizeof(bool));
}

//
// main propagation function
// Dijkstra method, breadth-first
// runs for a specified number of cycles,
//   or until it runs out of cells to update,
//   or until the Start cell is found (atStart = true)
//
/*********************************************************************
 *
输入参数：
 * costs: 每个网格单元格的代价，存储在一维数组中
 * start_x,start_y:起始位置的坐标
 * end_x,end_y:目标位置的坐标
 * cycles:算法运行的最大轮数
 * potential:每个单元格的潜力值,存储在一维数组中

输出参数：
 * bool类型,表示是否成功找到路径
 * 函数中的变量和数据结构：

 * cells_visited_:被访问的单元格数量
 * lethal_cost_:表示代表不可行走区域的代价
 * threshold_:表示当前轮次的代价阈值
 * buffer1_,buffer2_,buffer3_:存储扩展的单元格的优先级缓冲区
 * currentBuffer_,currentEnd_:表示当前轮次的优先级缓冲区和缓冲区中的单元格数量
 * nextBuffer_,nextEnd_:表示下一轮次的优先级缓冲区和缓冲区中的单元格数量
 * overBuffer_,overEnd_:表示溢出的优先级缓冲区和缓冲区中的单元格数量
 * pending_:一个bool类型的数组,记录单元格是否在优先级缓冲区中
 * potential:每个单元格的潜力值
 * precise_:一个bool类型的变量,表示算法是否精确计算单元格的潜力值
 * nx_:地图中单行单元格的数量

主要算法流程：
 * 初始化优先级缓冲区和单元格的潜力值。
 * 将起始单元格加入当前优先级缓冲区，并根据精度计算其潜力值。
 * 进行一定数量的循环，直到找到路径或达到指定的轮次：
 * 从当前优先级缓冲区中取出单元格，并更新其周围的单元格的潜力值。
 * 将被更新的单元格加入下一轮次的优先级缓冲区中。
 * 如果当前优先级缓冲区中的单元格已经全部更新，则增加代价阈值，将下一轮次的优先级缓冲区设置为溢出的优先级缓冲区，并交换当前和下一轮的优先级缓冲区。
 * 如果目标单元格的潜力值已经被更新，则表示找到路径。
 * 返回算法是否成功找到路径的结果。
 *       
 *********************************************************************/
bool DijkstraExpansion::calculatePotentials(unsigned char* costs, double start_x, double start_y, double end_x, double end_y,
                                           int cycles, float* potential) {
    cells_visited_ = 0; // 记录访问的单元格数量
    // priority buffers
    threshold_ = lethal_cost_;  // 优先级阈值初始值为致命代价值
    currentBuffer_ = buffer1_;  // 当前优先级缓冲区
    currentEnd_ = 0;            // 当前优先级缓冲区已存储元素个数
    nextBuffer_ = buffer2_; // 下一个优先级缓冲区
    nextEnd_ = 0;           // 下一个优先级缓冲区已存储元素个数
    overBuffer_ = buffer3_; // 溢出优先级缓冲区
    overEnd_ = 0;           // 溢出优先级缓冲区已存储元素个数
    memset(pending_, 0, ns_ * sizeof(bool));    // 初始化 pending_ 数组，表示单元格是否在优先级队列中
    std::fill(potential, potential + ns_, POT_HIGH);     // 初始化 potential 数组，表示各单元格的潜力值

     // 设置起点
    int k = toIndex(start_x, start_y);  // 起点在 costmap 中的下标

    if(precise_)    // 如果需要精确计算
    {
        double dx = start_x - (int)start_x, dy = start_y - (int)start_y;    // 计算 start_x 和 start_y 与其整数部分之差
        dx = floorf(dx * 100 + 0.5) / 100;  // 向下取整，// 对小数部分保留两位，四舍五入到小数点后第二位
        dy = floorf(dy * 100 + 0.5) / 100;
        potential[k] = neutral_cost_ * 2 * dx * dy; // 计算起点潜力值
        potential[k+1] = neutral_cost_ * 2 * (1-dx)*dy;
        potential[k+nx_] = neutral_cost_*2*dx*(1-dy);
        potential[k+nx_+1] = neutral_cost_*2*(1-dx)*(1-dy);//*/

        push_cur(k+2);  // 将起点周围的八个单元格加入当前优先级缓冲区
        push_cur(k-1);
        push_cur(k+nx_-1);
        push_cur(k+nx_+2);

        push_cur(k-nx_);
        push_cur(k-nx_+1);
        push_cur(k+nx_*2);
        push_cur(k+nx_*2+1);
    }else{                   // 如果不需要精确计算
        potential[k] = 0;
        push_cur(k+1);
        push_cur(k-1);
        push_cur(k-nx_);
        push_cur(k+nx_);
    }

    int nwv = 0;            // 最大优先级块大小
    int nc = 0;            // 放入优先级块的单元格数量
    int cycle = 0;        // which cycle we're on

    // 设置起始单元格
    int startCell = toIndex(end_x, end_y);

    for (; cycle < cycles; cycle++) // 进行这么多次循环，除非被打断
            {
        // 
        if (currentEnd_ == 0 && nextEnd_ == 0) // 优先缓冲区为空
            return false;

        // stats
        nc += currentEnd_;
        if (currentEnd_ > nwv)
            nwv = currentEnd_;

        // reset pending_ flags on current priority buffer
        // 重置当前优先级缓冲区上的 pending_ 标志
        int *pb = currentBuffer_;
        int i = currentEnd_;
        while (i-- > 0)
            pending_[*(pb++)] = false;

        // process current priority buffer
        pb = currentBuffer_;
        i = currentEnd_;
        while (i-- > 0)
            updateCell(costs, potential, *pb++);

        // 交换优先缓冲区 currentBuffer_ <=> nextBuffer_
        currentEnd_ = nextEnd_;
        nextEnd_ = 0;
        pb = currentBuffer_;        // 交换缓冲区
        currentBuffer_ = nextBuffer_;
        nextBuffer_ = pb;

        // 看看我们是否完成了这个优先级交换
        if (currentEnd_ == 0) {
            threshold_ += priorityIncrement_;    // 增加优先级阈值
            currentEnd_ = overEnd_;    // 将当前设置为溢出块
            overEnd_ = 0;
            pb = currentBuffer_;        // 交换缓冲区
            currentBuffer_ = overBuffer_;
            overBuffer_ = pb;
        }

        // 检查我们是否点击了开始单元格
        if (potential[startCell] < POT_HIGH)
            break;
    }
    //ROS_INFO("CYCLES %d/%d ", cycle, cycles);
    if (cycle < cycles)
        return true; // 完成
    else
        return false;
}

//
// Critical function: calculate updated potential value of a cell,
//   given its neighbors' values
// Planar-wave update calculation from two lowest neighbors in a 4-grid
// Quadratic approximation to the interpolated value
// No checking of bounds here, this function should be fast
//

#define INVSQRT2 0.707106781

// 更新网格图中的一个格点,用于更新单个格点的最短路径的权重
inline void DijkstraExpansion::updateCell(unsigned char* costs, float* potential, int n) {
    cells_visited_++;   // 更新被访问格点的计数器

    // do planar wave update
     // 计算当前格点的代价值
    float c = getCost(costs, n);
    // 如果当前格点的代价值大于等于障碍物的代价值，不进行更新
    if (c >= lethal_cost_)    // don't propagate into obstacles
        return;

    // 计算当前格点的代价值
    float pot = p_calc_->calculatePotential(potential, c, n);

    // // 将受影响的邻居添加到优先级队列中
    if (pot < potential[n]) {   // 如果当前格点的代价值比之前的小，说明需要更新
        float le = INVSQRT2 * (float)getCost(costs, n - 1);
        float re = INVSQRT2 * (float)getCost(costs, n + 1);
        float ue = INVSQRT2 * (float)getCost(costs, n - nx_);
        float de = INVSQRT2 * (float)getCost(costs, n + nx_);
        potential[n] = pot; // 更新当前格点的势能值

        //ROS_INFO("UPDATE %d %d %d %f", n, n%nx, n/nx, potential[n]);
        // 如果当前格点的代价值小于等于threshold,加入低成本缓存块
        if (pot < threshold_)    // low-cost buffer block(低成本缓冲块)
                {
            if (potential[n - 1] > pot + le)
                push_next(n-1);
            if (potential[n + 1] > pot + re)
                push_next(n+1);
            if (potential[n - nx_] > pot + ue)
                push_next(n-nx_);
            if (potential[n + nx_] > pot + de)
                push_next(n+nx_);
        } else            // overflow block(溢出块)
        {
            if (potential[n - 1] > pot + le)
                push_over(n-1);
            if (potential[n + 1] > pot + re)
                push_over(n+1);
            if (potential[n - nx_] > pot + ue)
                push_over(n-nx_);
            if (potential[n + nx_] > pot + de)
                push_over(n+nx_);
        }
    }
}

} //end namespace global_planner
