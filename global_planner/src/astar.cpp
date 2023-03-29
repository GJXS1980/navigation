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
#include<global_planner/astar.h>
#include<costmap_2d/cost_values.h>

namespace global_planner {
// 初始化 A* 扩展器
AStarExpansion::AStarExpansion(PotentialCalculator* p_calc, int xs, int ys) :
        Expander(p_calc, xs, ys) {
}

// 计算潜在代价值的函数(潜在代价是使用PotentialCalculator类的calculatePotential函数计算的)
// 输入：起始位置、目标位置、代价地图和运行次数
// 返回：布尔值，指示是否找到路径
bool AStarExpansion::calculatePotentials(unsigned char* costs, double start_x, double start_y, double end_x, double end_y,
                                        int cycles, float* potential) {
    queue_.clear(); //  清空队列
    int start_i = toIndex(start_x, start_y);    // 计算起点的索引
    queue_.push_back(Index(start_i, 0));    // 将起点加入队列

    std::fill(potential, potential + ns_, POT_HIGH);    // 将代价数组中所有元素初始化为 POT_HIGH
    potential[start_i] = 0; // 将起点的代价值设为 0

    // 索引计算公式:index = x + nx_ * y,其中nx:The x size of the map
    int goal_i = toIndex(end_x, end_y); // 计算终点的索引
    int cycle = 0;

    //   A* 搜索循环，直到找到终点或超出指定搜索次数
    while (queue_.size() > 0 && cycle < cycles) {
        Index top = queue_[0];  // 取出队列的第一个元素
        std::pop_heap(queue_.begin(), queue_.end(), greater1());    // 弹出堆顶元素
        queue_.pop_back();  // 将队列的最后一个元素弹出

        int i = top.i;
        if (i == goal_i)
            return true;
        // 将相邻的四个点添加到队列中
        add(costs, potential, potential[i], i + 1, end_x, end_y);
        add(costs, potential, potential[i], i - 1, end_x, end_y);
        add(costs, potential, potential[i], i + nx_, end_x, end_y);
        add(costs, potential, potential[i], i - nx_, end_x, end_y);

        cycle++;
    }

    return false;
}

//  添加点到队列的函数(计算相邻单元格的潜在代价，，并将其添加到队列中)
void AStarExpansion::add(unsigned char* costs, float* potential, float prev_potential, int next_i, int end_x,
                         int end_y) {
    //  如果超出索引范围则返回(ns = nx * ny,在potential_calculator.h中setSize函数求解)
    if (next_i < 0 || next_i >= ns_)
        return;
    //  如果代价值已经计算过则返回
    if (potential[next_i] < POT_HIGH)
        return;
    // 如果代价超过内切障碍值(253),则返回(在expander.h定义)
    // costmap_2d::NO_INFORMATION为未知空间Unknown Space(255)，unknown_(true)
    if(costs[next_i]>=lethal_cost_ && !(unknown_ && costs[next_i]==costmap_2d::NO_INFORMATION))
        return;

    potential[next_i] = p_calc_->calculatePotential(potential, costs[next_i] + neutral_cost_, next_i, prev_potential);
    int x = next_i % nx_, y = next_i / nx_;
    // 采用曼哈顿距离求解
    float distance = abs(end_x - x) + abs(end_y - y);

    //  计算F(n) = G(n) + H(n)
    queue_.push_back(Index(next_i, potential[next_i] + distance * neutral_cost_));
    std::push_heap(queue_.begin(), queue_.end(), greater1());
}

} //end namespace global_planner
