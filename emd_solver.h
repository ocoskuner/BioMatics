#pragma once
#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include "blosum62.h"

double calculateEMD(const std::vector<double>& dist1, const std::vector<double>& dist2) {
    const int N = 21;
    if (dist1.size() != N || dist2.size() != N)
        throw std::invalid_argument("Distributions must be size 21");

    std::vector<double> f1 = dist1, f2 = dist2;
    double sum1 = std::accumulate(f1.begin(), f1.end(), 0.0);
    double sum2 = std::accumulate(f2.begin(), f2.end(), 0.0);
    if (sum1 == 0 || sum2 == 0) return 1e6;
    for (auto& v : f1) v /= sum1;
    for (auto& v : f2) v /= sum2;

    auto cost = getRealBlosumLogOddsMatrix();  

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            cost[i][j] = std::max(cost[i][j], -3.5);  

    double total_cost = 0.0;
    std::vector<double> supply = f1;
    std::vector<double> demand = f2;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double flow = std::min(supply[i], demand[j]);
            total_cost += flow * cost[i][j];
            supply[i] -= flow;
            demand[j] -= flow;
        }
    }

    return std::tanh(total_cost);
}

