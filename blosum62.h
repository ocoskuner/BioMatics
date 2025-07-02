#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>

extern const std::vector<char> AA;
extern const int AA_SIZE;

constexpr double EPSILON = 1e-4;

std::vector<std::vector<double>> getRealBlosumLogOddsMatrix();
