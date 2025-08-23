#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>

/**
 * Amino-acid alphabet used throughout (20 residues + gap '-')
 */
extern const std::vector<char> AA;
extern const int AA_SIZE;

constexpr double EPSILON = 1e-4;

/**
 * Construct a BLOSUM-inspired log-odds cost matrix over the AA alphabet.
 * The matrix is derived from background frequencies and simple pair modeling,
 * including a gap symbol '-'. Values approximate dissimilarity (higher = less similar).
 * @return AA_SIZE x AA_SIZE matrix of log-odds costs
 */
std::vector<std::vector<double>> getRealBlosumLogOddsMatrix();
