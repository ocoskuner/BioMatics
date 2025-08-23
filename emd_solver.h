#pragma once
#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include "blosum62.h"

// Forward declarations for logging functions
void logDebug(const std::string& message);
void logError(const std::string& message);
void logInfo(const std::string& message);

// Sentinel value used to indicate an error/invalid distance
constexpr double EMD_ERROR = 1e6;
/**
 * Compute a simple Earth Mover's Distance (EMD) between two categorical
 * distributions over the amino-acid alphabet using a BLOSUM-derived cost.
 * Distributions are normalized internally. Uses a greedy transport and
 * returns tanh(total_cost) as a dissimilarity proxy.
 *
 * Contract:
 *  - Typical range is approximately [-1, 1); negative values may occur when
 *    BLOSUM-derived costs are negative for highly similar residues.
 *  - On error or invalid inputs, the function returns the sentinel EMD_ERROR (=1e6).
 *
 * @param dist1 First distribution (size 21; non-negative; can be unnormalized)
 * @param dist2 Second distribution (size 21; non-negative; can be unnormalized)
 * @return Dissimilarity proxy in roughly [-1, 1) or EMD_ERROR on failure
 * @throws std::invalid_argument if inputs are not size 21
 */
double calculateEMD(const std::vector<double>& dist1, const std::vector<double>& dist2) {
    try {
        const int N = 21;
        if (dist1.size() != N || dist2.size() != N) {
            logError("Invalid distribution sizes: dist1=" + std::to_string(dist1.size()) + ", dist2=" + std::to_string(dist2.size()) + ", expected=" + std::to_string(N));
            throw std::invalid_argument("Distributions must be size 21");
        }

        std::vector<double> f1 = dist1, f2 = dist2; // copy for normalization
        double sum1 = std::accumulate(f1.begin(), f1.end(), 0.0);
        double sum2 = std::accumulate(f2.begin(), f2.end(), 0.0);
        
        if (sum1 == 0) {
            logError("First distribution has zero sum");
            logDebug("Fallback used: returning EMD_ERROR due to zero-sum first distribution");
            return EMD_ERROR;
        }
        if (sum2 == 0) {
            logError("Second distribution has zero sum");
            logDebug("Fallback used: returning EMD_ERROR due to zero-sum second distribution");
            return EMD_ERROR;
        }
        
        // Normalize distributions
        for (auto& v : f1) v /= sum1;
        for (auto& v : f2) v /= sum2;

        auto cost = getRealBlosumLogOddsMatrix();  // dissimilarity proxy

        // Clamp cost values to prevent extreme values and neutralize invalids
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                cost[i][j] = std::max(cost[i][j], -3.5);
                if (std::isnan(cost[i][j]) || std::isinf(cost[i][j])) {
                    logError("Invalid cost value at (" + std::to_string(i) + "," + std::to_string(j) + ")");
                    cost[i][j] = 0.0; // neutral fallback
                    logDebug("Fallback used: neutral cost applied at (" + std::to_string(i) + "," + std::to_string(j) + ")");
                }
            }
        }

        double total_cost = 0.0;
        std::vector<double> supply = f1;
        std::vector<double> demand = f2;

        // Greedy EMD computation
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                // Greedy transport from i -> j (no allocation/throw expected here)
                double flow = std::min(supply[i], demand[j]);
                if (flow > 0) {
                    total_cost += flow * cost[i][j];
                    supply[i] -= flow;
                    demand[j] -= flow;
                }
            }
        }

        // Check for numerical issues
        if (std::isnan(total_cost) || std::isinf(total_cost)) {
            logError("Invalid total cost in EMD: " + std::to_string(total_cost));
            logDebug("Fallback used: returning EMD_ERROR due to invalid total cost");
            return EMD_ERROR;
        }

        double result = std::tanh(total_cost); // bound to (0,1)
        
        if (std::isnan(result) || std::isinf(result)) {
            logError("Invalid EMD result: " + std::to_string(result));
            logDebug("Fallback used: returning EMD_ERROR due to invalid result");
            return EMD_ERROR;
        }

        return result;
    }
    catch (const std::exception& e) {
        logError("Exception in calculateEMD: " + std::string(e.what()));
        logDebug("Fallback used: returning EMD_ERROR due to exception");
        return EMD_ERROR; // return high dissimilarity on error
    }
    catch (...) {
        logError("Unknown error in calculateEMD");
        logDebug("Fallback used: returning EMD_ERROR due to unknown error");
        return EMD_ERROR;
    }
}

