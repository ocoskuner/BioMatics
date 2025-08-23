#include "blosum62.h"
#include <iostream>

// Forward declarations for logging functions
void logDebug(const std::string& message);
void logError(const std::string& message);
void logInfo(const std::string& message);

/**
 * Global amino-acid alphabet: 20 standard residues plus a gap symbol '-'.
 */
const std::vector<char> AA = {
    'A','C','D','E','F','G','H','I','K','L','M',
    'N','P','Q','R','S','T','V','W','Y','-'
};

const int AA_SIZE = 21;

/**
 * Build a simple log-odds cost matrix inspired by BLOSUM-like modeling.
 * Increases self-pair probability, penalizes gap-interactions, then converts
 * to a log-odds cost: cost = -log2(P_ij / (P_i * P_j)).
 *
 * @return AA_SIZE x AA_SIZE matrix of costs (higher implies less likely/similar)
 */
std::vector<std::vector<double>> getRealBlosumLogOddsMatrix() {
    try {
        // Removed verbose debug logging to reduce output noise
        std::unordered_map<char, double> P = {
            {'A', 0.078}, {'C', 0.019}, {'D', 0.053}, {'E', 0.062},
            {'F', 0.040}, {'G', 0.073}, {'H', 0.023}, {'I', 0.053},
            {'K', 0.058}, {'L', 0.091}, {'M', 0.023}, {'N', 0.044},
            {'P', 0.052}, {'Q', 0.040}, {'R', 0.052}, {'S', 0.071},
            {'T', 0.058}, {'V', 0.065}, {'W', 0.014}, {'Y', 0.034},
            {'-', 0.020}
        };

        // Validate that all amino acids have probabilities
        for (int i = 0; i < AA_SIZE; ++i) {
            char aa = AA[i];
            if (P.find(aa) == P.end()) {
                logError(std::string("Missing probability for amino acid: ") + std::string(1, aa));
                throw std::runtime_error("Missing amino acid probability");
            }
        }

        std::vector<std::vector<double>> P_ij(AA_SIZE, std::vector<double>(AA_SIZE));
        std::vector<std::vector<double>> cost(AA_SIZE, std::vector<double>(AA_SIZE));

        for (int i = 0; i < AA_SIZE; ++i) {
            for (int j = 0; j < AA_SIZE; ++j) {
                try {
                    char a = AA[i], b = AA[j];

                    if (a == b)
                        P_ij[i][j] = P[a] * P[b] * 2.0;   // boost self-pairs
                    else if (a == '-' || b == '-')
                        P_ij[i][j] = P[a] * P[b] * 0.2;   // penalize gap interactions
                    else
                        P_ij[i][j] = P[a] * P[b];         // baseline co-occurrence

                    double denom = P[a] * P[b];
                    double ratio = (P_ij[i][j] + EPSILON) / (denom + EPSILON);
                    cost[i][j] = -log2(ratio);            // convert to log-odds cost
                    
                    // Validate the computed cost
                    if (std::isnan(cost[i][j]) || std::isinf(cost[i][j])) {
                        logError(std::string("Invalid cost computed for (") + std::string(1, a) + "," + std::string(1, b) + "): " + std::to_string(cost[i][j]));
                        cost[i][j] = 1.0; // fallback value
                    }
                }
                catch (const std::exception& e) {
                    logError("Error computing cost for (" + std::to_string(i) + "," + std::to_string(j) + "): " + e.what());
                    cost[i][j] = 1.0; // fallback value
                }
            }
        }
        
        // Log some sample values for verification
        
        return cost;
    }
    catch (const std::exception& e) {
        logError("Error in getRealBlosumLogOddsMatrix: " + std::string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error in getRealBlosumLogOddsMatrix");
        throw;
    }
}
