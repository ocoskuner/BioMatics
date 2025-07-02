#include "blosum62.h"

const std::vector<char> AA = {
    'A','C','D','E','F','G','H','I','K','L','M',
    'N','P','Q','R','S','T','V','W','Y','-'
};

const int AA_SIZE = 21;

std::vector<std::vector<double>> getRealBlosumLogOddsMatrix() {
    std::unordered_map<char, double> P = {
        {'A', 0.078}, {'C', 0.019}, {'D', 0.053}, {'E', 0.062},
        {'F', 0.040}, {'G', 0.073}, {'H', 0.023}, {'I', 0.053},
        {'K', 0.058}, {'L', 0.091}, {'M', 0.023}, {'N', 0.044},
        {'P', 0.052}, {'Q', 0.040}, {'R', 0.052}, {'S', 0.071},
        {'T', 0.058}, {'V', 0.065}, {'W', 0.014}, {'Y', 0.034},
        {'-', 0.020}
    };

    std::vector<std::vector<double>> P_ij(AA_SIZE, std::vector<double>(AA_SIZE));
    std::vector<std::vector<double>> cost(AA_SIZE, std::vector<double>(AA_SIZE));

    for (int i = 0; i < AA_SIZE; ++i) {
        for (int j = 0; j < AA_SIZE; ++j) {
            char a = AA[i], b = AA[j];

            if (a == b)
                P_ij[i][j] = P[a] * P[b] * 2.0;
            else if (a == '-' || b == '-')
                P_ij[i][j] = P[a] * P[b] * 0.2;
            else
                P_ij[i][j] = P[a] * P[b];

            double denom = P[a] * P[b];
            double ratio = (P_ij[i][j] + EPSILON) / (denom + EPSILON);
            cost[i][j] = -log2(ratio);
        }
    }

    return cost;
}
