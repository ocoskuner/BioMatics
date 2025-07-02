#pragma once
#include <fstream>
#include <vector>
#include <string>

inline void writeFasta(const std::vector<std::string>& aligned, const std::string& filename) {
    std::ofstream out(filename);
    for (size_t i = 0; i < aligned.size(); ++i) {
        out << ">Sequence_" << i + 1 << "\n" << aligned[i] << "\n";
    }
}

inline void writeCSV(const std::vector<std::vector<double>>& matrix, const std::vector<std::string>& names, const std::string& filename) {
    std::ofstream out(filename);
    size_t N = matrix.size();
    out << "Seq1,Seq2,EMD\n";
    for (size_t i = 0; i < N; ++i)
        for (size_t j = i + 1; j < N; ++j)
            out << names[i] << "," << names[j] << "," << matrix[i][j] << "\n";
}


