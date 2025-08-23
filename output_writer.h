#pragma once
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

// Forward declarations for logging functions
void logDebug(const std::string& message);
void logError(const std::string& message);
void logInfo(const std::string& message);
void logWarning(const std::string& message);

/**
 * Write aligned sequences to a FASTA file with generated headers.
 * @param aligned Aligned sequences (all of equal length)
 * @param filename Output file path
 * @throws std::runtime_error When there are no sequences to write or file cannot be created,
 *         or when sequence lengths are inconsistent
 */
inline void writeFasta(const std::vector<std::string>& aligned, const std::string& filename) {
    try {
        logDebug("Writing " + std::to_string(aligned.size()) + " sequences to " + filename);
        
        if (aligned.empty()) {
            logError("No sequences to write");
            throw std::runtime_error("No sequences to write");
        }
        
        std::ofstream out(filename);
        if (!out.is_open()) {
            logError("Cannot create output file: " + filename);
            throw std::runtime_error("Cannot create output file: " + filename);
        }
        
        // Validate sequence lengths
        if (!aligned.empty()) {
            size_t expected_length = aligned[0].length();
            for (size_t i = 1; i < aligned.size(); ++i) {
                if (aligned[i].length() != expected_length) {
                    logError("Sequence length mismatch: seq0=" + std::to_string(expected_length) + ", seq" + std::to_string(i) + "=" + std::to_string(aligned[i].length()));
                    throw std::runtime_error("Sequence length mismatch in output");
                }
            }
            logDebug("All " + std::to_string(aligned.size()) + " sequences have length " + std::to_string(expected_length));
        }
        
        for (size_t i = 0; i < aligned.size(); ++i) {
            out << ">Sequence_" << i + 1 << "\n" << aligned[i] << "\n";
        }
        
        out.close();
        logInfo("Successfully wrote alignment to " + filename);
    }
    catch (const std::exception& e) {
        logError("Error writing FASTA file: " + std::string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error writing FASTA file");
        throw;
    }
}

/**
 * Write a pairwise distance matrix to CSV as edge list rows (i<j).
 * @param matrix Symmetric NxN distances
 * @param names Sequence names for rows/columns
 * @param filename Output CSV file path
 * @throws std::runtime_error When the CSV file cannot be created or when matrix/names sizes mismatch
 */
inline void writeCSV(const std::vector<std::vector<double>>& matrix, const std::vector<std::string>& names, const std::string& filename) {
    try {
        logDebug("Writing distance matrix with names to " + filename);
        
        std::ofstream out(filename);
        if (!out.is_open()) {
            logError("Cannot create CSV file: " + filename);
            throw std::runtime_error("Cannot create CSV file: " + filename);
        }
        
        size_t N = matrix.size();
        if (N != names.size()) {
            logError("Matrix size (" + std::to_string(N) + ") doesn't match names size (" + std::to_string(names.size()) + ")");
            throw std::runtime_error("Matrix and names size mismatch");
        }
        
        out << "Seq1,Seq2,EMD\n";
        int pairCount = 0;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                out << names[i] << "," << names[j] << "," << matrix[i][j] << "\n"; // edge list
                pairCount++;
            }
        }
        
        out.close();
        logInfo("Successfully wrote " + std::to_string(pairCount) + " distance pairs to " + filename);
    }
    catch (const std::exception& e) {
        logError("Error writing CSV file: " + std::string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error writing CSV file");
        throw;
    }
}

/**
 * Write a pairwise distance matrix to CSV by numeric indices when names are unavailable.
 * @param matrix Symmetric NxN distances
 * @param filename Output CSV file path
 * @throws std::runtime_error When the CSV file cannot be created or when matrix is empty/non-square
 */
inline void writeCSV(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    try {
        logDebug("Writing distance matrix with indices to " + filename);
        
        std::ofstream out(filename);
        if (!out.is_open()) {
            logError("Cannot create CSV file: " + filename);
            throw std::runtime_error("Cannot create CSV file: " + filename);
        }
        
        size_t N = matrix.size();
        if (N == 0) {
            logError("Empty distance matrix");
            throw std::runtime_error("Empty distance matrix");
        }
        
        out << "i,j,EMD\n";
        int pairCount = 0;
        for (size_t i = 0; i < N; ++i) {
            if (matrix[i].size() != N) {
                logError("Non-square matrix: row " + std::to_string(i) + " has " + std::to_string(matrix[i].size()) + " columns, expected " + std::to_string(N));
                throw std::runtime_error("Non-square distance matrix");
            }
            for (size_t j = i + 1; j < N; ++j) {
                out << i << "," << j << "," << matrix[i][j] << "\n";
                pairCount++;
            }
        }
        
        out.close();
        logInfo("Successfully wrote " + std::to_string(pairCount) + " distance pairs to " + filename);
    }
    catch (const std::exception& e) {
        logError("Error writing CSV file: " + std::string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error writing CSV file");
        throw;
    }
}


