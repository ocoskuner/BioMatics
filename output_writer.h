#pragma once
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <cctype>

// Forward declarations for logging functions
void logDebug(const std::string& message);
void logError(const std::string& message);
void logInfo(const std::string& message);
void logWarning(const std::string& message);

/**
 * Write aligned sequences to a FASTA file using provided names as headers.
 * Falls back to generated names if sizes mismatch.
 * @param aligned Aligned sequences (all of equal length)
 * @param names Headers for each aligned sequence (same order and size as aligned)
 * @param filename Output file path
 * @throws std::runtime_error When there are no sequences to write or file cannot be created,
 *         or when sequence lengths are inconsistent
 */
inline void writeFasta(const std::vector<std::string>& aligned, const std::vector<std::string>& names, const std::string& filename) {
    try {
        logDebug("Writing " + std::to_string(aligned.size()) + " sequences to " + filename + " with provided names");
        
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
        size_t expected_length = aligned[0].length();
        for (size_t i = 1; i < aligned.size(); ++i) {
            if (aligned[i].length() != expected_length) {
                logError("Sequence length mismatch: seq0=" + std::to_string(expected_length) + ", seq" + std::to_string(i) + "=" + std::to_string(aligned[i].length()));
                throw std::runtime_error("Sequence length mismatch in output");
            }
        }
        logDebug("All " + std::to_string(aligned.size()) + " sequences have length " + std::to_string(expected_length));

        if (names.size() != aligned.size()) {
            logWarning("Names size (" + std::to_string(names.size()) + ") does not match aligned size (" + std::to_string(aligned.size()) + "). Falling back to generated headers.");
            for (size_t i = 0; i < aligned.size(); ++i) {
                out << ">Sequence_" << i + 1 << "\n" << aligned[i] << "\n";
            }
        } else {
            for (size_t i = 0; i < aligned.size(); ++i) {
                out << ">" << names[i] << "\n" << aligned[i] << "\n";
            }
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
 * Write aligned sequences to A3M format (HH-suite compatible).
 * Rules:
 *  - Use the first sequence as master.
 *  - Columns where master has a residue are MATCH columns: residues are uppercase; '-' kept as deletions.
 *  - Columns where master has '-' are INSERT columns: residues in other sequences become lowercase letters; '-' omitted.
 * @param aligned Aligned sequences (all of equal length)
 * @param filename Output file path (e.g., final_msa.a3m)
 */
inline void writeA3M(const std::vector<std::string>& aligned, const std::string& filename) {
    try {
        logDebug("Writing A3M for " + std::to_string(aligned.size()) + " sequences to " + filename);

        if (aligned.empty()) {
            logError("No sequences to write (A3M)");
            throw std::runtime_error("No sequences to write (A3M)");
        }

        const size_t L = aligned[0].size();
        for (size_t i = 1; i < aligned.size(); ++i) {
            if (aligned[i].size() != L) {
                logError("Sequence length mismatch in A3M: row 0 has " + std::to_string(L) + ", row " + std::to_string(i) + " has " + std::to_string(aligned[i].size()));
                throw std::runtime_error("Sequence length mismatch in A3M");
            }
        }

        const std::string& master = aligned[0];
        std::vector<bool> isMatch(L, false);
        for (size_t i = 0; i < L; ++i) {
            isMatch[i] = (master[i] != '-');
        }

        std::ofstream out(filename);
        if (!out.is_open()) {
            logError("Cannot create A3M file: " + filename);
            throw std::runtime_error("Cannot create A3M file: " + filename);
        }

        for (size_t idx = 0; idx < aligned.size(); ++idx) {
            const std::string& src = aligned[idx];
            std::string line;
            line.reserve(src.size());

            for (size_t i = 0; i < L; ++i) {
                char c = src[i];
                if (isMatch[i]) {
                    if (c == '-') {
                        line.push_back('-'); // deletion relative to master
                    } else {
                        line.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
                    }
                } else {
                    if (c != '-') {
                        line.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
                    }
                }
            }

            out << ">Sequence_" << (idx + 1) << "\n" << line << "\n";
        }

        out.close();
        logInfo("Successfully wrote A3M to " + filename);
    }
    catch (const std::exception& e) {
        logError("Error writing A3M file: " + std::string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error writing A3M file");
        throw;
    }
}

/**
 * Write aligned sequences to A3M format using provided names as headers.
 * Falls back to generated names if sizes mismatch.
 * @param aligned Aligned sequences
 * @param names Headers for each aligned sequence (same size/order as aligned)
 * @param filename Output file path (e.g., final_msa.a3m)
 */
inline void writeA3M(const std::vector<std::string>& aligned, const std::vector<std::string>& names, const std::string& filename) {
    try {
        logDebug("Writing A3M (with names) for " + std::to_string(aligned.size()) + " sequences to " + filename);

        if (aligned.empty()) {
            logError("No sequences to write (A3M)");
            throw std::runtime_error("No sequences to write (A3M)");
        }

        const size_t L = aligned[0].size();
        for (size_t i = 1; i < aligned.size(); ++i) {
            if (aligned[i].size() != L) {
                logError("Sequence length mismatch in A3M: row 0 has " + std::to_string(L) + ", row " + std::to_string(i) + " has " + std::to_string(aligned[i].size()));
                throw std::runtime_error("Sequence length mismatch in A3M");
            }
        }

        const std::string& master = aligned[0];
        std::vector<bool> isMatch(L, false);
        for (size_t i = 0; i < L; ++i) {
            isMatch[i] = (master[i] != '-');
        }

        std::ofstream out(filename);
        if (!out.is_open()) {
            logError("Cannot create A3M file: " + filename);
            throw std::runtime_error("Cannot create A3M file: " + filename);
        }

        bool useProvided = (names.size() == aligned.size());
        if (!useProvided) {
            logWarning("Names size (" + std::to_string(names.size()) + ") does not match aligned size (" + std::to_string(aligned.size()) + "). Falling back to generated headers for A3M.");
        }

        for (size_t idx = 0; idx < aligned.size(); ++idx) {
            const std::string& src = aligned[idx];
            std::string line;
            line.reserve(src.size());

            for (size_t i = 0; i < L; ++i) {
                char c = src[i];
                if (isMatch[i]) {
                    if (c == '-') {
                        line.push_back('-');
                    } else {
                        line.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
                    }
                } else {
                    if (c != '-') {
                        line.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
                    }
                }
            }

            if (useProvided) {
                out << ">" << names[idx] << "\n" << line << "\n";
            } else {
                out << ">Sequence_" << (idx + 1) << "\n" << line << "\n";
            }
        }

        out.close();
        logInfo("Successfully wrote A3M to " + filename);
    }
    catch (const std::exception& e) {
        logError("Error writing A3M file: " + std::string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error writing A3M file");
        throw;
    }
}

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


