#pragma once
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

// Forward declarations for debug functions
void logDebug(const std::string& message);
void logError(const std::string& message);
void logInfo(const std::string& message);

/**
 * Read sequences and their header names from a FASTA file.
 * Lines starting with '>' denote headers; subsequent non-space characters
 * are concatenated as sequence letters. Supports files with UTF-8 BOM.
 *
 * @param filename Path to the FASTA file
 * @return pair of (sequences, names) where names[i] belongs to sequences[i]
 * @throws std::runtime_error When the file cannot be opened, contains no sequences,
 *         contains empty sequences, or when I/O errors occur
 */
inline std::pair<std::vector<std::string>, std::vector<std::string>> readFastaWithNames(const std::string& filename) {
    try {
        logDebug("Opening FASTA file: " + filename);
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            logError("Cannot open FASTA file: " + filename);
            throw std::runtime_error("Cannot open FASTA file: " + filename);
        }
        
        std::vector<std::string> sequences;
        std::vector<std::string> names;
        std::string line, currentSeq, currentName;
        int lineCount = 0;
        int sequenceCount = 0;

        while (std::getline(file, line)) {
            lineCount++;
            if (line.empty()) continue;

            // Strip UTF-8 BOM if present
            if (line.compare(0, 3, "\xEF\xBB\xBF") == 0) {
                line = line.substr(3);
                logDebug("Stripped UTF-8 BOM from line " + std::to_string(lineCount));
            }

            if (line[0] == '>') {
                if (!currentSeq.empty()) {
                    sequences.push_back(currentSeq);
                    names.push_back(currentName);
                    sequenceCount++;
                    logDebug("Processed sequence " + std::to_string(sequenceCount) + " (" + currentName + ") - length: " + std::to_string(currentSeq.length()));
                    currentSeq.clear();
                }
                currentName = line.substr(1);  // drop '>'
                if (currentName.empty()) {
                    logDebug("Empty sequence name found at line " + std::to_string(lineCount));
                    currentName = "Unnamed_" + std::to_string(sequenceCount + 1);
                }
            } else {
                for (char c : line) {
                    if (!std::isspace(c)) {
                        currentSeq += c; // ignore whitespace within sequences
                    }
                }
            }
        }

        // Don't forget the last sequence
        if (!currentSeq.empty()) {
            sequences.push_back(currentSeq);
            names.push_back(currentName);
            sequenceCount++;
            logDebug("Processed final sequence " + std::to_string(sequenceCount) + " (" + currentName + ") - length: " + std::to_string(currentSeq.length()));
        }

        file.close();
        
        if (sequences.empty()) {
            logError("No sequences found in FASTA file");
            throw std::runtime_error("No sequences found in FASTA file");
        }
        
        logInfo("Successfully loaded " + std::to_string(sequences.size()) + " sequences from " + filename);
        
        // Validate sequences
        for (size_t i = 0; i < sequences.size(); ++i) {
            if (sequences[i].empty()) {
                logError("Empty sequence found: " + names[i]);
                throw std::runtime_error("Empty sequence found: " + names[i]);
            }
            
            // Check for invalid characters (basic validation)
            for (char c : sequences[i]) {
                if (c != '-' && (c < 'A' || c > 'Z') && (c < 'a' || c > 'z')) {
                    logDebug("Warning: Non-standard character '" + std::string(1, c) + "' in sequence " + names[i]);
                }
            }
        }

        return {sequences, names};
    }
    catch (const std::exception& e) {
        logError("Error reading FASTA file: " + std::string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error reading FASTA file");
        throw;
    }
}

