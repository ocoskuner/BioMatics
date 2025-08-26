/**
 * Multiple Sequence Alignment (MSA) tool for protein sequences.
 *
 * Pipeline:
 *  - Read sequences from FASTA with names
 *  - Build pairwise distance matrix using per-column distributions and EMD
 *  - Construct a UPGMA-like guide tree
 *  - Progressive profile–profile alignment with adaptive gap penalties
 *  - Iterative refinement and output final alignment
 */
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <unordered_map>
#include <stdexcept>
#include <iomanip>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "blosum62.h"
#include "emd_solver.h"
#include "fasta_reader.h"
#include "output_writer.h"

using namespace std;

// Global log file stream
ofstream logFile;

// Debug utility functions with file logging
void initializeLogFile() {
    logFile.open("msa_debug.log", ios::out | ios::trunc);
    if (logFile.is_open()) {
        logFile << "=== MSA Program Debug Log ===" << endl;
        logFile << "Started at: " << __DATE__ << " " << __TIME__ << endl;
        logFile << "=============================" << endl << endl;
    }
}

void closeLogFile() {
    if (logFile.is_open()) {
        logFile << endl << "=== End of Log ===" << endl;
        logFile.close();
    }
}

void logDebug(const std::string& message) {
    // Write to file (detailed)
    if (logFile.is_open()) {
        logFile << "[DEBUG] " << message << endl;
        logFile.flush();
    }
    // Don't print to console to reduce noise
}

void logError(const std::string& message) {
    // Write to file
    if (logFile.is_open()) {
        logFile << "[ERROR] " << message << endl;
        logFile.flush();
    }
    // Also print to console for critical errors
    std::cerr << "[ERROR] " << message << std::endl;
}

void logWarning(const std::string& message) {
    if (logFile.is_open()) {
        logFile << "[WARNING] " << message << endl;
        logFile.flush();
    }
    std::cerr << "[WARNING] " << message << std::endl;
}

void logInfo(const std::string& message) {
    // Write to file
    if (logFile.is_open()) {
        logFile << "[INFO] " << message << endl;
        logFile.flush();
    }
    // Print important info to console
    std::cout << "[INFO] " << message << std::endl;
}

void logProgress(const std::string& message) {
    // Write to file
    if (logFile.is_open()) {
        logFile << "[PROGRESS] " << message << endl;
        logFile.flush();
    }
    // Print progress to console
    std::cout << message << std::endl;
}

// Forward declaration to avoid E0020 diagnostic in some IDEs
vector<vector<double>> computeDistanceMatrix(const vector<string>& seqs);

unordered_map<char, int> AA_INDEX;

// Global cached full pairwise distance matrix for reuse (e.g., in refinement)
static vector<vector<double>> GLOBAL_DIST_MATRIX;
static bool GLOBAL_DIST_READY = false;

/**
 * Binary tree node used for the guide tree.
 */
struct TreeNode {
    int id;
    int size;
    TreeNode* left;
    TreeNode* right;
    TreeNode(int i) : id(i), size(1), left(nullptr), right(nullptr) {}
    TreeNode(TreeNode* l, TreeNode* r, int newid)
        : id(newid), left(l), right(r), size(l->size + r->size) {}
};

/**
 * Profile representation of an alignment.
 *  - freq: per-column categorical distributions over 21 symbols (20 AA + '-')
 *  - size: number of sequences represented by the profile
 *  - aligned: aligned sequences (strings of equal length)
 */
struct Profile {
    vector<vector<double>> freq;
    int size;
    vector<string> aligned;
    vector<string> names; // aligned sequence names in the same row order
};

/**
 * Reorder final alignment rows to match the original FASTA order.
 * If a name from the original list is not found, it is skipped.
 * If the reordering cannot cover all rows, returns the original order.
 */
pair<vector<string>, vector<string>> reorderAlignmentToOriginal(
    const vector<string>& aligned,
    const vector<string>& currentNames,
    const vector<string>& originalNames)
{
    try {
        if (aligned.size() != currentNames.size()) {
            logWarning("Cannot reorder: names/aligned size mismatch");
            return {aligned, currentNames};
        }

        unordered_map<string, vector<int>> nameToIndices;
        nameToIndices.reserve(currentNames.size());
        for (int i = 0; i < static_cast<int>(currentNames.size()); ++i) {
            nameToIndices[currentNames[i]].push_back(i);
        }

        unordered_map<string, size_t> usedCount;
        vector<string> orderedAligned;
        vector<string> orderedNames;
        orderedAligned.reserve(aligned.size());
        orderedNames.reserve(aligned.size());

        for (const string& nm : originalNames) {
            auto it = nameToIndices.find(nm);
            size_t used = usedCount[nm];
            if (it != nameToIndices.end() && used < it->second.size()) {
                int idx = it->second[used];
                usedCount[nm] = used + 1;
                orderedAligned.push_back(aligned[idx]);
                orderedNames.push_back(currentNames[idx]);
            } else {
                logWarning(string("Original name not found in final alignment: ") + nm);
            }
        }

        if (orderedAligned.size() != aligned.size()) {
            logWarning("Reordering incomplete; preserving final alignment order");
            return {aligned, currentNames};
        }

        return {orderedAligned, orderedNames};
    }
    catch (const exception& e) {
        logError(string("Error during reordering: ") + e.what());
        return {aligned, currentNames};
    }
}

/**
 * Pick the consensus symbol with the maximum weight in a column distribution.
 * @param profile Probability distribution of size 21 over the AA alphabet (including '-')
 * @return The consensus character
 */
char pickConsensus(const vector<double>& profile) {
    double max_weight = -1;
    char selected = '-';
    for (size_t i = 0; i < profile.size(); ++i)
        // Track the most frequent symbol in the column
        if (profile[i] > max_weight) {
            max_weight = profile[i];
            selected = AA[i];
        }
    return selected;
}

/**
 * Pick a consensus symbol if its probability exceeds a threshold; otherwise fall back to max.
 * @param profile Probability distribution of size 21 over the AA alphabet (including '-')
 * @param threshold Minimum fraction required to accept a strong consensus
 * @return The strong consensus character (or max if no symbol passes the threshold)
 */
char pickStrongConsensus(const vector<double>& profile, double threshold = 0.7) {
    double total = accumulate(profile.begin(), profile.end(), 0.0);
    for (size_t i = 0; i < profile.size(); ++i)
        // Accept a strong consensus if it exceeds the threshold
        if (profile[i] / total >= threshold)
            return AA[i];
    return pickConsensus(profile);
}

/**
 * Initialize global mapping from amino acid symbol to index in the 21-letter alphabet.
 */
void initAAIndex() {
    for (int i = 0; i < AA.size(); ++i)
        AA_INDEX[AA[i]] = i;
}

/**
 * Compute Shannon entropy of a categorical distribution.
 * @param p Probability distribution (elements should be non-negative and sum to 1)
 * @return Entropy in bits
 */
double calcEntropy(const vector<double>& p) {
    try {
        double H = 0;
        int non_zero_count = 0;
        for (double val : p) {
            if (val > 0) {
                H -= val * log2(val);
                non_zero_count++;
            }
        }
        // Only log unusually high entropy (very diverse columns)
        if (H > 4.2) {
            logDebug("High entropy detected: " + to_string(H) + " with " + to_string(non_zero_count) + " different amino acids");
        }
        return H;
    }
    catch (const exception& e) {
        logError("Error in calcEntropy: " + string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error in calcEntropy");
        throw;
    }
}

/**
 * Pretty-print the guide tree in ASCII, optionally showing merge distances for internal nodes.
 * @param node Root of the (sub)tree
 * @param names Sequence names for leaf IDs
 * @param nodeDistances Optional map from internal node ID to its merge distance
 * @param prefix Internal: indentation prefix for recursive printing
 * @param isLeft Internal: whether current branch is a left child
 */
void printAsciiTree(const TreeNode* node,
    const std::vector<std::string>& names,
    const std::unordered_map<int, double>& nodeDistances,
    const std::string& prefix = "", bool isLeft = true)
{
if (!node) return;

std::cout << prefix;
std::cout << (isLeft ? "|--" : "`--");

// If leaf, print the sequence name
if (!node->left && !node->right && node->id < names.size()) {
std::cout << names[node->id] << std::endl;
} else {
std::cout << "Node_" << node->id;
if (nodeDistances.count(node->id)) {
std::cout << " [dist = " << nodeDistances.at(node->id) << "]";
}
std::cout << std::endl;
}

if (node->left || node->right) {
printAsciiTree(node->left, names, nodeDistances, prefix + (isLeft ? "|   " : "    "), true);
printAsciiTree(node->right, names, nodeDistances, prefix + (isLeft ? "|   " : "    "), false);
}
}




/**
 * Compute adaptive gap open/extend penalties per column based on entropy.
 * Columns with lower entropy (more conserved) receive larger gap penalties.
 * @param profile Profile distributions (L x 21)
 * @param go_vec [out] Per-column gap-open penalties (size L)
 * @param ge_vec [out] Per-column gap-extend penalties (size L)
 * @param go_max Maximum gap-open penalty
 * @param go_min Minimum gap-open penalty
 * @param ge_max Maximum gap-extend penalty
 * @param ge_min Minimum gap-extend penalty
 */
void computeGapVectors(const vector<vector<double>>& profile, vector<double>& go_vec, vector<double>& ge_vec,
                       double go_max = 10.0, double go_min = 5.0,
                       double ge_max = 5.0, double ge_min = 2.5) {
    int L = profile.size();
    go_vec.resize(L);
    ge_vec.resize(L);
    vector<double> ent(L);
    for (int i = 0; i < L; ++i)
        // Per-column entropy as conservation proxy
        ent[i] = calcEntropy(profile[i]);

    double Hmin = *min_element(ent.begin(), ent.end());
    double Hmax = *max_element(ent.begin(), ent.end());

    for (int i = 0; i < L; ++i) {
        double norm = (Hmax > Hmin) ? (ent[i] - Hmin) / (Hmax - Hmin + 1e-6) : 0.5; // normalize to [0,1]
        norm = std::min(1.0, std::max(0.0, norm));
        // Lower entropy => higher penalties
        go_vec[i] = go_max - norm * (go_max - go_min);
        ge_vec[i] = ge_max - norm * (ge_max - ge_min);
    }
}

/**
 * Convert aligned sequences into per-column categorical distributions over the AA alphabet.
 * @param aligned Vector of aligned sequences (all strings must have the same length)
 * @return Matrix of size L x 21 with per-column frequencies normalized to 1
 * @throws std::invalid_argument When input is empty; may rethrow underlying errors
 */
vector<vector<double>> computeDistributions(const vector<string>& aligned) {
    logDebug("Computing distributions for " + to_string(aligned.size()) + " aligned sequences");
    try {
        if (aligned.empty()) {
            logError("Empty aligned sequences provided to computeDistributions");
            throw invalid_argument("Empty aligned sequences");
        }
        
        int L = aligned[0].size();
        logDebug("Sequence length: " + to_string(L));
        
        vector<vector<double>> dist(L, vector<double>(21, 0.0));
        vector<int> pureCols; // collect pure columns to summarize
        for (int i = 0; i < L; ++i) {
            unordered_map<char, int> count;
            for (size_t s = 0; s < aligned.size(); ++s)
                count[aligned[s][i]]++; // include gaps
            
            double total = 0;
            for (auto& kv : count)
                total += kv.second;
                
            // Log column composition for debugging
            if (count.size() == 1) {
                pureCols.push_back(i);
            } else if (count.size() > 10 && i < 5) { // High diversity columns
                logDebug("Column " + to_string(i) + " has high diversity: " + to_string(count.size()) + " different amino acids");
            }
            
            for (int j = 0; j < 21; ++j)
                dist[i][j] = count.count(AA[j]) ? count[AA[j]] / total : 0.0; // normalize
        }
        if (!pureCols.empty()) {
            // summarize pure columns as requested
            string idxs;
            const int maxShow = 50; // avoid excessively long lines
            for (size_t k = 0; k < pureCols.size(); ++k) {
                if (k) idxs += ",";
                idxs += to_string(pureCols[k]);
                if (k + 1 == maxShow && pureCols.size() > maxShow) { idxs += ",..."; break; }
            }
            logDebug("Columns [" + idxs + "] are pure (entropy will be 0)");
        }
        logDebug("Successfully computed distributions");
        return dist;
    }
    catch (const exception& e) {
        logError("Error in computeDistributions: " + string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error in computeDistributions");
        throw;
    }
}

/**
 * Compute a symmetric pairwise distance matrix between raw sequences.
 * Each sequence is represented as per-position distributions (single sequence),
 * and the distance is the mean EMD across corresponding positions.
 * @param seqs Input sequences (unaligned)
 * @return NxN distance matrix where N = seqs.size()
 * @throws std::invalid_argument When input is empty; may rethrow processing errors
 */
vector<vector<double>> computeDistanceMatrix(const vector<string>& seqs) {
    logInfo("Computing distance matrix for " + to_string(seqs.size()) + " sequences");
    try {
        if (seqs.empty()) {
            logError("Empty sequences provided to computeDistanceMatrix");
            throw invalid_argument("Empty sequences");
        }
        
        int N = seqs.size();
        vector<vector<double>> D(N, vector<double>(N, 0.0));
        
        // Precompute one-hot EMD lookup: tanh(max(cost[i][j], -3.5)) for AA indices
        static vector<vector<double>> emd_lookup;
        if (emd_lookup.empty()) {
            auto cost = getRealBlosumLogOddsMatrix();
            emd_lookup.assign(cost.size(), vector<double>(cost[0].size(), 0.0));
            for (size_t i = 0; i < cost.size(); ++i) {
                for (size_t j = 0; j < cost[i].size(); ++j) {
                    double c = cost[i][j];
                    if (std::isnan(c) || std::isinf(c)) c = 0.0;
                    c = std::max(c, -3.5);
                    emd_lookup[i][j] = std::tanh(c);
                }
            }
        }
        
        logDebug("Computing pairwise EMD distances (one-hot lookup)");
        // Build list of unique (i,j) pairs with i < j to parallelize safely
        vector<pair<int,int>> pairs;
        pairs.reserve((static_cast<long long>(N) * (N - 1)) / 2);
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                pairs.emplace_back(i, j);
            }
        }

        #pragma omp parallel for schedule(static)
        for (int p = 0; p < static_cast<int>(pairs.size()); ++p) {
            int i = pairs[p].first;
            int j = pairs[p].second;
            try {
                double sum_emd = 0.0;
                int L = static_cast<int>(min(seqs[i].size(), seqs[j].size()));
                for (int k = 0; k < L; ++k) {
                    char a = seqs[i][k];
                    char b = seqs[j][k];
                    auto ita = AA_INDEX.find(a);
                    auto itb = AA_INDEX.find(b);
                    if (ita == AA_INDEX.end() || itb == AA_INDEX.end()) {
                        sum_emd += EMD_ERROR; // identical behavior when distribution is invalid
                    } else {
                        sum_emd += emd_lookup[ita->second][itb->second];
                    }
                }
                double value = L ? (sum_emd / L) : EMD_ERROR;
                D[i][j] = value;
                D[j][i] = value;

                if (((static_cast<long long>(i) * N) + j) % 10 == 0) {
                    logDebug("Computed distance for pair (" + to_string(i) + "," + to_string(j) + "): " + to_string(value));
                }
            }
            catch (const exception& e) {
                logError("Error computing EMD for sequences " + to_string(i) + " and " + to_string(j) + ": " + string(e.what()));
                D[i][j] = D[j][i] = EMD_ERROR; // fallback distance
                logDebug("Fallback used: set D[" + to_string(i) + "][" + to_string(j) + "] to EMD_ERROR");
            }
        }
        
        logInfo("Successfully computed distance matrix");
        writeCSV(D, "emd_distance_matrix.csv"); // writes as edge list CSV (i<j)
        logInfo("Distance matrix saved to emd_distance_matrix.csv");
        return D;
    }
    catch (const exception& e) {
        logError("Error in computeDistanceMatrix: " + string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error in computeDistanceMatrix");
        throw;
    }
}

/**
 * Build a UPGMA-like guide tree by iteratively merging the closest clusters.
 * @param D Symmetric NxN distance matrix
 * @param nodeDistances Optional output: record merge distance for internal nodes (keyed by new node ID)
 * @return Root node of the resulting binary tree
 * @throws std::invalid_argument When matrix is empty; std::runtime_error if a minimum pair cannot be found
 */
TreeNode* buildGuideTree(const std::vector<std::vector<double>>& D,
    std::unordered_map<int, double>* nodeDistances = nullptr)
{
    logInfo("Building guide tree from distance matrix");
    try {
        int N = D.size();
        if (N == 0) {
            logError("Empty distance matrix provided to buildGuideTree");
            throw invalid_argument("Empty distance matrix");
        }
        
        logDebug("Distance matrix size: " + to_string(N) + "x" + to_string(N));
        
        std::vector<TreeNode*> nodes;
        std::vector<int> ids(N);
        for (int i = 0; i < N; ++i) {
            nodes.push_back(new TreeNode(i));
            ids[i] = i;
        }

        int nextId = N;
        std::vector<std::vector<double>> dist = D;

        int iteration = 0;
        while (nodes.size() > 1) {
            iteration++;
            logDebug("Guide tree iteration " + to_string(iteration) + ", remaining nodes: " + to_string(nodes.size()));
            
            int a = -1, b = -1;
            double minDist = 1e9;
            
            for (int i = 0; i < nodes.size(); ++i) {
                for (int j = i + 1; j < nodes.size(); ++j) {
                    if (dist[ids[i]][ids[j]] < minDist) {
                        minDist = dist[ids[i]][ids[j]];
                        a = i; 
                        b = j;
                    }
                }
            }
            if (a == -1 || b == -1) {
                logError("Could not find minimum distance pair in guide tree construction");
                throw runtime_error("Unable to find minimum distance pair");
            }
            logDebug("Merging nodes " + to_string(a) + " and " + to_string(b) + " with distance " + to_string(minDist));

            TreeNode* merged = new TreeNode(nodes[a], nodes[b], nextId);

            if (nodeDistances) {
                (*nodeDistances)[nextId] = minDist;
            }

            int sa = nodes[a]->size, sb = nodes[b]->size;
            std::vector<double> newDist(N + iteration + 1, 0.0);
            
            for (int k = 0; k < nodes.size(); ++k) {
                if (k == a || k == b) continue;
                double wa = dist[ids[a]][ids[k]];
                double wb = dist[ids[b]][ids[k]];
                newDist[ids[k]] = (wa * sa + wb * sb) / (sa + sb);
            }

            if (a > b) std::swap(a, b);
            nodes.erase(nodes.begin() + b);
            nodes.erase(nodes.begin() + a);
            ids.erase(ids.begin() + b);
            ids.erase(ids.begin() + a);

            nodes.push_back(merged);
            ids.push_back(dist.size());

            for (auto& row : dist) row.push_back(0.0);
            dist.push_back(newDist);

            nextId++;
        }

        logInfo("Successfully built guide tree");
        return nodes[0];
    }
    catch (const exception& e) {
        logError("Error in buildGuideTree: " + string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error in buildGuideTree");
        throw;
    }
}


/**
 * Align two profiles using affine gap penalties with per-column gap vectors.
 * Uses a BLOSUM-derived EMD dissimilarity plus a small consensus match bonus.
 * @param A Profile A distributions (m x 21)
 * @param B Profile B distributions (n x 21)
 * @param go_vec_A Gap-open penalties per column of A (size m)
 * @param ge_vec_A Gap-extend penalties per column of A (size m)
 * @param go_vec_B Gap-open penalties per column of B (size n)
 * @param ge_vec_B Gap-extend penalties per column of B (size n)
 * @return Pair of (alignment score, list of aligned consensus pairs)
 * @throws std::invalid_argument When any profile is empty; may rethrow underlying errors
 */
pair<double, vector<pair<char, char>>> align(const vector<vector<double>>& A, const vector<vector<double>>& B,
                                             const vector<double>& go_vec_A, const vector<double>& ge_vec_A,
                                             const vector<double>& go_vec_B, const vector<double>& ge_vec_B) {
    logDebug("Starting profile alignment - A size: " + to_string(A.size()) + ", B size: " + to_string(B.size()));
    
    try {
        if (A.empty() || B.empty()) {
            logError("Empty profile provided to align function");
            throw invalid_argument("Empty profile in align");
        }
        
        int m = A.size(), n = B.size();
        vector<vector<double>> M(m + 1, vector<double>(n + 1, -1e9)); // match/mismatch
        vector<vector<double>> X = M, Y = M; // affine gap states
        M[0][0] = 0;

        logDebug("Initializing gap penalties");
        for (int i = 1; i <= m; ++i)
            X[i][0] = -go_vec_A[i - 1] - (i - 1) * ge_vec_A[i - 1];
        for (int j = 1; j <= n; ++j)
            Y[0][j] = -go_vec_B[j - 1] - (j - 1) * ge_vec_B[j - 1];

        logDebug("Computing dynamic programming matrix");
        for (int k = 2; k <= m + n; ++k) {
            #pragma omp parallel for schedule(dynamic)
            for (int i = max(1, k - n); i <= min(k - 1, m); ++i) {
                int j = k - i;

                double go = (go_vec_A[i - 1] + go_vec_B[j - 1]) / 2.0;
                double ge = (ge_vec_A[i - 1] + ge_vec_B[j - 1]) / 2.0;

                double emd_raw = calculateEMD(A[i - 1], B[j - 1]);
                char a_con = pickStrongConsensus(A[i - 1]);
                char b_con = pickStrongConsensus(B[j - 1]);
                double consensus_bonus = (a_con == b_con) ? -3.0 : 0.0;

                double emd = -emd_raw + consensus_bonus;

                M[i][j] = emd + max({M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]});

                double entropy_i = calcEntropy(A[i - 1]);
                double entropy_j = calcEntropy(B[j - 1]);
                double entropy_factor = 1.0 + 0.75 * (1.0 - std::min(entropy_i, entropy_j));

                double adj_go = go * entropy_factor;
                double adj_ge = ge * entropy_factor;

                X[i][j] = max(M[i - 1][j] - adj_go, X[i - 1][j] - adj_ge);
                Y[i][j] = max(M[i][j - 1] - adj_go, Y[i][j - 1] - adj_ge);
            }    
        }

        logDebug("Performing traceback");
        vector<pair<char, char>> aln;
        int i = m, j = n;
        while (i > 0 || j > 0) { 
            if (i > 0 && j > 0 && M[i][j] >= X[i][j] && M[i][j] >= Y[i][j]) {
                aln.push_back({pickStrongConsensus(A[i - 1]), pickStrongConsensus(B[j - 1])});
                --i; --j;
            } else if (i > 0 && X[i][j] >= Y[i][j]) {
                aln.push_back({pickStrongConsensus(A[i - 1]), '-'});
                --i;
            } else if (j > 0){
                aln.push_back({'-', pickStrongConsensus(B[j - 1])});
                --j;
            } else {
                break;
            }    
        }

        reverse(aln.begin(), aln.end());
        double final_score = M[m][n];
        
        logDebug("Alignment completed - Score: " + to_string(final_score) + ", Length: " + to_string(aln.size()));
        return {final_score, aln};
    }
    catch (const exception& e) {
        logError("Error in align function: " + string(e.what()));
        throw;
    }
    catch (...) {
        logError("Unknown error in align function");
        throw;
    }
}

/**
 * Combine two profiles according to a column-wise alignment mapping.
 * @param A Left profile
 * @param B Right profile
 * @param aln Aligned consensus character pairs produced by align()
 * @return Merged profile including updated distributions and aligned sequences
 */
Profile combine(const Profile& A, const Profile& B, const vector<pair<char, char>>& aln) {
    int L = aln.size();
    vector<vector<double>> prof(L, vector<double>(21, 0.0));
    vector<string> aligned(A.aligned.size() + B.aligned.size(), string(L, '-'));
    vector<string> names;
    names.reserve(A.aligned.size() + B.aligned.size());
    // Row order: first all rows from A, then all rows from B
    for (const auto& n : A.names) names.push_back(n);
    for (const auto& n : B.names) names.push_back(n);
    int ai = 0, bi = 0;

    for (int i = 0; i < L; ++i) {
        char a = aln[i].first, b = aln[i].second;
        if (a != '-') for (int j = 0; j < 21; ++j) prof[i][j] += A.freq[ai][j] * A.size;
        if (b != '-') for (int j = 0; j < 21; ++j) prof[i][j] += B.freq[bi][j] * B.size;
        double total = (a != '-' ? A.size : 0) + (b != '-' ? B.size : 0);
        if (total > 0) for (int j = 0; j < 21; ++j) prof[i][j] /= total;

        for (int k = 0; k < A.aligned.size(); ++k)
            aligned[k][i] = (a != '-') ? A.aligned[k][ai] : '-';
        for (int k = 0; k < B.aligned.size(); ++k)
            aligned[k + A.aligned.size()][i] = (b != '-') ? B.aligned[k][bi] : '-';

        if (a != '-') ai++;
        if (b != '-') bi++;
    }

    return {prof, A.size + B.size, aligned, names};
}

/**
 * Recursively align profiles along a guide tree (progressive alignment).
 * @param node Root of the current subtree
 * @param profiles Leaf profiles for original sequences
 * @return Profile representing the alignment of sequences in this subtree
 * @throws std::out_of_range If node IDs are invalid; may rethrow underlying errors
 */
Profile alignFromTree(TreeNode* node, const vector<Profile>& profiles) {
    if (!node->left && !node->right)
        return profiles[node->id];

    Profile left = alignFromTree(node->left, profiles);
    Profile right = alignFromTree(node->right, profiles);

    vector<double> goL, geL, goR, geR;
    computeGapVectors(left.freq, goL, geL, 10.0, 5.0, 5.0, 2.5);
    computeGapVectors(right.freq, goR, geR, 10.0, 5.0, 5.0, 2.5);

    auto aln = align(left.freq, right.freq, goL, geL, goR, geR);
    return combine(left, right, aln.second);
}

/**
 * Compute a simple column-similarity score.
 * A column is counted if all residues are reasonably similar to the first
 * according to a BLOSUM-derived threshold.
 * @param aligned Final aligned sequences
 * @return Fraction of columns deemed similar
 */
double computeCSScore(const vector<string>& aligned) {
    int n = aligned.size();
    if (n == 0) return 0.0;
    int L = aligned[0].size();
    if (L == 0) return 0.0;

    auto blosum = getRealBlosumLogOddsMatrix();
    int match_count = 0;

    for (int i = 0; i < L; ++i) {
        char ref = aligned[0][i];
        if (ref == '-' || AA_INDEX.find(ref) == AA_INDEX.end()) continue;
        int ref_idx = AA_INDEX[ref];

        bool similar = true;
        for (int j = 1; j < n; ++j) {
            char cj = aligned[j][i];
            if (cj == '-' || AA_INDEX.find(cj) == AA_INDEX.end()) {
                similar = false;
                break;
            }
            int cj_idx = AA_INDEX[cj];
            double blosum_score = blosum[ref_idx][cj_idx];
            if (blosum_score < -3.0) {
                similar = false;
                break;
            }
        }

        if (similar) match_count++;
    }

    return static_cast<double>(match_count) / L;
}


/**
 * Compute the Sum-of-Pairs Score (SPS) variant using BLOSUM or strict identity.
 * @param aligned Aligned sequences
 * @param threshold Minimum BLOSUM log-odds to consider a pair matching
 * @param use_blosum If true, use BLOSUM-based similarity; otherwise use exact identity
 * @param debug If true, prints matched vs total pair counts
 * @return Fraction of matching residue pairs across all columns
 */
double computeSPSScore(const vector<string>& aligned, double threshold = -3.0, bool use_blosum = true, bool debug = false) {
    int n = aligned.size();
    if (n < 2) return 0.0;

    int L = aligned[0].size();
    if (L == 0) return 0.0;

    auto blosum = getRealBlosumLogOddsMatrix();
    int total_pairs = 0;
    int matching_pairs = 0;

    for (int col = 0; col < L; ++col) {
        for (int i = 0; i < n; ++i) {
            char a = aligned[i][col];
            if (a == '-' || AA_INDEX.find(a) == AA_INDEX.end()) continue;
            int a_idx = AA_INDEX[a];

            for (int j = i + 1; j < n; ++j) {
                char b = aligned[j][col];
                if (b == '-' || AA_INDEX.find(b) == AA_INDEX.end()) continue;
                int b_idx = AA_INDEX[b];

                total_pairs++;

                if (use_blosum) {
                    if (blosum[a_idx][b_idx] >= threshold) {
                        matching_pairs++;
                    }
                } else {
                    if (a == b) {
                        matching_pairs++;
                    }
                }
            }
        }
    }

    if (total_pairs == 0) return 0.0;
    if (debug) {
        cout << "Matched pairs: " << matching_pairs << ", Total pairs: " << total_pairs << endl;
    }

    return static_cast<double>(matching_pairs) / total_pairs;
}

/**
 * Identify columns where the dominant residue accounts for all but a small number of mismatches.
 * @param aligned Aligned sequences
 * @param mismatch_allow Maximum number of non-dominant non-gap residues allowed
 * @return Indices of strong columns
 */
vector<int> findStrongColumns(const vector<string>& aligned, int mismatch_allow) {
    int n = aligned.size(), L = aligned[0].size();
    vector<int> strongCols;

    for (int i = 0; i < L; ++i) {
        unordered_map<char, int> freq;
        for (int j = 0; j < n; ++j)
            if (aligned[j][i] != '-') freq[aligned[j][i]]++;

        int max_count = 0;
        for (const auto& kv : freq)
            max_count = max(max_count, kv.second);

        if (n - max_count <= mismatch_allow)
            strongCols.push_back(i);
    }
    return strongCols;
}

/**
 * Perform simple iterative refinement by re-adding each sequence to the alignment
 * built from the remaining sequences, keeping improvements and anchor consistency.
 * @param seqs Input sequences
 * @param initialProfile Starting MSA profile
 * @param iterations Number of refinement passes
 * @param verbose If true, prints progress details
 * @return Refined profile
 * @throws std::exception May rethrow underlying errors during refinement steps
 */
Profile iterativeRefinement(const vector<string>& seqs,
    const vector<string>& names,
    const Profile& initialProfile,
    int iterations = 3,
    bool verbose = false)
{
Profile best = initialProfile;
double bestCS = computeCSScore(best.aligned);
if (verbose) cout << "\n Initial CS Score: " << bestCS << endl;

vector<int> anchorColumns = findStrongColumns(best.aligned, 1);

for (int iter = 0; iter < iterations; ++iter) {
if (verbose) cout << "\n Iteration " << iter + 1 << endl;

for (int i = 0; i < seqs.size(); ++i) {
vector<string> rest;
vector<string> restNames;
vector<int> restIdx;
for (int j = 0; j < seqs.size(); ++j) {
if (i != j) { rest.push_back(seqs[j]); restNames.push_back(names[j]); restIdx.push_back(j); }
}

vector<Profile> partialProfiles;
for (int ri = 0; ri < rest.size(); ++ri) {
    const auto& s = rest[ri];
    Profile p;
    p.freq = computeDistributions({s});
    p.size = 1;
    p.aligned = {s};
    p.names = { restNames[ri] };
    partialProfiles.push_back(std::move(p));
}

vector<vector<double>> dist;
if (GLOBAL_DIST_READY) {
    int R = static_cast<int>(restIdx.size());
    dist.assign(R, vector<double>(R, 0.0));
    for (int a = 0; a < R; ++a) {
        for (int b = 0; b < R; ++b) {
            dist[a][b] = GLOBAL_DIST_MATRIX[restIdx[a]][restIdx[b]];
        }
    }
    try {
        writeCSV(dist, "emd_distance_matrix.csv"); // preserve previous side-effect
    } catch (...) {
        // Ignore write failures here to match robustness of original pipeline
    }
} else {
    dist = computeDistanceMatrix(rest);
}
TreeNode* tree = buildGuideTree(dist);
Profile alignedRest = alignFromTree(tree, partialProfiles);

Profile toAdd;
toAdd.freq = computeDistributions({seqs[i]});
toAdd.size = 1;
toAdd.aligned = {seqs[i]};
toAdd.names = { names[i] };
vector<double> go1, ge1, go2, ge2;
computeGapVectors(alignedRest.freq, go1, ge1, 10.0, 5.0, 5.0, 2.5);
computeGapVectors(toAdd.freq, go2, ge2, 10.0, 5.0, 5.0, 2.5);

auto aln = align(alignedRest.freq, toAdd.freq, go1, ge1, go2, ge2);
Profile combined = combine(alignedRest, toAdd, aln.second);

double cs = computeCSScore(combined.aligned);
vector<int> newStrong = findStrongColumns(combined.aligned, 1);

int anchor_overlap = 0;
for (int col : newStrong)
if (find(anchorColumns.begin(), anchorColumns.end(), col) != anchorColumns.end())
anchor_overlap++;

double anchor_fraction = static_cast<double>(anchor_overlap) / std::max(anchorColumns.size(), size_t(1));
if (verbose) {
cout << "  ➤ Anchor match: " << anchor_fraction << endl;
cout << "  ➤ New CS after adding seq " << i << ": " << cs << endl;
}

if (cs > bestCS || anchor_fraction > 0.7) {
best = combined;
bestCS = cs;
if (verbose) cout << "  CS updated!\n";
}
}
}

return best;
}


int main() {
    // Initialize log file first
    initializeLogFile();
    logInfo("Starting Multiple Sequence Alignment (MSA) program");
    logProgress("Initializing MSA program...");
    auto programStart = std::chrono::steady_clock::now();
    
    try {
        srand(static_cast<unsigned int>(time(0)));
        logDebug("Random seed initialized");
        
        initAAIndex();
        logDebug("Amino acid index initialized");

        logProgress("Reading FASTA file: input.fasta");
        auto [seqs, names] = readFastaWithNames("input.fasta");

        if (seqs.empty()) {
            logError("FASTA file empty or not found.");
            cerr << "FASTA file empty or not found.\n";
            closeLogFile();
            return 1;
        }
        
        logInfo("Successfully read " + to_string(seqs.size()) + " sequences");
        logProgress("Loaded " + to_string(seqs.size()) + " protein sequences");
        for (size_t i = 0; i < seqs.size() && i < 5; ++i) {
            logDebug("Sequence " + to_string(i+1) + " (" + names[i] + "): length " + to_string(seqs[i].length()));
        }
        if (seqs.size() > 5) {
            logDebug("... and " + to_string(seqs.size() - 5) + " more sequences");
        }

        logProgress("Computing pairwise distances...");
        vector<vector<double>> D;
        try {
            D = computeDistanceMatrix(seqs);
            GLOBAL_DIST_MATRIX = D;
            GLOBAL_DIST_READY = true;
        }
        catch (const exception& e) {
            logError("Failed to compute distance matrix: " + string(e.what()));
            closeLogFile();
            return 1;
        }

        logProgress("Creating sequence profiles...");
        vector<Profile> profiles;
        try {
            for (size_t idx = 0; idx < seqs.size(); ++idx) {
                const auto& s = seqs[idx];
                const auto& nm = names[idx];
                Profile p;
                p.freq = computeDistributions({s});
                p.size = 1;
                p.aligned = {s};
                p.names = {nm};
                profiles.push_back(std::move(p));
            }
            logDebug("Created " + to_string(profiles.size()) + " initial profiles");
        }
        catch (const exception& e) {
            logError("Failed to create initial profiles: " + string(e.what()));
            closeLogFile();
            return 1;
        }

        logProgress("Building guide tree...");
        std::unordered_map<int, double> nodeDistances;
        TreeNode* tree = nullptr;
        try {
            tree = buildGuideTree(D, &nodeDistances);
            logInfo("Guide tree construction completed");
        }
        catch (const exception& e) {
            logError("Failed to build guide tree: " + string(e.what()));
            closeLogFile();
            return 1;
        }
        
        cout << "\nGuide Tree:" << endl;
        try {
            printAsciiTree(tree, names, nodeDistances);
        }
        catch (const exception& e) {
            logError("Error printing guide tree: " + string(e.what()));
        }

        logProgress("Performing progressive alignment...");
        Profile initial;
        try {
            initial = alignFromTree(tree, profiles);
            logInfo("Progressive alignment completed");
        }
        catch (const exception& e) {
            logError("Failed in progressive alignment: " + string(e.what()));
            closeLogFile();
            return 1;
        }

        logProgress("Starting iterative refinement...");
        Profile final;
        try {
            final = iterativeRefinement(seqs, names, initial, 3, true);  
            logInfo("Iterative refinement completed");
        }
        catch (const exception& e) {
            logError("Failed in iterative refinement: " + string(e.what()));
            // Continue with initial alignment if refinement fails
            logInfo("Using initial alignment as final result");
            logDebug("Fallback used: iterative refinement failed, proceeding with initial alignment");
            final = initial;
        }

        // Reorder final alignment to original FASTA order so that A3M master is the input's first sequence
        auto reordered = reorderAlignmentToOriginal(final.aligned, final.names, names);
        const vector<string>& outAligned = reordered.first;
        const vector<string>& outNames = reordered.second;
        logInfo("Final alignment reordered to original FASTA order for output");

        cout << "\n=== FINAL MULTIPLE ALIGNMENT ===" << endl;
        if (outNames.size() == outAligned.size()) {
            for (size_t i = 0; i < outAligned.size(); ++i) {
                cout << ">" << outNames[i] << "\n" << outAligned[i] << endl;
            }
        } else {
            for (size_t i = 0; i < outAligned.size(); ++i) {
                cout << ">Sequence_" << i + 1 << "\n" << outAligned[i] << endl;
            }
        }

        logProgress("Writing results to file...");
        try {
            if (outNames.size() == outAligned.size()) {
                writeFasta(outAligned, outNames, "final_msa.fasta");
                writeA3M(outAligned, outNames, "final_msa.a3m");
            } else {
                writeFasta(outAligned, "final_msa.fasta");
                writeA3M(outAligned, "final_msa.a3m");
            }
            cout << "\nOutput written to final_msa.fasta" << endl;
            cout << "A3M written to final_msa.a3m" << endl;
            logInfo("Output file written successfully");
        }
        catch (const exception& e) {
            logError("Failed to write output file: " + string(e.what()));
            closeLogFile();
            return 1;
        }

        // Compute and log final alignment statistics
        try {
            double cs_score = computeCSScore(final.aligned);
            double sps_score = computeSPSScore(final.aligned, -3.0, true, false); // Disable verbose output
            logInfo("Final alignment statistics:");
            logInfo("- Column similarity score: " + to_string(cs_score));
            logInfo("- Sum-of-pairs score: " + to_string(sps_score));
            logInfo("- Number of sequences: " + to_string(final.aligned.size()));
            if (!final.aligned.empty()) {
                logInfo("- Alignment length: " + to_string(final.aligned[0].length()));
            }
            
            cout << "\nAlignment Statistics:" << endl;
            cout << "- Sequences: " << final.aligned.size() << endl;
            if (!final.aligned.empty()) {
                cout << "- Length: " << final.aligned[0].length() << endl;
            }
            cout << "- Column similarity: " << fixed << setprecision(3) << cs_score << endl;
            cout << "- Sum-of-pairs score: " << fixed << setprecision(3) << sps_score << endl;
        }
        catch (const exception& e) {
            logError("Error computing alignment statistics: " + string(e.what()));
        }

        logProgress("MSA program completed successfully!");
        auto programEnd = std::chrono::steady_clock::now();
        double elapsedSec = std::chrono::duration<double>(programEnd - programStart).count();
        logInfo("Total runtime (s): " + to_string(elapsedSec));
        cout << "Total runtime: " << fixed << setprecision(3) << elapsedSec << " s" << endl;
        logInfo("Check 'msa_debug.log' for detailed debug information");
        closeLogFile();
        return 0;
    }
    catch (const exception& e) {
        logError("Critical error in main: " + string(e.what()));
        cerr << "Critical error: " << e.what() << endl;
        closeLogFile();
        return 1;
    }
    catch (...) {
        logError("Unknown critical error in main");
        cerr << "Unknown critical error occurred" << endl;
        closeLogFile();
        return 1;
    }
}


