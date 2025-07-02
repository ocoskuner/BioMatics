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
#include "blosum62.h"
#include "emd_solver.h"
#include "fasta_reader.h"
#include "output_writer.h"
#include <omp.h>

using namespace std;


vector<vector<double>> computeDistanceMatrix(const vector<string>& seqs, const vector<string>& names);


unordered_map<char, int> AA_INDEX;

struct TreeNode {
    int id;
    int size;
    TreeNode* left;
    TreeNode* right;
    TreeNode(int i) : id(i), size(1), left(nullptr), right(nullptr) {}
    TreeNode(TreeNode* l, TreeNode* r, int newid)
        : id(newid), left(l), right(r), size(l->size + r->size) {}
};

struct Profile {
    vector<vector<double>> freq;
    int size;
    vector<string> aligned;
};

char pickConsensus(const vector<double>& profile) {
    double max_weight = -1;
    char selected = '-';
    for (size_t i = 0; i < profile.size(); ++i)
        if (profile[i] > max_weight) {
            max_weight = profile[i];
            selected = AA[i];
        }
    return selected;
}

char pickStrongConsensus(const vector<double>& profile, double threshold = 0.7) {
    double total = accumulate(profile.begin(), profile.end(), 0.0);
    for (size_t i = 0; i < profile.size(); ++i)
        if (profile[i] / total >= threshold)
            return AA[i];
    return pickConsensus(profile);
}

void initAAIndex() {
    for (int i = 0; i < AA.size(); ++i)
        AA_INDEX[AA[i]] = i;
}

double calcEntropy(const vector<double>& p) {
    double H = 0;
    for (double val : p)
        if (val > 0) H -= val * log2(val);
    return H;
}

void deleteTree(TreeNode* node) {
    if (!node) return;
    deleteTree(node->left);
    deleteTree(node->right);
    delete node;
}

void printAsciiTree(const TreeNode* node,
    const std::vector<std::string>& names,
    const std::unordered_map<int, double>& nodeDistances,
    const std::string& prefix = "", bool isLeft = true)
{
if (!node) return;

std::cout << prefix;
std::cout << (isLeft ? "|--" : "`--");


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




void computeGapVectors(const vector<vector<double>>& profile,
                       vector<double>& go_vec, vector<double>& ge_vec,
                       double go_max = 10.0, double go_min = 5.0,
                       double ge_max = 5.0, double ge_min = 2.5) {
    int L = profile.size();
    go_vec.resize(L);
    ge_vec.resize(L);
    vector<double> ent(L);

    #pragma omp parallel for
    for (int i = 0; i < L; ++i)
        ent[i] = calcEntropy(profile[i]);

    double Hmin = *min_element(ent.begin(), ent.end());
    double Hmax = *max_element(ent.begin(), ent.end());

    double length_factor = 1.0 + log(1 + L / 100.0);

    #pragma omp parallel for
    for (int i = 0; i < L; ++i) {
        double norm = (Hmax > Hmin) ? (ent[i] - Hmin) / (Hmax - Hmin + 1e-6) : 0.5;
        norm = std::min(1.0, std::max(0.0, norm));
        go_vec[i] = (go_max - norm * (go_max - go_min)) * length_factor;
        ge_vec[i] = (ge_max - norm * (ge_max - ge_min)) * length_factor;
    }
}

vector<vector<double>> computeDistributions(const vector<string>& aligned) {
    int L = aligned[0].size();
    vector<vector<double>> dist(L, vector<double>(21, 0.0));
    for (int i = 0; i < L; ++i) {
        unordered_map<char, int> count;
        for (size_t s = 0; s < aligned.size(); ++s)
            count[aligned[s][i]]++;
        double total = 0;
        for (auto& kv : count)
            total += kv.second;
        for (int j = 0; j < 21; ++j)
            dist[i][j] = count.count(AA[j]) ? count[AA[j]] / total : 0.0;
    }
    return dist;
}

#include <omp.h> 

vector<vector<double>> computeDistanceMatrix(const vector<string>& seqs, const vector<string>& names) {
    int N = seqs.size();
    vector<vector<double>> D(N, vector<double>(N, 0.0));
    vector<vector<vector<double>>> distros(N);

    for (int i = 0; i < N; ++i)
        distros[i] = computeDistributions({seqs[i]});

    #pragma omp parallel for collapse(2) shared(D, distros)
    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j) {
            double sum_emd = 0.0;
            int L = std::min(distros[i].size(), distros[j].size());
            for (int k = 0; k < L; ++k)
                sum_emd += calculateEMD(distros[i][k], distros[j][k]);
            D[i][j] = D[j][i] = L ? sum_emd / L : 1e6;
        }

    writeCSV(D, names, "emd_distance_matrix.csv");
    return D;
}



TreeNode* buildGuideTree(const std::vector<std::vector<double>>& D,
    std::unordered_map<int, double>* nodeDistances = nullptr)
{
int N = D.size();
std::vector<TreeNode*> nodes;
std::vector<int> ids(N);
for (int i = 0; i < N; ++i) {
nodes.push_back(new TreeNode(i));
ids[i] = i;
}

int nextId = N;
std::vector<std::vector<double>> dist = D;

while (nodes.size() > 1) {
int a = -1, b = -1;
double minDist = 1e9;
for (int i = 0; i < nodes.size(); ++i)
for (int j = i + 1; j < nodes.size(); ++j)
if (dist[ids[i]][ids[j]] < minDist) {
minDist = dist[ids[i]][ids[j]];
a = i; b = j;
}

TreeNode* merged = new TreeNode(nodes[a], nodes[b], nextId);

if (nodeDistances) {
(*nodeDistances)[nextId] = minDist;
}

int sa = nodes[a]->size, sb = nodes[b]->size;
std::vector<double> newDist(dist.size() + 1, 0.0);
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

return nodes[0];
}


pair<double, vector<pair<char, char>>> align(const vector<vector<double>>& A, const vector<vector<double>>& B,
    const vector<double>& go_vec_A, const vector<double>& ge_vec_A,
    const vector<double>& go_vec_B, const vector<double>& ge_vec_B) {
    int m = A.size(), n = B.size();
    vector<vector<double>> M(m + 1, vector<double>(n + 1, -1e9));
    vector<vector<double>> X = M, Y = M;
    vector<vector<int>> trace(m + 1, vector<int>(n + 1, -1));  // 0: M, 1: X, 2: Y
    M[0][0] = 0;

    for (int i = 1; i <= m; ++i)
        X[i][0] = -go_vec_A[i - 1] - (i - 1) * ge_vec_A[i - 1];
    for (int j = 1; j <= n; ++j)
        Y[0][j] = -go_vec_B[j - 1] - (j - 1) * ge_vec_B[j - 1];

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            double go = (go_vec_A[i - 1] + go_vec_B[j - 1]) / 2.0;
            double ge = (ge_vec_A[i - 1] + ge_vec_B[j - 1]) / 2.0;

            double emd_raw = calculateEMD(A[i - 1], B[j - 1]);
            char a_con = pickStrongConsensus(A[i - 1]);
            char b_con = pickStrongConsensus(B[j - 1]);
            double consensus_bonus = (a_con == b_con) ? -3.0 : 0.0;
            double emd = -emd_raw + consensus_bonus;

            double match_score = emd + max({ M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1] });
            M[i][j] = match_score;
            trace[i][j] = 0;

            double entropy_i = calcEntropy(A[i - 1]);
            double entropy_j = calcEntropy(B[j - 1]);
            double entropy_factor = 1.0 + 0.75 * (1.0 - std::min(entropy_i, entropy_j));
            double adj_go = go * entropy_factor;
            double adj_ge = ge * entropy_factor;

            double x_val = max(M[i - 1][j] - adj_go, X[i - 1][j] - adj_ge);
            double y_val = max(M[i][j - 1] - adj_go, Y[i][j - 1] - adj_ge);

            X[i][j] = x_val;
            Y[i][j] = y_val;

            if (x_val > M[i][j]) {
                M[i][j] = x_val;
                trace[i][j] = 1;
            }
            if (y_val > M[i][j]) {
                M[i][j] = y_val;
                trace[i][j] = 2;
            }
        }
    }

    vector<pair<char, char>> aln;
    int i = m, j = n;

    while (i > 0 || j > 0) {
        if (trace[i][j] == 0 && i > 0 && j > 0) {
            aln.push_back({ 'M', 'M' });  // Match
            --i; --j;
        }
        else if (trace[i][j] == 1 && i > 0) {
            aln.push_back({ 'M', '-' });  // Deletion in B
            --i;
        }
        else if (trace[i][j] == 2 && j > 0) {
            aln.push_back({ '-', 'M' });  // Insertion in B
            --j;
        }
        else {
            break;
        }
    }

    reverse(aln.begin(), aln.end());
    return { M[m][n], aln };
}

Profile combine(const Profile& A, const Profile& B, const vector<pair<char, char>>& aln) {
    int L = aln.size();
    vector<vector<double>> prof(L, vector<double>(21, 0.0));
    vector<string> aligned(A.aligned.size() + B.aligned.size(), string(L, '-'));
    int ai = 0, bi = 0;

    for (int i = 0; i < L; ++i) {
        char a = aln[i].first;
        char b = aln[i].second;

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

    return { prof, A.size + B.size, aligned };
}

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

Profile iterativeRefinement(const vector<string>& seqs,
    const vector<string>& names,
    const Profile& initialProfile,
    int iterations = 3,
    bool verbose = false)

{
    Profile best = initialProfile;
    double bestCS = computeCSScore(best.aligned);
    if (verbose) cout << "\n🔁 Initial CS Score: " << bestCS << endl;

    vector<int> anchorColumns = findStrongColumns(best.aligned, 1);
   
    if ((bestCS < 0.08 && seqs.size() > 5) ||
        (best.aligned[0].size() > 250) ||
        (best.aligned[0].size() > 200 && seqs.size() > 6)) {
        if (verbose) std::cout << "⚠️ Refinement skipped: ";
        if (bestCS < 0.1 && seqs.size() > 5)
            std::cout << "CS too low and seq count high.\n";
        else if (best.aligned[0].size() > 250)
            std::cout << "Alignment too long.\n";
        else
            std::cout << "Long alignment + many sequences.\n";
        return best;
    }

    for (int iter = 0; iter < iterations; ++iter) {
        if (verbose) cout << "\n🌀 Iteration " << iter + 1 << endl;

        for (int i = 0; i < seqs.size(); ++i) {
            vector<string> rest;
            vector<string> names_rest;
            for (int j = 0; j < seqs.size(); ++j) {
                if (i != j) {
                    rest.push_back(seqs[j]);
                    names_rest.push_back(names[j]);
                }
            }

            vector<Profile> partialProfiles;
            for (const auto& s : rest)
                partialProfiles.push_back({ computeDistributions({s}), 1, {s} });

            vector<vector<double>> dist = computeDistanceMatrix(rest, names_rest);
            TreeNode* tree = buildGuideTree(dist);
            Profile alignedRest = alignFromTree(tree, partialProfiles);
            deleteTree(tree);  

            Profile toAdd = { computeDistributions({seqs[i]}), 1, {seqs[i]} };
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
    srand(static_cast<unsigned int>(time(0)));
    initAAIndex();

    auto [seqs, names] = readFastaWithNames("input.fasta");

    if (seqs.empty()) {
        cerr << " FASTA file empty or not found.\n";
        return 1;
    }

    cout << "✅ Input read. Number of sequences: " << seqs.size() << endl;

    vector<vector<double>> D = computeDistanceMatrix(seqs, names);

    cout << "✅ Distance matrix computed.\n";

    vector<Profile> profiles;
    for (const auto& s : seqs)
        profiles.push_back({computeDistributions({s}), 1, {s}});

    std::unordered_map<int, double> nodeDistances;
    TreeNode* tree = buildGuideTree(D, &nodeDistances);
    cout << "✅ Guide tree built.\n";

    Profile initial = alignFromTree(tree, profiles);
    cout << "✅ Initial progressive alignment done.\n";

    Profile final = iterativeRefinement(seqs, names, initial, 3, true);

cout << "?? Final alignment has " << final.aligned.size() << " sequences." << endl;
cout << "?? First aligned sequence length: " << (final.aligned.empty() ? 0 : final.aligned[0].size()) << endl;

    cout << "\n--- FINAL MULTIPLE ALIGNMENT ---\n";
    for (const auto& row : final.aligned)
        cout << row << endl;

    writeFasta(final.aligned, "final_msa.fasta");
    cout << "\n📄 Output written to final_msa.fasta\n";

    return 0;
}

