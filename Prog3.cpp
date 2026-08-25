/*
 Name: KAZI ISLAM
===============================================================================
Program: Median String and Motif (Consensus) Finder
-------------------------------------------------------------------------------
Description:
    This program finds the top 5 median strings of a given motif length (L)
    from a set of DNA sequences provided in a FASTA file. It then computes
    the consensus (motif) strings and their scores.

Inputs:
    1. Input FASTA file path (string)
    2. Motif length L (integer)
    3. Output file path (string)

Outputs:
    - The best 5 median strings with their total distances (tot_dist)
    - Their motif consensus strings with consensus scores
    - Motif positions and substrings from the original sequences

Methods / Algorithms Used:
    - Hamming distance to measure dissimilarity between motifs
    - Brute-force generation of all possible L-mers (4^L combinations)
    - For each L-mer, compute total distance to all sequences
    - Select the top 5 lowest total-distance L-mers (median strings)
    - Build consensus motifs by counting the most frequent base at each position

Compilation:
    g++ -o median_finder Prog3.cpp


Execution Example:
    ./median_finder HMP-part.fa 6 results.txt

Notes:
    - Handles invalid/messy FASTA input by replacing invalid characters with 'N'
    - 'N' in sequences is treated as a mismatch during Hamming distance
    - Supports sequences shorter than motif length (treated with placeholders)
===============================================================================
*/

/**
 * @brief Compute Hamming distance between two DNA strings.
 * @param a First DNA string.
 * @param b Second DNA string.
 * @return Integer distance (number of mismatched bases).
 *
 * Notes:
 *  - 'N' is treated as a mismatch with any character.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <queue>
#include <cctype>
#include <limits>

using namespace std;


static inline int hammingDistance(const string& a, const string& b) {
    int d = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        char x = a[i], y = b[i];
        if (x == 'N' || y == 'N' || x != y) ++d; // 'N' treated as mismatch
    }
    return d;
}

/**
 * @brief Read DNA sequences from a FASTA file.
 * @param filename Path to FASTA input file.
 * @return Vector of sequences as strings.
 *
 * Notes:
 *  - Non-ACGT characters are replaced with 'N'.
 *  - Lines starting with '>' are treated as sequence headers.
 */

vector<string> readFasta(const string& filename) {
    ifstream in(filename);
    vector<string> seqs;
    if (!in) {
        cerr << "ERROR: cannot open input file: " << filename << "\n";
        return seqs;
    }
    string line, cur;
    while (getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            if (!cur.empty()) { seqs.push_back(cur); cur.clear(); }
        } else {
            for (char ch : line) {
                if (isalpha(static_cast<unsigned char>(ch))) {
                    char x = toupper(static_cast<unsigned char>(ch));
                    if (x=='A'||x=='C'||x=='G'||x=='T') cur.push_back(x);
                    else cur.push_back('N');
                }
            }
        }
    }
    if (!cur.empty()) seqs.push_back(cur);
    return seqs;
}

/**
 * @brief Generate all possible L-length DNA strings (L-mers).
 * @param L Motif length.
 * @return Vector containing all combinations of A, C, G, T of length L.
 */

vector<string> generateLmers(int L) {
    vector<string> out;
    if (L <= 0) return out;
    size_t total = 1;
    for (int i = 0; i < L; ++i) total *= 4; // exact integer count
    const char bases[4] = {'A','C','G','T'};
    out.reserve(total);
    for (size_t idx = 0; idx < total; ++idx) {
        size_t x = idx;
        string s(L, 'A');
        for (int p = L - 1; p >= 0; --p) {
            s[p] = bases[x % 4];
            x /= 4;
        }
        out.push_back(std::move(s));
    }
    return out;
}

/**
 * @brief Compute total distance of an L-mer to a set of sequences.
 * @param lmer The L-mer being evaluated.
 * @param seqs Vector of DNA sequences.
 * @param L Length of motif.
 * @return Total integer distance across all sequences.
 */

int totalDistanceToSet(const string& lmer, const vector<string>& seqs, int L) {
    int total = 0;
    for (const auto& seq : seqs) {
        if ((int)seq.size() < L) { total += L; continue; } // edge case
        int best = L;
        for (size_t i = 0; i + L <= seq.size(); ++i) {
            int d = hammingDistance(lmer, seq.substr(i, L));
            if (d < best) {
                best = d;
                if (best == 0) break;
            }
        }
        total += best;
    }
    return total;
}
/**
 * @brief Structure holding a median string and its total distance.
 */
struct MedianItem { string med; int tot; };

/**
 * @brief Find the best 5 median strings with smallest total distances.
 * @param seqs Vector of DNA sequences.
 * @param L Motif length.
 * @return Vector of the best five MedianItem results.
 */

vector<MedianItem> bestFiveMedians(const vector<string>& seqs, int L) {
    vector<MedianItem> all;
    all.reserve(1u << (2*L)); // up to 4^L

    for (const auto& lmer : generateLmers(L)) {
        int tot = totalDistanceToSet(lmer, seqs, L);
        all.push_back({lmer, tot});
    }

    sort(all.begin(), all.end(), [](const MedianItem& a, const MedianItem& b){
        if (a.tot != b.tot) return a.tot < b.tot; // smaller total first
        return a.med < b.med;                     // then lexicographic
    });

    if (all.size() > 5) all.resize(5);
    return all;
}
/**
 * @brief Structure holding consensus results.
 */

struct ConsensusOut {
    string consensus;
    int score = 0;
    vector<pair<int,string>> pos_and_str; // (position 0-based, motif string)
};

/**
 * @brief Build consensus motif for a median string from sequences.
 * @param median The median string.
 * @param seqs Vector of sequences.
 * @param L Motif length.
 * @return ConsensusOut object with consensus string, score, and positions.
 */

ConsensusOut consensusForMedian(const string& median, const vector<string>& seqs, int L) {
    ConsensusOut out;
    if (L <= 0) return out;

    vector<int> A(L,0), C(L,0), G(L,0), T(L,0);

    for (const auto& seq : seqs) {
        if ((int)seq.size() < L) { // keep order with a placeholder
            out.pos_and_str.push_back({0, string(L, '-')});
            continue;
        }
        int bestDist = L+1, bestPos = 0;
        for (size_t i = 0; i + L <= seq.size(); ++i) {
            int d = hammingDistance(median, seq.substr(i, L));
            if (d < bestDist) { bestDist = d; bestPos = (int)i; if (d==0) break; }
        }
        string win = seq.substr(bestPos, L);
        out.pos_and_str.push_back({bestPos, win});
        for (int k = 0; k < L; ++k) {
            switch (win[k]) {
                case 'A': ++A[k]; break;
                case 'C': ++C[k]; break;
                case 'G': ++G[k]; break;
                case 'T': ++T[k]; break;
            }
        }
    }

    out.consensus.assign(L, 'A');
    out.score = 0;
    for (int k = 0; k < L; ++k) {
        int mx = max(max(A[k],C[k]), max(G[k],T[k]));
        out.score += mx;
        if (mx == A[k]) out.consensus[k] = 'A';
        else if (mx == C[k]) out.consensus[k] = 'C';
        else if (mx == G[k]) out.consensus[k] = 'G';
        else out.consensus[k] = 'T';
    }
    return out;
}

/**
 * @brief Main function: manages input/output and overall workflow.
 * @param argc Number of command-line arguments.
 * @param argv Argument values: [program, input_file, L, output_file]
 * @return 0 if successful, non-zero for errors.
 */

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <input_fasta> <L> <output_txt>\n";
        cerr << "Example: " << argv[0] << " HMP-part.fa 6 results.txt\n";
        return 1;
    }
    string input_file = argv[1];
    int L = stoi(argv[2]);
    string output_path = argv[3];
    if (L <= 0) {
        cerr << "ERROR: L must be a positive integer.\n";
        return 1;
    }

    vector<string> sequences = readFasta(input_file);
    if (sequences.empty()) {
        cerr << "ERROR: no sequences read from " << input_file
             << ". Check working directory and file path.\n";
        return 1;
    }

    auto top5 = bestFiveMedians(sequences, L);

    ofstream out(output_path);
    if (!out) {
        cerr << "ERROR: cannot open output file: " << output_path << "\n";
        return 1;
    }

    for (const auto& item : top5) {
        auto cons = consensusForMedian(item.med, sequences, L);
        out << "median string: " << item.med << " (tot_dist = " << item.tot << ")\n";
        out << "motif consensus string: " << cons.consensus
            << " (consensus_score = " << cons.score << ")\n";
        out << "motif positions/strings S=(s1..st):  \n";
        for (size_t i = 0; i < cons.pos_and_str.size(); ++i) {
            if (i == 0)
                out << "    ";  // keep spacing before first item only
            else
                out << ", ";     // only comma and space between items

            out << cons.pos_and_str[i].first
                << "(" << cons.pos_and_str[i].second << ")";
        }
        out << "\n\n";
    }
    out.close();
    return 0;
}
