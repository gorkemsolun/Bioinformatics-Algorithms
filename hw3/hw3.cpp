// Center star alinment. Output in PHYLIP format.
// CS481 Görkem Kadir Solun 22003214

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <limits>

using namespace std;

#define INT_MIN -2000000000 // Redefine a minimum integer value for initialization as we may use scores < 0.

// Affine-gap global alignment using three matrices: M, F, and E
// https://docs.google.com/presentation/d/15urWu7TXmdQ5ocRs-k2O5zGPXaHu-oehhpGND3S_SbE/edit?slide=id.p38#slide=id.p38
// M: match/mismatch matrix (merged version of V and G matrices in the slides)
// F: matrix for vertical gaps
// E: matrix for horizontal gaps
pair<string, string> affine_alignment(const string& string1, const string& string2,
    int match, int mismatch,
    int gap_open, int gap_extend,
    int* alignment_score_out = nullptr) {

    vector<vector<int>> M(string1.size() + 1, vector<int>(string2.size() + 1, INT_MIN));  // M: match/mismatch (merged version of V and G)
    vector<vector<int>> F(string1.size() + 1, vector<int>(string2.size() + 1, INT_MIN));  // F: vertical gap (gap in string2)
    vector<vector<int>> E(string1.size() + 1, vector<int>(string2.size() + 1, INT_MIN));  // E: horizontal gap (gap in string1)

    // For M: 0 indicates coming from M, 1 from F, 2 from E.
    vector<vector<int>> traceM(string1.size() + 1, vector<int>(string2.size() + 1, -1));
    // For F: 0 indicates coming from M (vertical gap starting), 1 from F (gap extension).
    vector<vector<int>> traceF(string1.size() + 1, vector<int>(string2.size() + 1, -1));
    // For E: 0 indicates coming from M (horizontal gap starting), 1 from E (gap extension).
    vector<vector<int>> traceE(string1.size() + 1, vector<int>(string2.size() + 1, -1));


    M[0][0] = 0;
    F[0][0] = E[0][0] = INT_MIN;
    for (int i = 1; i <= string1.size(); ++i) {
        M[i][0] = INT_MIN;
        F[i][0] = gap_open + gap_extend * (i - 1);
        traceF[i][0] = (i == 1 ? 0 : 1);
        E[i][0] = INT_MIN;
    }
    for (int j = 1; j <= string2.size(); ++j) {
        M[0][j] = INT_MIN;
        E[0][j] = gap_open + gap_extend * (j - 1);
        traceE[0][j] = (j == 1 ? 0 : 1);
        F[0][j] = INT_MIN;
    }

    for (int i = 1; i <= string1.size(); ++i) {
        for (int j = 1; j <= string2.size(); ++j) {
            int score_sub = (string1[i - 1] == string2[j - 1]) ? match : mismatch;

            M[i][j] = M[i - 1][j - 1] + score_sub;
            traceM[i][j] = 0;
            if (F[i - 1][j - 1] + score_sub > M[i][j]) {
                M[i][j] = F[i - 1][j - 1] + score_sub;
                traceM[i][j] = 1;
            }
            if (E[i - 1][j - 1] + score_sub > M[i][j]) {
                M[i][j] = E[i - 1][j - 1] + score_sub;
                traceM[i][j] = 2;
            }

            F[i][j] = M[i - 1][j] + gap_open + gap_extend;
            traceF[i][j] = 0;
            if (F[i - 1][j] + gap_extend > F[i][j]) {
                F[i][j] = F[i - 1][j] + gap_extend;
                traceF[i][j] = 1;
            }

            E[i][j] = M[i][j - 1] + gap_open + gap_extend;
            traceE[i][j] = 0;
            if (E[i][j - 1] + gap_extend > E[i][j]) {
                E[i][j] = E[i][j - 1] + gap_extend;
                traceE[i][j] = 1;
            }
        }
    }

    int finalState = 0; // 0: M, 1: F, 2: E
    int bestScore = M[string1.size()][string2.size()];
    if (F[string1.size()][string2.size()] > bestScore) {
        bestScore = F[string1.size()][string2.size()];
        finalState = 1;
    }
    if (E[string1.size()][string2.size()] > bestScore) {
        bestScore = E[string1.size()][string2.size()];
        finalState = 2;
    }
    if (alignment_score_out) {
        *alignment_score_out = bestScore;
    }

    int i = string1.size(), j = string2.size();
    string alignment1 = "", alignment2 = "";
    int state = finalState;
    while (i > 0 || j > 0) {
        if (state == 0) { // M
            int prev = traceM[i][j];
            alignment1.push_back(string1[i - 1]);
            alignment2.push_back(string2[j - 1]);
            --i; --j;
            state = prev;
        } else if (state == 1) { // F
            if (traceF[i][j] == 0) {
                state = 0;
            } else {
                state = 1;
            }
            alignment1.push_back(string1[i - 1]);
            alignment2.push_back('-');
            --i;
        } else { // E 
            if (traceE[i][j] == 0) {
                state = 0;
            } else {
                state = 2;
            }
            alignment1.push_back('-');
            alignment2.push_back(string2[j - 1]);
            --j;
        }
    }

    reverse(alignment1.begin(), alignment1.end());
    reverse(alignment2.begin(), alignment2.end());
    return { alignment1, alignment2 };
}

// Read a FASTA file and return a vector of (name, sequence) pairs.
vector<pair<string, string>> readFASTA(const string& filename) {
    ifstream in(filename);
    if (!in) {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }

    vector<pair<string, string>> sequences;
    string line, header, sequence;
    while (getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            if (!header.empty()) {
                sequences.push_back({ header, sequence });
                sequence = "";
            }
            header = line.substr(1); // Remove '>' character
        } else {
            for (char c : line)
                if (!isspace(c)) {
                    sequence.push_back(c);
                }
        }
    }
    if (!header.empty()) {
        sequences.push_back({ header, sequence });
    }
    return sequences;
}

// Given an aligned center string and the original center,
// compute the gap pattern as a vector<int> for each gap slot.
vector<int> compute_gap_pattern(const string& alignedCenter, const string& center) {
    vector<int> pattern(center.size() + 1, 0);
    int position = 0;
    for (size_t i = 0; i < alignedCenter.size(); ++i) {
        if (alignedCenter[i] == '-') {
            ++pattern[position];
        } else {
            ++position;
        }
    }
    return pattern;
}

// From a pairwise alignment, extract the mapping for the other sequence.
vector<char> compute_other_mapping(const string& alignedCenter, const string& alignedOther, const string& center) {
    vector<char> mapping;
    int position = 0;
    for (size_t i = 0; i < alignedCenter.size(); ++i) {
        if (alignedCenter[i] != '-') {
            mapping.push_back(alignedOther[i]);
            ++position;
        }
    }
    return mapping;
}

// Merge a sequence with its gap pattern to produce the final aligned string.
string merge_with_master(const string& sequence, const vector<int>& pattern, const vector<int>& masterGap) {
    string result;
    for (int k = 0; k <= sequence.size(); ++k) {
        int extra = masterGap[k] - pattern[k];
        result.append(extra, '-');
        if (k < sequence.size()) {
            result.push_back(sequence[k]);
        }
    }
    return result;
}

// Merge the mapping for sequences aligned to the center with the master gap pattern.
string merge_with_master_mapping(const vector<char>& mapping, const vector<int>& pattern, const vector<int>& masterGap) {
    string result;
    for (int k = 0; k <= mapping.size(); ++k) {
        int extra = masterGap[k] - pattern[k];
        result.append(extra, '-');
        if (k < mapping.size()) {
            result.push_back(mapping[k]);
        }
    }
    return result;
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        cerr << "Usage: " << argv[0] << " -i input.fasta -o output.phy -s match:mismatch:gap_open:gap_extend" << endl;
        return 1;
    }

    string inputFile, outputFile, scoresStr;
    // Command-line argument parsing.
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            inputFile = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            outputFile = argv[++i];
        } else if (arg == "-s" && i + 1 < argc) {
            scoresStr = argv[++i];
        } else {
            cerr << "Unknown argument: " << arg << endl;
            return 1;
        }
    }

    // Parse the score string, expected in the format "match:mismatch:gap_open:gap_extend"
    int match, mismatch, gap_open, gap_extend;
    {
        vector<int> vals;
        stringstream ss(scoresStr);
        string token;
        while (getline(ss, token, ':')) {
            vals.push_back(stoi(token));
        }
        if (vals.size() != 4) {
            cerr << "Error: Score must have four values separated by ':'" << endl;
            return 1;
        }
        match = vals[0];
        mismatch = vals[1];
        gap_open = vals[2];
        gap_extend = vals[3];
    }

    // Read the FASTA sequences.
    vector<pair<string, string>> fasta = readFASTA(inputFile);
    if (fasta.size() == 0) {
        cerr << "No sequences found in " << inputFile << endl;
        return 1;
    }

    // If there is only one sequence, output it in PHYLIP format.
    if (fasta.size() == 1) {
        ofstream out(outputFile);
        string sequence = fasta[0].second;
        out << "1 " << sequence.size() << "\n";
        out << fasta[0].first << " " << sequence << "\n";
        out.close();
        return 0;
    }

    // Step 1: Select the center sequence using pairwise alignment scores.
    vector<int> sumScores(fasta.size(), 0);
    for (int i = 0; i < fasta.size(); ++i) {
        for (int j = i + 1; j < fasta.size(); ++j) {
            int score = 0;
            affine_alignment(fasta[i].second, fasta[j].second, match, mismatch, gap_open, gap_extend, &score);
            sumScores[i] += score;
            sumScores[j] += score;
        }
    }

    int centerIndex = 0;
    int bestSum = sumScores[0];
    for (int i = 1; i < fasta.size(); ++i) {
        if (sumScores[i] > bestSum) {
            bestSum = sumScores[i];
            centerIndex = i;
        }
    }

    string centerSeq = fasta[centerIndex].second;
    int centerLength = centerSeq.size();

    // Step 2: Align all sequences to the center.
    vector< vector<int> > gapPatterns(fasta.size(), vector<int>(centerLength + 1, 0));
    vector< vector<char> > mappings(fasta.size(), vector<char>());

    // For the center itself.
    mappings[centerIndex].resize(centerLength);
    for (int k = 0; k < centerLength; ++k) {
        mappings[centerIndex][k] = centerSeq[k];
    }

    // For every other sequence, perform pairwise alignment with the center.
    for (int i = 0; i < fasta.size(); ++i) {
        if (i == centerIndex) {
            continue;
        }
        auto alignmentPair = affine_alignment(centerSeq, fasta[i].second,
            match, mismatch, gap_open, gap_extend);
        string alignedCenter = alignmentPair.first;
        string alignedOther = alignmentPair.second;
        gapPatterns[i] = compute_gap_pattern(alignedCenter, centerSeq);
        mappings[i] = compute_other_mapping(alignedCenter, alignedOther, centerSeq);
    }

    // Step 3: Merge gap patterns from all pairwise alignments to create the master gap pattern.
    vector<int> masterGap(centerLength + 1, 0);
    for (int i = 0; i < fasta.size(); ++i) {
        for (int k = 0; k <= centerLength; ++k) {
            masterGap[k] = max(masterGap[k], gapPatterns[i][k]);
        }
    }

    // Step 4: Build the final aligned strings for all sequences using the master gap pattern.
    vector<string> finalAlign(fasta.size(), "");
    finalAlign[centerIndex] = merge_with_master(centerSeq, gapPatterns[centerIndex], masterGap);
    for (int i = 0; i < fasta.size(); ++i) {
        if (i == centerIndex) { continue; }
        finalAlign[i] = merge_with_master_mapping(mappings[i], gapPatterns[i], masterGap);
    }

    int finalLength = finalAlign[centerIndex].size();

    // Step 5: Write the multiple sequence alignment in PHYLIP format.
    ofstream out(outputFile);
    if (!out) {
        cerr << "Error: Could not open output file " << outputFile << endl;
        return 1;
    }
    out << fasta.size() << " " << finalLength << "\n";
    for (int i = 0; i < fasta.size(); ++i) {
        out << fasta[i].first << " " << finalAlign[i] << "\n";
    }
    out.close();

    return 0;
}
