// Center star alinment
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

// Affine-gap global alignment using three matrices: M, F, and E
// M: match/mismatch matrix
// F: matrix for vertical gaps (i.e. gaps in the second sequence)
// E: matrix for horizontal gaps (i.e. gaps in the first sequence)
// The function returns a pair of aligned strings (first: alignment for s1, second: alignment for s2).
// Optionally, the final alignment score is stored in the provided pointer.
pair<string, string> affine_alignment(const string& s1, const string& s2,
    int match, int mismatch,
    int gap_open, int gap_extend,
    int* alignment_score_out = nullptr) {
    int m = s1.size(), n = s2.size();
    // Initialize the three matrices with dimensions (m+1) x (n+1)
    vector<vector<int>> M(m + 1, vector<int>(n + 1, INT_MIN));
    vector<vector<int>> F(m + 1, vector<int>(n + 1, INT_MIN));  // F: vertical gap (gap in s2)
    vector<vector<int>> E(m + 1, vector<int>(n + 1, INT_MIN));  // E: horizontal gap (gap in s1)

    // Traceback matrices:
    // For M: 0 indicates coming from M, 1 from F, 2 from E.
    vector<vector<int>> traceM(m + 1, vector<int>(n + 1, -1));
    // For F: 0 indicates coming from M (vertical gap starting), 1 from F (gap extension).
    vector<vector<int>> traceF(m + 1, vector<int>(n + 1, -1));
    // For E: 0 indicates coming from M (horizontal gap starting), 1 from E (gap extension).
    vector<vector<int>> traceE(m + 1, vector<int>(n + 1, -1));

    // Initialization
    M[0][0] = 0;
    F[0][0] = E[0][0] = INT_MIN;

    // Initialize first column for F (vertical gaps) and M for i > 0.
    for (int i = 1; i <= m; i++) {
        M[i][0] = INT_MIN;
        F[i][0] = gap_open + gap_extend * (i - 1);
        traceF[i][0] = (i == 1 ? 0 : 1); // For i==1, coming from M[0][0]; afterwards, extension.
        E[i][0] = INT_MIN;
    }
    // Initialize first row for E (horizontal gaps) and M for j > 0.
    for (int j = 1; j <= n; j++) {
        M[0][j] = INT_MIN;
        E[0][j] = gap_open + gap_extend * (j - 1);
        traceE[0][j] = (j == 1 ? 0 : 1);
        F[0][j] = INT_MIN;
    }

    // Fill DP matrices
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int score_sub = (s1[i - 1] == s2[j - 1]) ? match : mismatch;

            // Compute M[i][j]: from diagonal (match/mismatch) move from M, F, or E.
            int score_from_M = M[i - 1][j - 1] + score_sub;
            int score_from_F = F[i - 1][j - 1] + score_sub;
            int score_from_E = E[i - 1][j - 1] + score_sub;
            M[i][j] = score_from_M;
            traceM[i][j] = 0;
            if (score_from_F > M[i][j]) {
                M[i][j] = score_from_F;
                traceM[i][j] = 1;
            }
            if (score_from_E > M[i][j]) {
                M[i][j] = score_from_E;
                traceM[i][j] = 2;
            }

            // Compute F[i][j]: vertical gap (gap in s2)
            int scoreF_from_M = M[i - 1][j] + gap_open + gap_extend;
            int scoreF_from_F = F[i - 1][j] + gap_extend;
            F[i][j] = scoreF_from_M;
            traceF[i][j] = 0;
            if (scoreF_from_F > F[i][j]) {
                F[i][j] = scoreF_from_F;
                traceF[i][j] = 1;
            }

            // Compute E[i][j]: horizontal gap (gap in s1)
            int scoreE_from_M = M[i][j - 1] + gap_open + gap_extend;
            int scoreE_from_E = E[i][j - 1] + gap_extend;
            E[i][j] = scoreE_from_M;
            traceE[i][j] = 0;
            if (scoreE_from_E > E[i][j]) {
                E[i][j] = scoreE_from_E;
                traceE[i][j] = 1;
            }
        }
    }

    // Choose the best final score among M, F, and E at (m, n)
    int finalState = 0; // 0: M, 1: F, 2: E
    int bestScore = M[m][n];
    if (F[m][n] > bestScore) { bestScore = F[m][n]; finalState = 1; }
    if (E[m][n] > bestScore) { bestScore = E[m][n]; finalState = 2; }
    if (alignment_score_out)
        *alignment_score_out = bestScore;

    // Traceback to recover the alignment
    int i = m, j = n;
    string align1 = "", align2 = "";
    int state = finalState;
    while (i > 0 || j > 0) {
        if (state == 0) { // M matrix: diagonal move
            int prev = traceM[i][j];
            align1.push_back(s1[i - 1]);
            align2.push_back(s2[j - 1]);
            i--; j--;
            state = prev; // previous state (could be 0, 1, or 2)
        } else if (state == 1) { // F matrix: vertical gap (gap in s2)
            if (traceF[i][j] == 0)
                state = 0;
            else
                state = 1;
            align1.push_back(s1[i - 1]);
            align2.push_back('-');
            i--;
        } else { // state == 2, E matrix: horizontal gap (gap in s1)
            if (traceE[i][j] == 0)
                state = 0;
            else
                state = 2;
            align1.push_back('-');
            align2.push_back(s2[j - 1]);
            j--;
        }
    }

    reverse(align1.begin(), align1.end());
    reverse(align2.begin(), align2.end());
    return { align1, align2 };
}

// Read a FASTA file and return a vector of (name, sequence) pairs.
vector<pair<string, string>> readFASTA(const string& filename) {
    ifstream in(filename);
    if (!in) {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }

    vector<pair<string, string>> seqs;
    string line, header, seq;
    while (getline(in, line)) {
        if (line.empty())
            continue;
        if (line[0] == '>') {
            if (!header.empty()) {
                seqs.push_back({ header, seq });
                seq = "";
            }
            header = line.substr(1); // Remove '>' character
        } else {
            for (char c : line)
                if (!isspace(c))
                    seq.push_back(c);
        }
    }
    if (!header.empty())
        seqs.push_back({ header, seq });
    return seqs;
}

// Given an aligned center string and the original center,
// compute the gap pattern as a vector<int> for each gap slot.
vector<int> compute_gap_pattern(const string& alignedCenter, const string& center) {
    int n = center.size();
    vector<int> pattern(n + 1, 0);
    int pos = 0;
    for (size_t i = 0; i < alignedCenter.size(); i++) {
        if (alignedCenter[i] == '-') {
            pattern[pos]++;
        } else {
            pos++;
        }
    }
    return pattern;
}

// From a pairwise alignment, extract the mapping for the other sequence.
vector<char> compute_other_mapping(const string& alignedCenter, const string& alignedOther, const string& center) {
    vector<char> mapping;
    int pos = 0;
    for (size_t i = 0; i < alignedCenter.size(); i++) {
        if (alignedCenter[i] != '-') {
            mapping.push_back(alignedOther[i]);
            pos++;
        }
    }
    return mapping;
}

// Merge a sequence with its gap pattern to produce the final aligned string.
string merge_with_master(const string& seq, const vector<int>& pattern, const vector<int>& masterGap) {
    string result;
    int n = seq.size();
    for (int k = 0; k <= n; k++) {
        int extra = masterGap[k] - pattern[k];
        result.append(extra, '-');
        if (k < n)
            result.push_back(seq[k]);
    }
    return result;
}

// Merge the mapping for sequences aligned to the center with the master gap pattern.
string merge_with_master_mapping(const vector<char>& mapping, const vector<int>& pattern, const vector<int>& masterGap) {
    string result;
    int n = mapping.size();
    for (int k = 0; k <= n; k++) {
        int extra = masterGap[k] - pattern[k];
        result.append(extra, '-');
        if (k < n)
            result.push_back(mapping[k]);
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
    for (int i = 1; i < argc; i++) {
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
    int numSeq = fasta.size();
    if (numSeq == 0) {
        cerr << "No sequences found in " << inputFile << endl;
        return 1;
    }

    // If there is only one sequence, output it in PHYLIP format.
    if (numSeq == 1) {
        ofstream out(outputFile);
        string seq = fasta[0].second;
        out << "1 " << seq.size() << "\n";
        out << fasta[0].first << " " << seq << "\n";
        out.close();
        return 0;
    }

    // Step 1: Select the center sequence using pairwise alignment scores.
    vector<int> sumScores(numSeq, 0);
    for (int i = 0; i < numSeq; i++) {
        for (int j = i + 1; j < numSeq; j++) {
            int score = 0;
            affine_alignment(fasta[i].second, fasta[j].second, match, mismatch, gap_open, gap_extend, &score);
            sumScores[i] += score;
            sumScores[j] += score;
        }
    }

    int centerIndex = 0;
    int bestSum = sumScores[0];
    for (int i = 1; i < numSeq; i++) {
        if (sumScores[i] > bestSum) {
            bestSum = sumScores[i];
            centerIndex = i;
        }
    }

    string centerSeq = fasta[centerIndex].second;
    int centerLength = centerSeq.size();

    // Step 2: Align all sequences to the center.
    vector< vector<int> > gapPatterns(numSeq, vector<int>(centerLength + 1, 0));
    vector< vector<char> > mappings(numSeq, vector<char>());

    // For the center itself.
    mappings[centerIndex].resize(centerLength);
    for (int k = 0; k < centerLength; k++) {
        mappings[centerIndex][k] = centerSeq[k];
    }

    // For every other sequence, perform pairwise alignment with the center.
    for (int i = 0; i < numSeq; i++) {
        if (i == centerIndex) continue;
        auto alignmentPair = affine_alignment(centerSeq, fasta[i].second,
            match, mismatch, gap_open, gap_extend);
        string alignedCenter = alignmentPair.first;
        string alignedOther = alignmentPair.second;
        gapPatterns[i] = compute_gap_pattern(alignedCenter, centerSeq);
        mappings[i] = compute_other_mapping(alignedCenter, alignedOther, centerSeq);
    }

    // Step 3: Merge gap patterns from all pairwise alignments to create the master gap pattern.
    vector<int> masterGap(centerLength + 1, 0);
    for (int i = 0; i < numSeq; i++) {
        for (int k = 0; k <= centerLength; k++) {
            masterGap[k] = max(masterGap[k], gapPatterns[i][k]);
        }
    }

    // Step 4: Build the final aligned strings for all sequences using the master gap pattern.
    vector<string> finalAlign(numSeq, "");
    finalAlign[centerIndex] = merge_with_master(centerSeq, gapPatterns[centerIndex], masterGap);
    for (int i = 0; i < numSeq; i++) {
        if (i == centerIndex) continue;
        finalAlign[i] = merge_with_master_mapping(mappings[i], gapPatterns[i], masterGap);
    }

    int finalLength = finalAlign[centerIndex].size();

    // Step 5: Write the multiple sequence alignment in PHYLIP format.
    ofstream out(outputFile);
    if (!out) {
        cerr << "Error: Could not open output file " << outputFile << endl;
        return 1;
    }
    out << numSeq << " " << finalLength << "\n";
    for (int i = 0; i < numSeq; i++) {
        out << fasta[i].first << " " << finalAlign[i] << "\n";
    }
    out.close();

    return 0;
}
