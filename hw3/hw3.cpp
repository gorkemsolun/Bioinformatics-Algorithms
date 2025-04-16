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

#define INT_MIN_Redefined INT_MIN + 1000000000 // Redefine a minimum integer value for initialization as we may use scores < 0.

// Affine-gap global alignment using three matrices: V, F, and E
// https://docs.google.com/presentation/d/15urWu7TXmdQ5ocRs-k2O5zGPXaHu-oehhpGND3S_SbE/edit?slide=id.p38#slide=id.p38
// V: match/mismatch matrix (merged version of V and G matrices in the slides)
// F: matrix for vertical gaps
// E: matrix for horizontal gaps
void affine_alignment(const string& string1, const string& string2,
    int matchScore, int mismatchScore,
    int gapOpeningScore, int gapExtensionScore,
    int* alignmentScore = nullptr, string* alignmentString1 = nullptr, string* alignmentString2 = nullptr) {

    vector<vector<int>> V(string1.size() + 1, vector<int>(string2.size() + 1, INT_MIN_Redefined));  // V: match/mismatch (merged version of V and G)
    vector<vector<int>> F(string1.size() + 1, vector<int>(string2.size() + 1, INT_MIN_Redefined));  // F: vertical gap (gap in string2)
    vector<vector<int>> E(string1.size() + 1, vector<int>(string2.size() + 1, INT_MIN_Redefined));  // E: horizontal gap (gap in string1)

    // 0 indicates coming from V, 1 from F, 2 from E.
    vector<vector<int>> traceV(string1.size() + 1, vector<int>(string2.size() + 1, -1));
    // 0 indicates coming from V (vertical gap starting), 1 from F (gap extension).
    vector<vector<int>> traceF(string1.size() + 1, vector<int>(string2.size() + 1, -1));
    // 0 indicates coming from V (horizontal gap starting), 1 from E (gap extension).
    vector<vector<int>> traceE(string1.size() + 1, vector<int>(string2.size() + 1, -1));


    V[0][0] = 0;
    F[0][0] = E[0][0] = INT_MIN_Redefined;
    for (size_t i = 1; i <= string1.size(); ++i) {
        V[i][0] = INT_MIN_Redefined;
        F[i][0] = gapOpeningScore + gapExtensionScore * (i - 1);
        traceF[i][0] = (i == 1 ? 0 : 1);
        E[i][0] = INT_MIN_Redefined;
    }
    for (size_t j = 1; j <= string2.size(); ++j) {
        V[0][j] = INT_MIN_Redefined;
        E[0][j] = gapOpeningScore + gapExtensionScore * (j - 1);
        traceE[0][j] = (j == 1 ? 0 : 1);
        F[0][j] = INT_MIN_Redefined;
    }

    for (size_t i = 1; i <= string1.size(); ++i) {
        for (size_t j = 1; j <= string2.size(); ++j) {
            int score_sub = (string1[i - 1] == string2[j - 1]) ? matchScore : mismatchScore;

            V[i][j] = V[i - 1][j - 1] + score_sub;
            traceV[i][j] = 0;
            if (F[i - 1][j - 1] + score_sub > V[i][j]) {
                V[i][j] = F[i - 1][j - 1] + score_sub;
                traceV[i][j] = 1;
            }
            if (E[i - 1][j - 1] + score_sub > V[i][j]) {
                V[i][j] = E[i - 1][j - 1] + score_sub;
                traceV[i][j] = 2;
            }

            F[i][j] = V[i - 1][j] + gapOpeningScore + gapExtensionScore;
            traceF[i][j] = 0;
            if (F[i - 1][j] + gapExtensionScore > F[i][j]) {
                F[i][j] = F[i - 1][j] + gapExtensionScore;
                traceF[i][j] = 1;
            }

            E[i][j] = V[i][j - 1] + gapOpeningScore + gapExtensionScore;
            traceE[i][j] = 0;
            if (E[i][j - 1] + gapExtensionScore > E[i][j]) {
                E[i][j] = E[i][j - 1] + gapExtensionScore;
                traceE[i][j] = 1;
            }
        }
    }

    int finalState = 0; // 0: V, 1: F, 2: E
    int bestScore = V[string1.size()][string2.size()];
    if (F[string1.size()][string2.size()] > bestScore) {
        bestScore = F[string1.size()][string2.size()];
        finalState = 1;
    }
    if (E[string1.size()][string2.size()] > bestScore) {
        bestScore = E[string1.size()][string2.size()];
        finalState = 2;
    }
    if (alignmentScore) {
        *alignmentScore = bestScore;
    }

    if (alignmentString1 == nullptr || alignmentString2 == nullptr) {
        return;
    }
    int i = string1.size(), j = string2.size();
    int state = finalState;
    while (i > 0 || j > 0) {
        if (state == 0) { // V
            int prev = traceV[i][j];
            alignmentString1->push_back(string1[i - 1]);
            alignmentString2->push_back(string2[j - 1]);
            --i; --j;
            state = prev;
        } else if (state == 1) { // F
            if (traceF[i][j] == 0) {
                state = 0;
            } else {
                state = 1;
            }
            alignmentString1->push_back(string1[i - 1]);
            alignmentString2->push_back('-');
            --i;
        } else { // E 
            if (traceE[i][j] == 0) {
                state = 0;
            } else {
                state = 2;
            }
            alignmentString1->push_back('-');
            alignmentString2->push_back(string2[j - 1]);
            --j;
        }
    }

    reverse(alignmentString1->begin(), alignmentString1->end());
    reverse(alignmentString2->begin(), alignmentString2->end());
}

vector<pair<string, string>> readFASTA(const string& filename) {
    ifstream in(filename);
    if (!in) {
        cout << "Error: Could not open file " << filename << endl;
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

int main(int argc, char* argv[]) {
    if (argc < 7) {
        cout << "Usage: " << argv[0] << " -i input.fasta -o output.phy -s matchScore:mismatchScore:gapOpeningScore:gapExtensionScore" << endl;
        return 0;
    }

    string inputFile, outputFile, scoresStr;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            inputFile = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            outputFile = argv[++i];
        } else if (arg == "-s" && i + 1 < argc) {
            scoresStr = argv[++i];
        } else {
            cout << "Unknown argument: " << arg << endl;
            return 0;
        }
    }

    int matchScore, mismatchScore, gapOpeningScore, gapExtensionScore;
    {
        vector<int> values;
        stringstream ss(scoresStr);
        string token;
        while (getline(ss, token, ':')) {
            values.push_back(stoi(token));
        }
        if (values.size() != 4) {
            cout << "Error: Score must have four values separated by ':'" << endl;
            return 0;
        }
        matchScore = values[0];
        mismatchScore = values[1];
        gapOpeningScore = values[2];
        gapExtensionScore = values[3];
    }

    /* string inputFile = "input.fasta";
    string outputFile = "output.phy";
    int matchScore = 5;
    int mismatchScore = -4;
    int gapOpeningScore = -16;
    int gapExtensionScore = -4; */

    vector<pair<string, string>> fasta = readFASTA(inputFile);
    if (fasta.size() == 0) {
        cout << "No sequences found in " << inputFile << endl;
        return 0;
    }

    // if there is only one sequence
    if (fasta.size() == 1) {
        ofstream out(outputFile);
        string sequence = fasta[0].second;
        out << "1 " << sequence.size() << "\n";
        out << fasta[0].first << " " << sequence << "\n";
        out.close();
        return 0;
    }

    // Compute sim(Si, Sj) for every pair (i,j)
    vector<int> sumScores(fasta.size(), 0);
    for (size_t i = 0; i < fasta.size(); ++i) {
        for (size_t j = i + 1; j < fasta.size(); ++j) {
            int alignmentScore = 0;
            affine_alignment(fasta[i].second, fasta[j].second, matchScore, mismatchScore, gapOpeningScore, gapExtensionScore, &alignmentScore);
            // Compute star_score(i) for every i
            sumScores[i] += alignmentScore;
            sumScores[j] += alignmentScore;
        }
    }

    // Choose the index c that maximizes star_score(c)  and make it the center of the star
    size_t centerSequenceIndex = 0;
    int bestSum = sumScores[0];
    for (size_t i = 1; i < fasta.size(); ++i) {
        if (sumScores[i] > bestSum) {
            bestSum = sumScores[i];
            centerSequenceIndex = i;
        }
    }

    // Produce a multiple alignment such that, for every i, the induced pairwise alignment of Sc and Si is the same as the optimum alignment of Sc and Si.
    vector<vector<int>> gapPatterns(fasta.size(), vector<int>(fasta[centerSequenceIndex].second.size() + 1, 0));
    vector<string*> alignedOthers(fasta.size(), nullptr);

    // Others
    for (size_t i = 0; i < fasta.size(); ++i) {
        if (i == centerSequenceIndex) {
            continue;
        }
        string* alignedCenter = new string;
        string* alignedOther = new string;
        affine_alignment(fasta[centerSequenceIndex].second, fasta[i].second,
            matchScore, mismatchScore, gapOpeningScore, gapExtensionScore, nullptr, alignedCenter, alignedOther);
        cout << "Alignment of " << fasta[centerSequenceIndex].first << " and " << fasta[i].first << ":" << endl;
        cout << *alignedCenter << "\n" << *alignedOther << endl;

        int position = 0;
        for (size_t j = 0; j < alignedCenter->size(); ++j) {
            if ((*alignedCenter)[j] == '-') {
                ++gapPatterns[i][position];
            } else {
                ++position;
            }
        }

        delete alignedCenter;
        alignedOthers[i] = alignedOther;
    }

    // Merge gap patterns from all pairwise alignments
    vector<int> mergedGap(fasta[centerSequenceIndex].second.size() + 1, 0);
    for (size_t i = 0; i < fasta.size(); ++i) {
        for (size_t k = 0; k <= fasta[centerSequenceIndex].second.size(); ++k) {
            mergedGap[k] = max(mergedGap[k], gapPatterns[i][k]);
        }
    }

    // Final alignment for each sequence
    // Final alignment for the center sequence
    vector<string> finalAlignment(fasta.size(), "");
    for (size_t k = 0; k <= fasta[centerSequenceIndex].second.size(); ++k) {
        finalAlignment[centerSequenceIndex].append(mergedGap[k] - gapPatterns[centerSequenceIndex][k], '-');
        if (k < fasta[centerSequenceIndex].second.size()) {
            finalAlignment[centerSequenceIndex].push_back(fasta[centerSequenceIndex].second[k]);
        }
    }
    // Combining the aligned sequences with the merged gap pattern
    for (size_t i = 0; i < fasta.size(); ++i) {
        if (i == centerSequenceIndex) {
            continue;
        }
        for (size_t k = 0; k <= alignedOthers[i]->size(); ++k) {
            if (k < alignedOthers[i]->size()) {
                finalAlignment[i].push_back((*alignedOthers[i])[k]);
            }
            finalAlignment[i].append(mergedGap[k] - gapPatterns[i][k], '-');
        }
    }

    // Swap the center sequence with the first sequence in the output
    swap(finalAlignment[0], finalAlignment[centerSequenceIndex]);
    swap(fasta[0], fasta[centerSequenceIndex]);

    ofstream out(outputFile);
    if (!out) {
        cout << "Error: Could not open output file " << outputFile << endl;
        return 0;
    }

    out << fasta.size() << " " << finalAlignment[centerSequenceIndex].size() << "\n";
    for (size_t i = 0; i < fasta.size(); ++i) {
        string sequenceId = fasta[i].first;
        if (sequenceId.size() > 10) {
            sequenceId = sequenceId.substr(0, 10);
        } else if (sequenceId.size() < 10) {
            sequenceId.append(10 - sequenceId.size(), ' ');
        }
        out << sequenceId;

        // Add a space every 10 characters
        out << finalAlignment[i][0];
        for (size_t j = 1; j < finalAlignment[i].size(); ++j) {
            if (j % 10 == 0) {
                out << " ";
            }
            out << finalAlignment[i][j];
        }
        out << "\n";
    }
    out.close();

    return 0;
}
