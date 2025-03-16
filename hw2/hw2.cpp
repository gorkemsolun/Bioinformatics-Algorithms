// Görkem Kadir Solun 22003214 CS481 HW2 
// This program reads two FASTA files (patterns and texts) and performs global or local alignment between each pair.
// The alignment is done using the Needleman-Wunsch (global) or Smith-Waterman (local) algorithms.
// The program outputs the alignment with the highest score (local) or longest overlap (global) to a file.
// The alignment results include the aligned sequences, alignment score, CIGAR string, and MD:Z string.
// The MD:Z string is a simple representation of mismatches and deletions in the alignment.

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
using namespace std;

struct AlignmentResult {
    int score;
    string alignedPattern;
    string alignedReference;
    string cigar;
    string mdz;
};

vector<string> readFasta(const string& filename) {
    vector<string> sequences;
    ifstream input_file(filename.c_str());
    if (!input_file) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    string line, sequence = "";
    while (getline(input_file, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            if (!sequence.empty()) {
                sequences.push_back(sequence);
                sequence = "";
            }
            continue;
        } else {
            sequence += line;
        }
    }
    if (!sequence.empty()) {
        sequences.push_back(sequence);
    }
    input_file.close();
    return sequences;
}

string prepareCigarString(const vector<char>& alignmentResult) {
    if (alignmentResult.empty()) {
        return "";
    }
    string cigar = "";
    int count = 1;
    char current = alignmentResult.back();
    // NOTE: Start from last element -1, as the input is reversed.
    for (int i = alignmentResult.size() - 2; i >= 0; --i) {
        if (alignmentResult[i] == current) {
            ++count;
        } else {
            cigar += to_string(count) + current;
            current = alignmentResult[i];
            count = 1;
        }
    }
    cigar += to_string(count) + current;
    return cigar;
}

// TODO: Fix this function to handle deletions correctly.
string prepareMDZString(const string& alignedPattern, const string& alignedReference, vector<char>& tracebacks) {
    reverse(tracebacks.begin(), tracebacks.end());
    cout << "tracebacks: " << string(tracebacks.begin(), tracebacks.end()) << endl;
    string mdz = "";
    int matchCount = 0;
    size_t i = 0;
    while (i < tracebacks.size()) {
        if (tracebacks[i] == 'M') {
            // Diagonal move, match or mismatch.
            if (alignedPattern[i] == alignedReference[i]) {
                ++matchCount;
            } else {
                // Mismatch
                mdz += to_string(matchCount);
                mdz.push_back(alignedReference[i]);
                matchCount = 0;
            }
            ++i;
        } else if (tracebacks[i] == 'D') {
            // Left move, deletion.
            mdz += to_string(matchCount);
            mdz.push_back('^');
            matchCount = 0;
            string deletedBases = "";
            while (i < tracebacks.size() && tracebacks[i] == 'D') {
                deletedBases.push_back(alignedReference[i]);
                ++i;
            }
            mdz += deletedBases;
        } else {
            // Up move, insertion or unexpected character in CIGAR/traceback.
            ++i;
        }
    }
    mdz += to_string(matchCount);
    return mdz;
}

AlignmentResult* globalAlignmentNeedlemanWunsch(const string& patterns, const string& references, int matchScore, int mismatchScore, int gapPenalty) {
    vector<vector<int>> dp(patterns.size() + 1, vector<int>(references.size() + 1, 0));
    vector<vector<char>> traceback(patterns.size() + 1, vector<char>(references.size() + 1, ' '));

    AlignmentResult* result = new AlignmentResult();

    // For global alignment, we need to initialize the first row and column with gap penalties.
    for (size_t i = 0; i <= patterns.size(); ++i) {
        dp[i][0] = i * gapPenalty;
        if (i > 0) {
            traceback[i][0] = 'u';
        }
    }
    for (size_t j = 0; j <= references.size(); ++j) {
        dp[0][j] = j * gapPenalty;
        if (j > 0) {
            traceback[0][j] = 'l';
        }
    }

    for (size_t i = 1; i <= patterns.size(); ++i) {
        for (size_t j = 1; j <= references.size(); ++j) {
            int up = dp[i - 1][j] + gapPenalty;
            int left = dp[i][j - 1] + gapPenalty;
            dp[i][j] = dp[i - 1][j - 1] + (patterns[i - 1] == references[j - 1] ? matchScore : mismatchScore); // Diagonal move

            // Record
            traceback[i][j] = 'd';
            if (left > dp[i][j]) {
                dp[i][j] = left;
                traceback[i][j] = 'l';
            }
            if (up > dp[i][j]) {
                dp[i][j] = up;
                traceback[i][j] = 'u';
            }

        }
    }

    size_t trace_i = patterns.size(), trace_j = references.size(); // For global alignment, start from the end. (bottom-right)
    result->alignedPattern = "";
    result->alignedReference = "";
    vector<char> tracebacks;

    while (trace_i > 0 || trace_j > 0) {
        if (trace_i > 0 && trace_j > 0 && traceback[trace_i][trace_j] == 'd') {
            result->alignedPattern.push_back(patterns[trace_i - 1]);
            result->alignedReference.push_back(references[trace_j - 1]);
            tracebacks.push_back('M');
            --trace_i;
            --trace_j;
        } else if (trace_i > 0 && traceback[trace_i][trace_j] == 'u') {
            result->alignedPattern.push_back(patterns[trace_i - 1]);
            result->alignedReference.push_back('-');
            tracebacks.push_back('I');
            --trace_i;
        } else if (trace_j > 0 && traceback[trace_i][trace_j] == 'l') {
            result->alignedPattern.push_back('-');
            result->alignedReference.push_back(references[trace_j - 1]);
            tracebacks.push_back('D');
            --trace_j;
        }
    }

    reverse(result->alignedPattern.begin(), result->alignedPattern.end());
    reverse(result->alignedReference.begin(), result->alignedReference.end());

    result->score = dp[patterns.size()][references.size()];
    result->cigar = prepareCigarString(tracebacks); // Prepare CIGAR string from traceback(reverse).
    result->mdz = prepareMDZString(result->alignedPattern, result->alignedReference, tracebacks);
    return result;
}

AlignmentResult* localAlignmentSmithWaterman(const string& patterns, const string& references, int matchScore, int mismatchScore, int gapPenalty) {
    vector<vector<int>> dp(patterns.size() + 1, vector<int>(references.size() + 1, 0));
    vector<vector<char>> traceback(patterns.size() + 1, vector<char>(references.size() + 1, ' '));

    // For local alignment, we need to initialize the first row and column with 0.
    // The rest of the matrix is initialized with -infinity. (used 0 for simplicity)

    AlignmentResult* result = new AlignmentResult();

    // For local
    result->score = 0; // Keep track of the maximum score in local alignment as we need to find the best alignment.
    size_t trace_i = 0, trace_j = 0; // Keep track of the position of the maximum score to start the traceback.

    for (size_t i = 1; i <= patterns.size(); ++i) {
        for (size_t j = 1; j <= references.size(); ++j) {
            // Compute scores for possible moves. 0, diagonal, up, left respectively.
            int diagonal = dp[i - 1][j - 1] + (patterns[i - 1] == references[j - 1] ? matchScore : mismatchScore);
            int up = dp[i - 1][j] + gapPenalty;
            int left = dp[i][j - 1] + gapPenalty;
            dp[i][j] = max(0, max(diagonal, max(up, left)));

            // Record
            if (dp[i][j] == 0) {
                traceback[i][j] = '0';
            } else if (dp[i][j] == diagonal) {
                traceback[i][j] = 'd';
            } else if (dp[i][j] == up) {
                traceback[i][j] = 'u';
            } else {
                traceback[i][j] = 'l';
            }

            // Update the maximum score and its position.
            if (dp[i][j] > result->score) {
                result->score = dp[i][j];
                trace_i = i;
                trace_j = j;
            }
        }
    }



    // Traceback
    result->alignedPattern = "";
    result->alignedReference = "";
    vector<char> tracebacks;
    while (trace_i > 0 && trace_j > 0 && dp[trace_i][trace_j] != 0) {
        if (traceback[trace_i][trace_j] == 'd') {
            result->alignedPattern.push_back(patterns[trace_i - 1]);
            result->alignedReference.push_back(references[trace_j - 1]);
            tracebacks.push_back('M');
            --trace_i;
            --trace_j;
        } else if (traceback[trace_i][trace_j] == 'u') {
            result->alignedPattern.push_back(patterns[trace_i - 1]);
            result->alignedReference.push_back('-');
            tracebacks.push_back('I');
            --trace_i;
        } else if (traceback[trace_i][trace_j] == 'l') {
            result->alignedPattern.push_back('-');
            result->alignedReference.push_back(references[trace_j - 1]);
            tracebacks.push_back('D');
            --trace_j;
        }
    }

    reverse(result->alignedPattern.begin(), result->alignedPattern.end());
    reverse(result->alignedReference.begin(), result->alignedReference.end());

    result->cigar = prepareCigarString(tracebacks); // Prepare CIGAR string from traceback(reverse).
    result->mdz = prepareMDZString(result->alignedPattern, result->alignedReference, tracebacks);
    return result;
}

int overlapLongestExactMatch(const string& alignedPattern, const string& alignedReference) {
    int maxMatch = 0, current = 0;
    for (size_t i = 0; i < alignedPattern.size(); ++i) {
        if (alignedPattern[i] != '-' && alignedReference[i] != '-' && alignedPattern[i] == alignedReference[i]) {
            ++current;
            maxMatch = max(maxMatch, current);
        } else {
            current = 0;
        }
    }
    return maxMatch;
}

int main(int argc, char* argv[]) {
    if (argc < 9) {
        cerr << "Usage: " << argv[0] << " -g|-l -p <patterns.fasta> -t <texts.fasta> -o <output.txt> -s <match> <mismatch> <gap>" << endl;
        return 1;
    }

    bool global = false, local = false;
    string patternFile, referenceFile, outputFileName;
    int matchScore = 0, mismatchScore = 0, gapPenalty = 0;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-g") {
            global = true;
        } else if (arg == "-l") {
            local = true;
        } else if (arg == "-p" && i + 1 < argc) {
            patternFile = argv[++i];
        } else if (arg == "-t" && i + 1 < argc) {
            referenceFile = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            outputFileName = argv[++i];
        } else if (arg == "-s" && i + 3 < argc) {
            matchScore = atoi(argv[++i]);
            mismatchScore = atoi(argv[++i]);
            gapPenalty = atoi(argv[++i]);
        }
    }


    /*
    // TODO: FOR MANUAL TESTING
    bool global = true, local = false;
    string patternFile = "patterns.fasta", referenceFile = "texts.fasta", outputFileName = "output.txt";
    int matchScore = 1, mismatchScore = -1, gapPenalty = -1;
    */

    vector<string> patterns = readFasta(patternFile);
    vector<string> references = readFasta(referenceFile);
    if (patterns.size() != references.size()) {
        cerr << "Error: Number of patterns and references do not match." << endl;
        return 1;
    }

    vector<AlignmentResult*> results;
    // Score is used for local alignment, overlap is used as score for global alignment.
    int bestScore = -1000000, bestIndex = -1;

    for (size_t i = 0; i < patterns.size(); ++i) {
        AlignmentResult* result;

        if (global) {
            result = globalAlignmentNeedlemanWunsch(patterns[i], references[i], matchScore, mismatchScore, gapPenalty);
        } else {
            result = localAlignmentSmithWaterman(patterns[i], references[i], matchScore, mismatchScore, gapPenalty);
        }

        results.push_back(result);
    }

    AlignmentResult* bestResult = results.size() ? results[0] : nullptr;

    for (size_t i = 0; i < results.size(); ++i) {
        if (global) {
            int overlap = overlapLongestExactMatch(results[i]->alignedPattern, results[i]->alignedReference);
            if (overlap > bestScore) {
                bestScore = overlap;
                bestResult = results[i];
                bestIndex = i;
            }

        } else {
            if (results[i]->score > bestScore) {
                bestScore = results[i]->score;
                bestResult = results[i];
                bestIndex = i;
            }
        }

        // TODO: REMOVE THIS BLOCK FOR FINAL VERSION
        cout << "Pattern " << i + 1 << ":" << endl;
        cout << "pattern: " << patterns[i] << endl;
        cout << "reference: " << references[i] << endl;
        cout << "pattern size: " << patterns[i].size() << endl;
        cout << "reference size: " << references[i].size() << endl;
        cout << "cigar: " << results[i]->cigar << endl;
        cout << "aligned pat: " << results[i]->alignedPattern << endl;
        cout << "aligned ref: " << results[i]->alignedReference << endl;
        cout << "mdz: " << results[i]->mdz << endl;
        cout << "score: " << results[i]->score << endl;
        cout << endl;
    }

    ofstream outputFile(outputFileName.c_str());
    if (!outputFile) {
        cerr << "Error: Cannot open output file " << outputFileName << endl;
        return 1;
    }

    if (global && bestResult != nullptr && bestIndex >= 0) {
        outputFile << "Longest overlap:" << endl;
        outputFile << "pattern=" << patterns[bestIndex] << endl;
        outputFile << "reference=" << references[bestIndex] << endl;
        outputFile << "Score =" << bestResult->score << endl;
        outputFile << "CIGAR =" << bestResult->cigar << endl;
        outputFile << "MD:Z=" << bestResult->mdz << endl;
    } else if (local && bestResult != nullptr && bestIndex >= 0) {
        outputFile << "Highest local alignment score:" << endl;
        outputFile << "pattern=" << patterns[bestIndex] << endl;
        outputFile << "reference=" << references[bestIndex] << endl;
        outputFile << "Score =" << bestResult->score << endl;
        outputFile << "CIGAR =" << bestResult->cigar << endl;
        outputFile << "MD:Z=" << bestResult->mdz << endl;
    }

    for (size_t i = 0; i < results.size(); ++i) {
        if (results[i]) {
            delete results[i];
        }
    }

    outputFile.close();
    return 0;
}
