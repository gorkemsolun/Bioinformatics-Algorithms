// Görkem Kadir Solun 22003214
// UPGMA to construct alignment1 phylogenetic relationship. Using the Needleman - Wunsch algorithm, 
// you will build alignment1 distance matrix for the DNA sequences of different species in alignment1 FASTA file.


#include <bits/stdc++.h>

using namespace std;

struct Cluster {
    int size;
    double height;
    string newick;
};

pair<string*, string*> needleman_wunsch(const string& sequence1, const string& sequence2,
    int matchScore, int mismatchScore, int gapScore) {
    vector<vector<int>> dpTable(sequence1.size() + 1, vector<int>(sequence2.size() + 1));
    vector<vector<char>> traceback(sequence1.size() + 1, vector<char>(sequence2.size() + 1));

    for (size_t i = 1; i <= sequence1.size(); ++i) {
        dpTable[i][0] = dpTable[i - 1][0] + gapScore;
        traceback[i][0] = 'U';
    }
    for (size_t j = 1; j <= sequence2.size(); ++j) {
        dpTable[0][j] = dpTable[0][j - 1] + gapScore;
        traceback[0][j] = 'L';
    }

    for (size_t i = 1; i <= sequence1.size(); ++i) {
        for (size_t j = 1; j <= sequence2.size(); ++j) {
            int scoreUp = dpTable[i - 1][j] + gapScore;
            int scoreLeft = dpTable[i][j - 1] + gapScore;

            // diagonal score
            dpTable[i][j] = dpTable[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchScore);
            traceback[i][j] = 'D';

            if (scoreUp > dpTable[i][j]) {
                dpTable[i][j] = scoreUp;
                traceback[i][j] = 'U';
            }
            if (scoreLeft > dpTable[i][j]) {
                dpTable[i][j] = scoreLeft;
                traceback[i][j] = 'L';
            }
        }
    }

    int i = sequence1.size(), j = sequence2.size();
    string* alignment1 = new string(), * alignment2 = new string();
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && traceback[i][j] == 'D') {
            alignment1->push_back(sequence1[i - 1]);
            alignment2->push_back(sequence2[j - 1]);
            --i; --j;
        } else if (i > 0 && traceback[i][j] == 'U') {
            alignment1->push_back(sequence1[i - 1]);
            alignment2->push_back('-');
            --i;
        } else {
            alignment1->push_back('-');
            alignment2->push_back(sequence2[j - 1]);
            --j;
        }
    }

    reverse(alignment1->begin(), alignment1->end());
    reverse(alignment2->begin(), alignment2->end());

    return { alignment1, alignment2 };
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        cerr << "Usage: " << argv[0] << " -i <input.fasta> -t <tree.txt> -s <match> <mismatch> <gap>\n";
        return 1;
    }

    string inputFile, outputFile;
    int matchScore = 1, mismatchScore = -1, gapScore = -1;

    for (int i = 1; i < argc; ++i) {
        string opt = argv[i];
        if (opt == "-i" && i + 1 < argc) {
            inputFile = argv[++i];
        } else if (opt == "-t" && i + 1 < argc) {
            outputFile = argv[++i];
        } else if (opt == "-s" && i + 3 < argc) {
            matchScore = stoi(argv[++i]);
            mismatchScore = stoi(argv[++i]);
            gapScore = stoi(argv[++i]);
        } else {
            cerr << "Unknown option: " << opt << '\n';
            return 1;
        }
    }


    ifstream input_file_stream(inputFile);
    if (!input_file_stream) {
        cerr << "Error opening input file: " << inputFile << '\n';
        return 1;
    }

    vector<pair<string, string>> sequences;
    string line, currentId, currentSequence;

    while (getline(input_file_stream, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            if (!currentId.empty()) {
                sequences.emplace_back(currentId, currentSequence);
            }

            currentId = line.substr(1);
            size_t position = currentId.find_last_of("\n");
            if (position != string::npos) {
                currentId = currentId.substr(0, position);
            }

            currentSequence.clear();
        } else {
            currentSequence += line;
        }
    }
    if (!currentId.empty()) {
        sequences.emplace_back(currentId, currentSequence);
    }

    input_file_stream.close();

    vector<vector<double>> distances(sequences.size(), vector<double>(sequences.size(), 0.0));
    for (size_t i = 0; i < sequences.size(); ++i) {
        for (size_t j = i + 1; j < sequences.size(); ++j) {
            pair<string*, string*> alignmentResult = needleman_wunsch(sequences[i].second, sequences[j].second,
                matchScore, mismatchScore, gapScore);

            // distance = number of mismatches + number of gaps

            for (size_t k = 0; k < alignmentResult.first->size(); ++k) {
                if ((*alignmentResult.first)[k] == '-' || (*alignmentResult.second)[k] == '-') {
                    ++distances[j][i];
                } else if ((*alignmentResult.first)[k] != (*alignmentResult.second)[k]) {
                    ++distances[j][i];
                }
            }

            delete alignmentResult.first;
            delete alignmentResult.second;

            distances[i][j] = distances[j][i];
        }
    }


    vector<Cluster> clusters(sequences.size());
    for (size_t i = 0; i < sequences.size(); ++i) {
        clusters[i].size = 1;
        clusters[i].height = 0.0;
        clusters[i].newick = sequences[i].first;
    }

    // UPGMA
    while (clusters.size() > 1) {
        double minDistance = numeric_limits<double>::infinity();

        size_t iMin = 0, jMin = 0;

        for (size_t i = 0; i < clusters.size(); ++i) {
            for (size_t j = i + 1; j < clusters.size(); ++j) {
                if (distances[i][j] < minDistance) {
                    minDistance = distances[i][j];
                    iMin = i;
                    jMin = j;
                }
            }
        }

        Cluster merged;
        merged.size = clusters[iMin].size + clusters[jMin].size;
        merged.height = minDistance / 2.0;

        merged.newick = "(" + clusters[iMin].newick + ":" + to_string(abs(merged.height - clusters[iMin].height))
            + "," + clusters[jMin].newick + ":" + to_string(abs(merged.height - clusters[jMin].height)) + ")";
        cout << merged.newick << endl;


        vector<Cluster> newClusters;
        vector<vector<double>> newDistances;
        for (size_t i = 0; i < clusters.size(); ++i) {
            if (i == iMin || i == jMin) {
                continue;
            }
            newClusters.push_back(clusters[i]);
        }

        newClusters.push_back(merged);

        newDistances.assign(newClusters.size(), vector<double>(newClusters.size(), 0.0));

        size_t index = 0;
        for (size_t i = 0; i < clusters.size(); ++i) {
            if (i == iMin || i == jMin) {
                continue;
            }
            size_t index2 = 0;
            for (size_t j = 0; j < clusters.size(); ++j) {
                if (j == iMin || j == jMin) {
                    continue;
                }
                newDistances[index][index2] = distances[i][j];
                ++index2;
            }

            newDistances[index][newClusters.size() - 1] = newDistances[newClusters.size() - 1][index] = (distances[iMin][i] * clusters[iMin].size + distances[jMin][i] * clusters[jMin].size) / merged.size;
            ++index;
        }

        clusters.swap(newClusters);
        distances.swap(newDistances);
    }

    string tree = clusters[0].newick + ":0.0;";

    ofstream outputFileStream(outputFile);
    if (!outputFileStream) {
        cerr << "Error opening output file: " << outputFile << '\n';
        return 1;
    }

    outputFileStream << tree << endl;
    cout << tree;
    cout << clusters.size() << endl;
    outputFileStream.close();


    return 0;
}
