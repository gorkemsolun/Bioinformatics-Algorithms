/*
 * Görkem Kadir Solun
 *
 * To prepare the executable: make
 * To clean the executable: make clean
 * To run the executable: ./hw1 -r <reference_fasta> -p <pattern_fasta> -o <output_prefix>
 *
 *
 * A complete solution for the multiple pattern matching assignment.
 * This implementation builds a single suffix tree that holds all reference sequences
 * (each separated by a unique terminator) and then searches for each pattern.
 *
 * Command-line arguments:
 *   -r <reference_fasta> : FASTA file with reference sequences.
 *   -p <pattern_fasta>   : FASTA file with patterns.
 *   -o <output_prefix>   : The prefix for output files (a .txt file is always produced,
 *                          and if -d is specified, a DOT file is produced).
 *   -d                   : Optional; if provided, outputs the suffix tree in DOT format.
 *
 * The implementation uses Ukkonen's algorithm for suffix tree construction.
 * It reports pattern occurrences per reference (using 0-indexed positions).
 *
 * The node leaf naming logic is slightly modified as the following to be suitable for single suffix tree
 * <reference_name>:<1-based_start_index_in_that_reference>:<global_concat_index>
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <functional>
#include <tuple>
using namespace std;

struct SuffixTreeNode {
    // Each edge is stored in a map from character to child node.
    map<char, SuffixTreeNode*> children;
    SuffixTreeNode* suffixLink; // Suffix link for Ukkonen's algorithm (not used in search)
    int start, * end; // For leaves, all share the same end pointer (leafEnd)
    int suffixIndex; // For leaves, stores starting index in the global text; for internal nodes, -1.

    SuffixTreeNode(int start, int* end)
        : suffixLink(nullptr), start(start), end(end), suffixIndex(-1) {
    }
};

class SuffixTree {
public:
    string text;       // The combined text (all references concatenated with unique terminators)
    vector<tuple<int, int, string>> referenceBoundaries; // each tuple = (startIndexInConcatenated, endIndexInConcatenated, referenceName)

    SuffixTreeNode* root;
    SuffixTreeNode* lastNewNode;
    SuffixTreeNode* activeNode;
    int activeEdge; // Index of the current character on the active edge
    int activeLength; // Length of the active prefix
    int remainingSuffixCount; // Number of suffixes yet to be added in the current phase
    int leafEnd; // Global end for all leaves (updated in each phase)
    int* rootEnd; // Pointer to the end of the root's edge (nullptr) 
    int* splitEnd; // Pointer to the end of the split edge (updated in each phase)

    SuffixTree(const string& text, const vector<tuple<int, int, string>>& referenceBoundaries) : text(text), referenceBoundaries(referenceBoundaries), lastNewNode(nullptr),
        activeEdge(-1), activeLength(0), remainingSuffixCount(0), leafEnd(-1), rootEnd(nullptr), splitEnd(nullptr) {
        root = new SuffixTreeNode(-1, new int(-1)); // Create root node with start = -1 and end = nullptr.
        activeNode = root;
    }

    ~SuffixTree() {
        freeTree(root);
        if (rootEnd) { delete rootEnd; }
    }

    // Iteratively add each character
    void build() {
        for (int i = 0; i < (int) text.size(); ++i) {
            extendSuffixTree(i);
        }
    }

    // Ukkonen extension for phase position
    void extendSuffixTree(int position) {
        leafEnd = position;
        ++remainingSuffixCount;
        lastNewNode = nullptr;

        while (remainingSuffixCount > 0) {
            if (activeLength == 0) {
                activeEdge = position; // current character index
            }

            // Current character is activeEdge-th character on the activeNode
            // If there is no outgoing edge starting with text[activeEdge] from activeNode
            if (activeNode->children.find(text[activeEdge]) == activeNode->children.end()) {
                // Create a new leaf edge
                activeNode->children[text[activeEdge]] = new SuffixTreeNode(position, &leafEnd);

                // Add suffix link from last created internal node (if any)
                if (lastNewNode != nullptr) {
                    lastNewNode->suffixLink = activeNode;
                    lastNewNode = nullptr;
                }
            } else {
                // There is an outgoing edge starting with text[activeEdge].
                SuffixTreeNode* next = activeNode->children[text[activeEdge]];
                int edgeLen = edgeLength(next);
                if (activeLength >= edgeLen) {
                    activeEdge += edgeLen;
                    activeLength -= edgeLen;
                    activeNode = next;
                    continue;
                }
                // Check if the current character on the edge is equal to text[position]
                if (text[next->start + activeLength] == text[position]) {
                    if (lastNewNode != nullptr && activeNode != root) {
                        lastNewNode->suffixLink = activeNode;
                        lastNewNode = nullptr;
                    }
                    ++activeLength;
                    break;
                }
                // Need to split the edge
                splitEnd = new int(next->start + activeLength - 1);
                SuffixTreeNode* split = new SuffixTreeNode(next->start, splitEnd);
                activeNode->children[text[activeEdge]] = split;
                // New leaf for current character
                split->children[text[position]] = new SuffixTreeNode(position, &leafEnd);
                next->start += activeLength;
                split->children[text[next->start]] = next;

                if (lastNewNode != nullptr) {
                    lastNewNode->suffixLink = split;
                }
                lastNewNode = split;
            }
            --remainingSuffixCount;
            if (activeNode == root && activeLength > 0) {
                --activeLength;
                activeEdge = position - remainingSuffixCount + 1;
            } else if (activeNode != root) {
                activeNode = (activeNode->suffixLink != nullptr) ? activeNode->suffixLink : root;
            }
        }
    }

    // Return the length of the edge (for a given node)
    int edgeLength(SuffixTreeNode* node) {
        return *(node->end) - node->start + 1;
    }

    // After tree construction, set suffixIndex for leaves using DFS.
    void setSuffixIndexByDFS(SuffixTreeNode* node, int labelHeight) {
        if (node == nullptr) return;
        bool leaf = true;
        for (auto& childPair : node->children) {
            leaf = false;
            setSuffixIndexByDFS(childPair.second, labelHeight + edgeLength(childPair.second));
        }
        if (leaf) {
            node->suffixIndex = text.size() - labelHeight;
        }
    }

    // Search for occurrences of a pattern in the suffix tree.
    // Returns a vector of global starting indices (in text) where the pattern occurs.
    vector<int> search(const string& pattern) {
        vector<int> result;
        SuffixTreeNode* current = root;
        int i = 0;
        while (i < (int) pattern.size()) {
            if (current->children.find(pattern[i]) != current->children.end()) {
                SuffixTreeNode* edge = current->children[pattern[i]];
                int edgeLen = edgeLength(edge);
                int j = 0;
                while (j < edgeLen && i < (int) pattern.size() && text[edge->start + j] == pattern[i]) {
                    ++i; ++j;
                }
                if (j == edgeLen) {
                    current = edge;
                } else {
                    if (i == (int) pattern.size()) {
                        collectLeafIndices(edge, result);
                    }
                    return result; // mismatch: pattern not found
                }
            } else {
                return result; // pattern not found
            }
        }
        // Pattern completely matched; now collect all leaf indices below current node.
        collectLeafIndices(current, result);
        return result;
    }

    // Helper that maps a global index to(referenceName, localPosition, globalIndex)
    // If not found, returns ("", -1, -1).
    tuple<string, int, int> mapIndexToRef(int globalIdx) const {
        for (auto& rb : referenceBoundaries) {
            int start = get<0>(rb), end = get<1>(rb);
            string name = get<2>(rb);
            // end - start is the length (minus the terminator)
            // But be consistent with how you stored boundaries.
            if (globalIdx >= start && globalIdx < end) {
                return make_tuple(name, globalIdx - start, globalIdx);
            }
        }
        return make_tuple(string(""), -1, -1);
    }

    // Print the suffix tree in DOT format (for visualization).
    void printDot(const string& filename) {
        ofstream dotFile(filename);
        dotFile << "digraph suffix_tree {" << endl;
        // Assign unique integer ids to nodes via DFS.
        map<SuffixTreeNode*, int> nodeId;
        int idCounter = 0;
        function<void(SuffixTreeNode*)> dfs = [&](SuffixTreeNode* node) {
            int curId = idCounter++;
            nodeId[node] = curId;

            // Build label for this node
            string label = "";
            if (node->suffixIndex != -1) {
                // It's a leaf
                string referenceName;
                int localPosition, globalPosition;
                tie(referenceName, localPosition, globalPosition) = mapIndexToRef(node->suffixIndex);
                if (!referenceName.empty() && localPosition != -1) {
                    // e.g. "chr2:5:10"
                    label = referenceName + ":" + to_string(localPosition) + ":" + to_string(globalPosition);
                } else {
                    // fallback if not found
                    label = to_string(node->suffixIndex);
                }
            }
            dotFile << "node" << curId << " [label=\"" << label << "\"];" << endl;

            // Recur for children
            for (auto& ch : node->children) {
                SuffixTreeNode* child = ch.second;
                dfs(child);
                int childId = nodeId[child];
                // Edge label is the substring from child->start to *(child->end)
                string edgeLabel = text.substr(child->start, *(child->end) - child->start + 1);
                dotFile << "node" << curId << " -> node" << childId
                    << " [label=\"" << edgeLabel << "\"];" << endl;
            }
            };

        dfs(root);
        dotFile << "}" << endl;
        dotFile.close();
    }

private:
    // Helper: recursively free nodes.
    void freeTree(SuffixTreeNode* node) {
        if (node == nullptr) { return; }
        for (auto& childPair : node->children) {
            freeTree(childPair.second);
        }
        // Only delete end pointer if it was allocated (i.e. not pointing to leafEnd)
        if (node->end != &leafEnd) {
            delete node->end;
        }
        delete node;
    }

    // Helper: collect suffix indices from all leaves in the subtree rooted at node.
    void collectLeafIndices(SuffixTreeNode* node, vector<int>& result) {
        if (node == nullptr) { return; }
        if (node->suffixIndex != -1) {
            result.push_back(node->suffixIndex);
            return;
        }
        for (auto& childPair : node->children) {
            collectLeafIndices(childPair.second, result);
        }
    }
};

vector<pair<string, string>>* readSequences(const string& filename) {
    vector<pair<string, string>>* sequences = new vector<pair<string, string>>();
    ifstream file(filename);
    if (!file) {
        cerr << "Cannot open file: " << filename << endl;
        return sequences;
    }
    string line, header, seq;
    while (getline(file, line)) {
        if (line.empty()) {
            break;
        }

        // Trim leading and trailing whitespaces
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        line.erase(0, line.find_first_not_of(" \t\r\n"));

        if (line[0] == '>') {
            if (!header.empty()) {
                sequences->push_back({ header, seq });
                seq = "";
            }
            header = line.substr(1); // drop the '>' character
        } else {
            seq += line;
        }
    }
    if (!header.empty()) {
        sequences->push_back({ header, seq });
    }
    file.close();
    return sequences;
}

//
// Main function: parses command-line arguments, reads input FASTA files, builds the suffix tree,
// searches for each pattern, and writes output files.
int main(int argc, char** argv) {
    string referenceFile, patternFile, outputPrefix;
    bool dotFlag = false;

    // Simple argument parsing
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-r" && i + 1 < argc) {
            referenceFile = argv[++i];
        } else if (arg == "-p" && i + 1 < argc) {
            patternFile = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            outputPrefix = argv[++i];
        } else if (arg == "-d") {
            dotFlag = true;
        }
    }
    if (referenceFile.empty() || patternFile.empty() || outputPrefix.empty()) {
        cerr << "Usage: " << argv[0] << " -r <reference.fasta> -p <patterns.fasta> -o <output_prefix> [-d]" << endl;
        return 1;
    }

    // Read the reference FASTA file
    vector<pair<string, string>>* references = readSequences(referenceFile); // pair: header, sequence
    // Read the patterns FASTA file
    vector<pair<string, string>>* patterns = readSequences(patternFile);

    // Build a combined text from all reference sequences.
    // Use a unique terminator for each reference (alphabet is {A,C,G,T} so any other char is safe).
    string combinedText = "";
    vector<tuple<int, int, string>> referenceBoundaries; // (start index, end index, header)
    vector<int> referenceStartIndices;
    vector<char> terminators;
    // Use a set of available unique terminators.
    string available = "$#@%^&!";
    if (references->size() > available.size()) {
        // If more references than available symbols, generate symbols from printable ASCII excluding A, C, G, T.
        for (char c = 33; c < 127 && terminators.size() < references->size(); ++c) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                terminators.push_back(c);
            }
        }
    } else {
        for (int i = 0; i < (int) references->size(); ++i) {
            terminators.push_back(available[i]);
        }
    }

    for (int i = 0; i < (int) references->size(); ++i) {
        int startIndex = combinedText.size();
        combinedText += (*references)[i].second;
        combinedText.push_back(terminators[i]);// Append a unique terminator.
        int endIndex = combinedText.size() - 1; // includes terminator
        referenceBoundaries.push_back(make_tuple(startIndex, endIndex, (*references)[i].first));
        referenceStartIndices.push_back(startIndex);
    }

    // Build the suffix tree from the combined text.
    SuffixTree st(combinedText, referenceBoundaries);
    st.build();
    st.setSuffixIndexByDFS(st.root, 0); // Set suffix indices for leaves using DFS.

    // A helper lambda to map a global index in combinedText to a reference header and a position.
    auto mapIndexToRef = [&](int idx) -> pair<string, int> {
        for (int i = 0; i < (int) referenceBoundaries.size(); ++i) {
            int start, end;
            string refHeader;
            tie(start, end, refHeader) = referenceBoundaries[i];
            // reference length = end - start, terminator is at the end
            if (idx >= start && idx < end) {
                return { refHeader, idx - start };
            }
        }
        return { "", -1 };
        };

    // Open the output file for writing pattern occurrences.
    ofstream outputFile(outputPrefix + ".txt");
    if (!outputFile) {
        cerr << "Cannot open output file: " << outputPrefix + ".txt" << endl;
        return 1;
    }

    // For each pattern, search in the suffix tree and group the results by reference.
    for (auto& p : (*patterns)) {
        vector<int> occurrences = st.search(p.second);
        // Group occurrences by reference header.
        map<string, vector<int>> occurrenceMap;
        for (int position : occurrences) {
            auto info = mapIndexToRef(position);
            if (info.first != "") {
                occurrenceMap[info.first].push_back(info.second);
            }
        }
        for (auto& entry : occurrenceMap) {
            sort(entry.second.begin(), entry.second.end());
        }

        outputFile << "(" << p.first << ") - ";
        bool isFirstReference = true;
        for (auto& entry : occurrenceMap) {
            if (!isFirstReference) {
                outputFile << ", ";
            }
            isFirstReference = false;
            outputFile << entry.first << ":";
            for (size_t i = 0; i < entry.second.size(); ++i) {
                if (i > 0) {
                    outputFile << ",";
                }
                outputFile << entry.second[i];
            }
        }
        outputFile << endl;
    }
    outputFile.close();

    // If the -d switch was provided, output the DOT file for the suffix tree.
    if (dotFlag) {
        st.printDot(outputPrefix + ".dot");
    }

    delete references;
    delete patterns;

    return 0;
}
