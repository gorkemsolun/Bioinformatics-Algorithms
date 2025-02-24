/*
 * main.cpp
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
 * It reports pattern occurrences per reference (using 1-indexed positions).
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
    SuffixTreeNode* suffixLink;
    int start;
    int* end; // For leaves, all share the same end pointer (leafEnd)
    int suffixIndex; // For leaves, stores starting index in the global text; for internal nodes, -1.

    SuffixTreeNode(int start, int* end)
        : suffixLink(nullptr), start(start), end(end), suffixIndex(-1) {
    }
};

class SuffixTree {
public:
    string text;       // The combined text (all references concatenated with unique terminators)
    SuffixTreeNode* root;
    SuffixTreeNode* lastNewNode;
    SuffixTreeNode* activeNode;
    int activeEdge;
    int activeLength;
    int remainingSuffixCount;
    int leafEnd;       // Global end for all leaves (updated in each phase)
    int* rootEnd;
    int* splitEnd;
    int size;          // Length of text

    SuffixTree(const string& text) : text(text), root(nullptr), lastNewNode(nullptr),
        activeEdge(-1), activeLength(0), remainingSuffixCount(0), leafEnd(-1), rootEnd(nullptr), splitEnd(nullptr) {
        size = text.size();
        // Create root node with start = -1 and end = nullptr.
        root = new SuffixTreeNode(-1, new int(-1));
        activeNode = root;
    }

    ~SuffixTree() {
        freeTree(root);
        if (rootEnd) delete rootEnd;
    }

    // Main build function: iteratively add each character.
    void build() {
        for (int i = 0; i < size; i++) {
            extendSuffixTree(i);
        }
    }

    // Ukkonen extension for phase pos.
    void extendSuffixTree(int pos) {
        leafEnd = pos;
        remainingSuffixCount++;
        lastNewNode = nullptr;

        while (remainingSuffixCount > 0) {
            if (activeLength == 0)
                activeEdge = pos; // current character index

            char curChar = text[activeEdge];
            // If there is no outgoing edge starting with curChar from activeNode
            if (activeNode->children.find(curChar) == activeNode->children.end()) {
                // Create a new leaf edge
                activeNode->children[curChar] = new SuffixTreeNode(pos, &leafEnd);

                // Add suffix link from last created internal node (if any)
                if (lastNewNode != nullptr) {
                    lastNewNode->suffixLink = activeNode;
                    lastNewNode = nullptr;
                }
            } else {
                // There is an outgoing edge starting with curChar.
                SuffixTreeNode* next = activeNode->children[curChar];
                int edgeLen = edgeLength(next);
                if (activeLength >= edgeLen) {
                    activeEdge += edgeLen;
                    activeLength -= edgeLen;
                    activeNode = next;
                    continue;
                }
                // Check if the current character on the edge is equal to text[pos]
                if (text[next->start + activeLength] == text[pos]) {
                    if (lastNewNode != nullptr && activeNode != root) {
                        lastNewNode->suffixLink = activeNode;
                        lastNewNode = nullptr;
                    }
                    activeLength++;
                    break;
                }
                // Need to split the edge
                splitEnd = new int;
                *splitEnd = next->start + activeLength - 1;
                SuffixTreeNode* split = new SuffixTreeNode(next->start, splitEnd);
                activeNode->children[curChar] = split;
                // New leaf for current character
                split->children[text[pos]] = new SuffixTreeNode(pos, &leafEnd);
                next->start += activeLength;
                split->children[text[next->start]] = next;

                if (lastNewNode != nullptr)
                    lastNewNode->suffixLink = split;
                lastNewNode = split;
            }
            remainingSuffixCount--;
            if (activeNode == root && activeLength > 0) {
                activeLength--;
                activeEdge = pos - remainingSuffixCount + 1;
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
            node->suffixIndex = size - labelHeight;
        }
    }

    // Search for occurrences of a pattern in the suffix tree.
    // Returns a vector of global starting indices (in text) where the pattern occurs.
    vector<int> search(const string& pattern) {
        vector<int> result;
        SuffixTreeNode* current = root;
        int i = 0;
        while (i < pattern.size()) {
            if (current->children.find(pattern[i]) != current->children.end()) {
                SuffixTreeNode* edge = current->children[pattern[i]];
                int edgeLen = edgeLength(edge);
                int j = 0;
                while (j < edgeLen && i < pattern.size() && text[edge->start + j] == pattern[i]) {
                    i++; j++;
                }
                if (j == edgeLen)
                    current = edge;
                else {
                    if (i == pattern.size())
                        collectLeafIndices(edge, result);
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

    // Print the suffix tree in DOT format (for visualization).
    void printDot(const string& filename) {
        ofstream dotFile(filename);
        dotFile << "digraph suffix_tree {\n";
        // Assign unique integer ids to nodes via DFS.
        map<SuffixTreeNode*, int> nodeId;
        int idCounter = 0;
        function<void(SuffixTreeNode*)> dfs = [&](SuffixTreeNode* node) {
            int curId = idCounter++;
            nodeId[node] = curId;
            string label = "";
            if (node->suffixIndex != -1)
                label = to_string(node->suffixIndex);
            dotFile << "node" << curId << " [label=\"" << label << "\"];\n";
            for (auto& childPair : node->children) {
                SuffixTreeNode* child = childPair.second;
                dfs(child);
                int childId = nodeId[child];
                string edgeLabel = text.substr(child->start, *(child->end) - child->start + 1);
                dotFile << "node" << curId << " -> node" << childId << " [label=\"" << edgeLabel << "\"];\n";
            }
            };
        dfs(root);
        dotFile << "}\n";
        dotFile.close();
    }

private:
    // Helper: recursively free nodes.
    void freeTree(SuffixTreeNode* node) {
        if (node == nullptr) return;
        for (auto& childPair : node->children) {
            freeTree(childPair.second);
        }
        // Only delete end pointer if it was allocated (i.e. not pointing to leafEnd)
        if (node->end != &leafEnd)
            delete node->end;
        delete node;
    }

    // Helper: collect suffix indices from all leaves in the subtree rooted at node.
    void collectLeafIndices(SuffixTreeNode* node, vector<int>& result) {
        if (node == nullptr) return;
        if (node->suffixIndex != -1) {
            result.push_back(node->suffixIndex);
            return;
        }
        for (auto& childPair : node->children)
            collectLeafIndices(childPair.second, result);
    }
};

//
// Main function: parses command-line arguments, reads input FASTA files, builds the suffix tree,
// searches for each pattern, and writes output files.
//
int main(int argc, char** argv) {
    string refFile, patFile, outPrefix;
    bool dotFlag = false;

    // Simple argument parsing
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "-r" && i + 1 < argc) {
            refFile = argv[++i];
        } else if (arg == "-p" && i + 1 < argc) {
            patFile = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            outPrefix = argv[++i];
        } else if (arg == "-d") {
            dotFlag = true;
        }
    }
    if (refFile.empty() || patFile.empty() || outPrefix.empty()) {
        cerr << "Usage: " << argv[0] << " -r <reference.fasta> -p <patterns.fasta> -o <output_prefix> [-d]" << endl;
        return 1;
    }

    // Read the reference FASTA file (multi-line sequences)
    vector<pair<string, string>> references; // pair: header, sequence
    ifstream refIn(refFile);
    if (!refIn) {
        cerr << "Cannot open reference file: " << refFile << endl;
        return 1;
    }
    string line, header, seq;
    while (getline(refIn, line)) {
        if (line.empty())
            continue;
        if (line[0] == '>') {
            if (!header.empty()) {
                references.push_back({ header, seq });
                seq = "";
            }
            header = line.substr(1); // drop the '>' character
        } else {
            seq += line;
        }
    }
    if (!header.empty())
        references.push_back({ header, seq });
    refIn.close();

    // Read the patterns FASTA file (assume header then pattern on next line)
    vector<pair<string, string>> patterns;
    ifstream patIn(patFile);
    if (!patIn) {
        cerr << "Cannot open patterns file: " << patFile << endl;
        return 1;
    }
    header = "";
    string pat;
    bool headerFlag = false;
    while (getline(patIn, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            header = line.substr(1);
            headerFlag = true;
        } else if (headerFlag) {
            pat = line;
            patterns.push_back({ header, pat });
            headerFlag = false;
        }
    }
    patIn.close();

    // Build a combined text from all reference sequences.
    // Use a unique terminator for each reference (alphabet is {A,C,G,T} so any other char is safe).
    string combinedText = "";
    vector<tuple<int, int, string>> refBoundaries; // (start index, end index, header)
    vector<int> refStartIndices;
    int numRefs = references.size();
    vector<char> terminators;
    // Use a set of available unique terminators.
    string available = "$#@%^&!";
    if (numRefs > available.size()) {
        // If more references than available symbols, generate symbols from printable ASCII excluding A, C, G, T.
        for (char c = 33; c < 127 && terminators.size() < numRefs; c++) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                terminators.push_back(c);
        }
    } else {
        for (int i = 0; i < numRefs; i++) {
            terminators.push_back(available[i]);
        }
    }
    for (int i = 0; i < numRefs; i++) {
        int startIndex = combinedText.size();
        combinedText += references[i].second;
        int seqLen = references[i].second.size();
        // Append a unique terminator.
        combinedText.push_back(terminators[i]);
        int endIndex = combinedText.size() - 1; // includes terminator
        refBoundaries.push_back(make_tuple(startIndex, endIndex, references[i].first));
        refStartIndices.push_back(startIndex);
    }

    // Build the suffix tree from the combined text.
    SuffixTree st(combinedText);
    st.build();
    st.setSuffixIndexByDFS(st.root, 0);

    // A helper lambda to map a global index in combinedText to a reference header and a 1-indexed position.
    auto mapIndexToRef = [&](int idx) -> pair<string, int> {
        for (int i = 0; i < refBoundaries.size(); i++) {
            int start, end;
            string refHeader;
            tie(start, end, refHeader) = refBoundaries[i];
            int refLen = end - start; // terminator is at the end
            if (idx >= start && idx < start + refLen) {
                return { refHeader, idx - start + 1 };
            }
        }
        return { "", -1 };
        };

    // Open the output file for writing pattern occurrences.
    ofstream outFile(outPrefix + ".txt");
    if (!outFile) {
        cerr << "Cannot open output file: " << outPrefix + ".txt" << endl;
        return 1;
    }

    // For each pattern, search in the suffix tree and group the results by reference.
    for (auto& p : patterns) {
        string patHeader = p.first;
        string patSeq = p.second;
        vector<int> occ = st.search(patSeq);
        // Group occurrences by reference header.
        map<string, vector<int>> occMap;
        for (int pos : occ) {
            auto info = mapIndexToRef(pos);
            if (info.first != "")
                occMap[info.first].push_back(info.second);
        }
        for (auto& entry : occMap)
            sort(entry.second.begin(), entry.second.end());

        outFile << "(" << patHeader << ") - ";
        bool firstRef = true;
        for (auto& entry : occMap) {
            if (!firstRef)
                outFile << ", ";
            firstRef = false;
            outFile << entry.first << ":";
            for (size_t i = 0; i < entry.second.size(); i++) {
                if (i > 0)
                    outFile << ",";
                outFile << entry.second[i];
            }
        }
        outFile << "\n";
    }
    outFile.close();

    // If the -d switch was provided, output the DOT file for the suffix tree.
    if (dotFlag) {
        st.printDot(outPrefix + ".dot");
    }

    return 0;
}
