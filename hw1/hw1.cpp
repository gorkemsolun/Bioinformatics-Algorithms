#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <cstdlib>
#include <string>
#include <algorithm>

using namespace std;

// Node for the Aho–Corasick trie.
struct ACNode {
    unordered_map<char, ACNode*> children;
    ACNode* failure;
    vector<int> patterns_ending_here;
    ACNode() : failure(nullptr) {}
};

class AhoCorasick {
public:
    ACNode* root;
    vector<string> patterns; // store patterns by index

    AhoCorasick() {
        root = new ACNode();
    }

    // Insert a pattern into the trie.
    void addPattern(const string& pattern, int index) {
        ACNode* node = root;
        for (char c : pattern) {
            if (node->children.find(c) == node->children.end())
                node->children[c] = new ACNode();
            node = node->children[c];
        }
        node->patterns_ending_here.push_back(index);
    }

    // Build failure links using BFS.
    void buildAutomaton() {
        queue<ACNode*> q;
        // For each child of the root, set failure link to root.
        for (auto& parent : root->children) {
            parent.second->failure = root;
            q.push(parent.second);
        }
        while (!q.empty()) {
            ACNode* curr = q.front();
            q.pop();
            for (auto& parent : curr->children) {
                char c = parent.first;
                ACNode* child = parent.second;
                ACNode* failure = curr->failure;
                while (failure && failure->children.find(c) == failure->children.end()) {
                    failure = failure->failure;
                }
                child->failure = (failure ? failure->children[c] : root);
                // Propagate patterns ending here from failure link.
                child->patterns_ending_here.insert(child->patterns_ending_here.end(),
                    child->failure->patterns_ending_here.begin(),
                    child->failure->patterns_ending_here.end());
                q.push(child);
            }
        }
    }

    // Search the given text (which in our case is the global concatenated references)
    // Returns, for each pattern, a vector of global (0-indexed) match positions.
    vector<vector<int>> search(const string& text) {
        vector<vector<int>> matches(patterns.size());
        ACNode* node = root;
        for (int i = 0; i < text.size(); i++) {
            char c = text[i];
            // Follow failure links if there is no edge for c.
            while (node && node->children.find(c) == node->children.end()) {
                node = node->failure;
            }
            if (!node) {
                node = root;
                continue;
            }
            node = node->children[c];
            for (int patternIndex : node->patterns_ending_here) {
                // Report the starting global position.
                int pos = i - patterns[patternIndex].size() + 1;
                matches[patternIndex].push_back(pos);
            }
        }
        return matches;
    }

    // Free nodes recursively.
    void freeNode(ACNode* node) {
        for (auto& parent : node->children) {
            freeNode(parent.second);
        }
        delete node;
    }

    ~AhoCorasick() {
        freeNode(root);
    }
};

struct FastaSequence {
    string name;
    string sequence;
};

vector<FastaSequence> readFastaSequences(const string& filename) {
    vector<FastaSequence> refs;
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening reference file: " << filename << endl;
        exit(1);
    }

    string line;
    FastaSequence current;
    bool inSequence = false;
    while (getline(infile, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            if (inSequence) {
                refs.push_back(current);
                current = FastaSequence();
            }
            current.name = line.substr(1);
            current.sequence = "";
            inSequence = true;
        } else {
            current.sequence += line;
        }
    }

    if (inSequence) {
        refs.push_back(current);
    }

    infile.close();
    return refs;
}

// This structure records where each reference’s sequence appears in the global text.
struct ReferenceBoundary {
    string name;
    int start; // global starting index (0-indexed) of the sequence
    int end;   // global ending index (inclusive) of the sequence (i.e. before the delimiter)
};

// Each node in the suffix tree.
struct STNode {
    // Children map: key is the first character of the edge.
    unordered_map<char, STNode*> children;
    // The edge label is represented by [start, *end] indices into the text.
    int start, * end;
    // For leaves, store the starting index (global position) of the corresponding suffix.
    vector<int> suffixIndices;
    // Node ID (for DOT output).
    int id;
    STNode(int start, int* end, int id) : start(start), end(end), id(id) {}
};

int globalNodeId = 0;

class SuffixTree {
public:
    string text;
    STNode* root;

    SuffixTree(const string& text) : text(text) {
        root = new STNode(-1, new int(-1), globalNodeId++);
        build();
    }

    ~SuffixTree() {
        freeNode(root);
    }

    // Naively insert every suffix.
    void build() {
        int n = text.size();
        for (int i = 0; i < n; i++) {
            insertSuffix(i);
        }
        // Compress the tree by merging nonbranching nodes.
        compress(root);
    }

    void insertSuffix(int index) {
        STNode* current = root;
        int n = text.size();
        int i = index;
        while (i < n) {
            char c = text[i];
            if (current->children.find(c) == current->children.end()) {
                int* leafEnd = new int(n - 1);
                STNode* leaf = new STNode(i, leafEnd, globalNodeId++);
                leaf->suffixIndices.push_back(index);
                current->children[c] = leaf;
                return;
            } else {
                STNode* child = current->children[c];
                int edgeStart = child->start;
                int edgeEnd = *(child->end);
                int j = edgeStart;
                while (j <= edgeEnd && i < n && text[j] == text[i]) {
                    j++; i++;
                }
                if (j > edgeEnd) {
                    current = child;
                    continue;
                }
                // Mismatch: split the edge.
                int* splitEnd = new int(j - 1);
                STNode* split = new STNode(child->start, splitEnd, globalNodeId++);
                current->children[c] = split;
                child->start = j;
                split->children[text[j]] = child;
                int* leafEnd = new int(n - 1);
                STNode* leaf = new STNode(i, leafEnd, globalNodeId++);
                leaf->suffixIndices.push_back(index);
                split->children[text[i]] = leaf;
                return;
            }
        }
    }

    // Recursively compress the tree by merging nodes with a single child.
    void compress(STNode* node) {
        for (auto& childPair : node->children) {
            STNode* child = childPair.second;
            compress(child);
            // Merge if the child is non–leaf and has exactly one child.
            while (child->children.size() == 1) {
                auto it = child->children.begin();
                STNode* grandchild = it->second;
                // Extend child's edge label to include the grandchild’s edge.
                child->end = grandchild->end;
                // Merge children.
                child->children = grandchild->children;
                child->suffixIndices.insert(child->suffixIndices.end(),
                    grandchild->suffixIndices.begin(),
                    grandchild->suffixIndices.end());
            }
        }
    }

    // Output the suffix tree in DOT language.
    void outputDot(const string& filename) {
        ofstream out(filename);
        out << "digraph suffix_tree {" << endl;
        outputNodeDot(root, out);
        out << "}" << endl;
        out.close();
    }

    void outputNodeDot(STNode* node, ofstream& out) {
        if (node == root) {
            out << "node_" << node->id << " [label=\"\"];" << endl;
        } else {
            string label = text.substr(node->start, *(node->end) - node->start + 1);
            out << "node_" << node->id << " [label=\"" << label << "\"];" << endl;
        }
        for (auto& childPair : node->children) {
            STNode* child = childPair.second;
            string edgeLabel = text.substr(child->start, *(child->end) - child->start + 1);
            out << "node_" << node->id << " -> node_" << child->id
                << " [label=\"" << edgeLabel << "\"];" << endl;
            outputNodeDot(child, out);
        }
    }

    // Free nodes recursively.
    void freeNode(STNode* node) {
        for (auto& child : node->children) {
            freeNode(child.second);
        }
        if (node->end) { delete node->end; }
        delete node;
    }
};

// Given a global position and pattern length, determine in which reference (if any)
// the entire pattern lies. Returns a pair of (reference name, local position [1-indexed]).
pair<string, int> getReferenceForPosition(int pos, int patLen, const vector<ReferenceBoundary>& boundaries) {
    for (const auto& b : boundaries) {
        if (pos >= b.start && pos + patLen - 1 <= b.end) {
            int localPos = pos - b.start + 1;
            return make_pair(b.name, localPos);
        }
    }
    return make_pair("", -1);
}

int main(int argc, char* argv[]) {
    string refFile, patFile, outPrefix;
    bool dotFlag = false;

    // Simple command–line argument parsing.
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
        cerr << "Usage: " << argv[0]
            << " -r <reference file> -p <pattern file> -o <output prefix> [-d]" << endl;
        return 1;
    }

    // Read references and patterns.
    vector<FastaSequence> refs = readFastaSequences(refFile);
    vector<FastaSequence> pats = readFastaSequences(patFile);

    // Build a global text by concatenating all references.
    // For each reference, record its global boundary.
    // We use '#' as a delimiter (assuming it does not appear in any sequence).
    string globalText = "";
    vector<ReferenceBoundary> boundaries;
    for (auto& r : refs) {
        ReferenceBoundary rb;
        rb.name = r.name;
        rb.start = globalText.size();
        globalText += r.sequence;
        rb.end = globalText.size() - 1; // last index of the sequence

        globalText.push_back('#');// Append the delimiter.
        boundaries.push_back(rb);
    }
    // globalText now contains all reference sequences separated by '#' characters

    // Build the Aho–Corasick automaton using the patterns.
    AhoCorasick ac;
    for (int i = 0; i < pats.size(); i++) {
        ac.patterns.push_back(pats[i].sequence);
        ac.addPattern(pats[i].sequence, i);
    }
    ac.buildAutomaton();

    // Search the global text.
    vector<vector<int>> globalMatches = ac.search(globalText);

    // Map matches back to their reference (and convert to local positions, 1-indexed).
    // For each pattern, we record a mapping: reference name -> vector of local positions.
    vector<unordered_map<string, vector<int>>> patternMatches(pats.size());
    for (int patIndex = 0; patIndex < globalMatches.size(); patIndex++) {
        int patLen = ac.patterns[patIndex].size();
        for (int pos : globalMatches[patIndex]) {
            pair<string, int> refInfo = getReferenceForPosition(pos, patLen, boundaries);
            if (refInfo.first != "") {
                patternMatches[patIndex][refInfo.first].push_back(refInfo.second);
            }
        }
    }

    // Write the output file (<outPrefix>.txt).
    // Each line is in the format:
    // (patternHeader) - ref1:pos1,pos2, ... , ref2:pos1, ...
    ofstream outFile(outPrefix + ".txt");
    for (int i = 0; i < pats.size(); i++) {
        outFile << "(" << pats[i].name << ") - ";
        bool firstRef = true;
        for (auto& entry : patternMatches[i]) {
            if (!firstRef) {
                outFile << ", ";
            }
            outFile << entry.first << ":";
            for (int j = 0; j < entry.second.size(); j++) {
                for (int i = 0; i < pats.size(); i++) {
                    if (j < entry.second.size() - 1) {
                        outFile << ",";
                    }
                }
                firstRef = false;
            }
            outFile << endl;
        }
        outFile.close();

        // If -d is specified, build the suffix tree for the global text and output its DOT file.
        if (dotFlag) {
            SuffixTree st(globalText);
            st.outputDot(outPrefix + ".dot");
        }

        return 0;
    }
