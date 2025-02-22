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
    // output holds indices of patterns that end at this node.
    vector<int> output;
    ACNode() : failure(nullptr) {}
};

class AhoCorasick {
public:
    ACNode* root;
    vector<string> patterns; // store patterns by index

    AhoCorasick() {
        root = new ACNode();
    }

    // Insert pattern into trie and associate it with its index.
    void addPattern(const string& pat, int index) {
        ACNode* node = root;
        for (char c : pat) {
            if (node->children.find(c) == node->children.end())
                node->children[c] = new ACNode();
            node = node->children[c];
        }
        node->output.push_back(index);
    }

    // Build failure links using BFS.
    void buildAutomaton() {
        queue<ACNode*> q;
        // For all children of root, failure is root.
        for (auto& p : root->children) {
            p.second->failure = root;
            q.push(p.second);
        }
        while (!q.empty()) {
            ACNode* curr = q.front();
            q.pop();
            for (auto& p : curr->children) {
                char c = p.first;
                ACNode* child = p.second;
                ACNode* f = curr->failure;
                while (f && f->children.find(c) == f->children.end())
                    f = f->failure;
                if (f)
                    child->failure = f->children[c];
                else
                    child->failure = root;
                // Merge outputs so that if failure node is terminal, we also output that pattern.
                child->output.insert(child->output.end(),
                    child->failure->output.begin(),
                    child->failure->output.end());
                q.push(child);
            }
        }
    }

    // Search the text; returns a vector (for each pattern) of 1-indexed positions where the pattern occurs.
    vector<vector<int>> search(const string& text) {
        vector<vector<int>> matches(patterns.size());
        ACNode* node = root;
        for (int i = 0; i < text.size(); i++) {
            char c = text[i];
            // Follow failure links if no matching child.
            while (node && node->children.find(c) == node->children.end())
                node = node->failure;
            if (!node) {
                node = root;
                continue;
            }
            node = node->children[c];
            for (int patIndex : node->output) {
                // Record match; positions are output 1-indexed.
                int pos = i - patterns[patIndex].size() + 1;
                matches[patIndex].push_back(pos + 1);
            }
        }
        return matches;
    }

    // Free memory recursively.
    void freeNode(ACNode* node) {
        for (auto& p : node->children)
            freeNode(p.second);
        delete node;
    }

    ~AhoCorasick() {
        freeNode(root);
    }
};

struct Reference {
    string name;
    string sequence;
};

vector<Reference> readFastaReference(const string& filename) {
    vector<Reference> refs;
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening reference file: " << filename << endl;
        exit(1);
    }

    string line;
    Reference current;
    bool inSequence = false;
    while (getline(infile, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            if (inSequence) {
                refs.push_back(current);
                current = Reference();
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

struct Pattern {
    string name;
    string pattern;
};

vector<Pattern> readFastaPatterns(const string& filename) {
    vector<Pattern> pats;
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening pattern file: " << filename << endl;
        exit(1);
    }
    string line;
    Pattern current;
    bool expectPattern = false;
    while (getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            current.name = line.substr(1);
            expectPattern = true;
        } else {
            if (expectPattern) {
                current.pattern = line;
                pats.push_back(current);
                expectPattern = false;
            }
        }
    }
    infile.close();
    return pats;
}

// Each node in the suffix tree.
struct STNode {
    // Children map: key is the first character of the edge.
    unordered_map<char, STNode*> children;
    // Edge label is represented by [start, end] indices into the text.
    int start, * end;
    // For leaves, store the starting index of the suffix.
    vector<int> suffixIndices;
    // Node ID for DOT output.
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

    // Naively insert all suffixes into the tree.
    void build() {
        int n = text.size();
        for (int i = 0; i < n; i++) {
            insertSuffix(i);
        }
        // Compress the tree (merge chains of nodes with single child).
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
            // Merge if child is not a leaf and has exactly one child.
            while (child->children.size() == 1) {
                auto it = child->children.begin();
                STNode* grandchild = it->second;
                // Extend child's edge to include grandchild’s edge.
                child->end = grandchild->end;
                // Merge the children.
                child->children = grandchild->children;
                // Append any suffix indices.
                child->suffixIndices.insert(child->suffixIndices.end(),
                    grandchild->suffixIndices.begin(),
                    grandchild->suffixIndices.end());
                // Continue merging if possible.
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
        // Print node label.
        if (node == root)
            out << "node_" << node->id << " [label=\"\"];" << endl;
        else {
            string label = text.substr(node->start, *(node->end) - node->start + 1);
            out << "node_" << node->id << " [label=\"" << label << "\"];" << endl;
        }
        // Print edges.
        for (auto& childPair : node->children) {
            STNode* child = childPair.second;
            string edgeLabel = text.substr(child->start, *(child->end) - child->start + 1);
            out << "node_" << node->id << " -> node_" << child->id
                << " [label=\"" << edgeLabel << "\"];" << endl;
            outputNodeDot(child, out);
        }
    }

    // Free the nodes recursively.
    void freeNode(STNode* node) {
        for (auto& child : node->children)
            freeNode(child.second);
        if (node->end) delete node->end;
        delete node;
    }
};

int main(int argc, char* argv[]) {
    string refFile, patFile, outPrefix;
    bool dotFlag = false;

    // Simple command-line argument parsing.
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "-r" && i + 1 < argc)
            refFile = argv[++i];
        else if (arg == "-p" && i + 1 < argc)
            patFile = argv[++i];
        else if (arg == "-o" && i + 1 < argc)
            outPrefix = argv[++i];
        else if (arg == "-d")
            dotFlag = true;
    }

    if (refFile.empty() || patFile.empty() || outPrefix.empty()) {
        cerr << "Usage: " << argv[0] << " -r <reference file> -p <pattern file> -o <output prefix> [-d]" << endl;
        return 1;
    }

    // Read reference and pattern files.
    vector<Reference> references = readFastaReference(refFile);
    vector<Pattern> patterns = readFastaPatterns(patFile);

    // Build Aho–Corasick automaton.
    AhoCorasick ac;
    for (int i = 0; i < patterns.size(); i++) {
        ac.patterns.push_back(patterns[i].pattern);
        ac.addPattern(patterns[i].pattern, i);
    }
    ac.buildAutomaton();

    // For each reference, run the automaton and store matches.
    // We'll record for each pattern a mapping from reference name to a list of match positions.
    vector<unordered_map<string, vector<int>>> patternMatches(patterns.size());

    for (auto& ref : references) {
        string text = ref.sequence;
        // For Aho–Corasick, we use the raw sequence.
        vector<vector<int>> matches = ac.search(text);
        for (int i = 0; i < matches.size(); i++) {
            if (!matches[i].empty())
                patternMatches[i][ref.name] = matches[i];
        }
    }

    // Write output file <prefix>.txt
    ofstream outFile(outPrefix + ".txt");
    for (int i = 0; i < patterns.size(); i++) {
        outFile << "(" << patterns[i].name << ") - ";
        bool firstRef = true;
        for (auto& entry : patternMatches[i]) {
            if (!firstRef)
                outFile << ", ";
            outFile << entry.first << ":";
            for (int j = 0; j < entry.second.size(); j++) {
                outFile << entry.second[j];
                if (j < entry.second.size() - 1)
                    outFile << ",";
            }
            firstRef = false;
        }
        outFile << endl;
    }
    outFile.close();

    // If the -d flag is provided, output a DOT file for each reference’s suffix tree.
    if (dotFlag) {
        for (auto& ref : references) {
            string text = ref.sequence;
            // Append terminator '$' if not present.
            if (text.empty() || text.back() != '$')
                text.push_back('$');
            SuffixTree st(text);
            // Create a DOT file; you might wish to sanitize ref.name if it has spaces.
            string dotFilename = outPrefix + "_" + ref.name + ".dot";
            st.outputDot(dotFilename);
        }
    }

    return 0;
}
