// Görkem Kadir Solun 22003214

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <string>

using namespace std;

struct Node {
    // Node structure to store the row, column, and value of a matrix element
    // Also, it has pointers to the next and previous nodes in the circular linked list
    int row, col, value;
    Node* next, * prev;
    Node(int row, int col, int value) : row(row), col(col), value(value), next(nullptr), prev(nullptr) {}
};

class SparseMatrix {
private:
    Node* head;
    int size; // dimension of the matrix, n x n

    void clear() {
        if (!head) {
            return;
        }

        Node* current = head;
        do {
            Node* temp = current;
            current = current->next;
            delete temp;
        } while (current != head);
        head = nullptr;
    }

public:
    SparseMatrix(int n) : head(nullptr), size(n) {}
    ~SparseMatrix() {
        clear();
    }

    void insert(int row, int col, int value) {
        // If the value is 0, do not insert
        if (value == 0) {
            return;
        }

        Node* newNode = new Node(row, col, value);

        if (!head) {
            head = newNode;
            head->next = head;
            head->prev = head;
        } else {
            // Find the correct position to insert the new node to make the list sorted by row and column indices lower to higher
            // insert in reverse order to keep it efficient
            Node* a = head->prev;
            while (a != head && (a->row > row || (a->row == row && a->col > col))) {
                a = a->prev;
            }
            newNode->next = a->next;
            newNode->prev = a;
            a->next->prev = newNode;
            a->next = newNode;
        }
    }

    static SparseMatrix* createRandomMatrix(int n) {
        // Create an n x n sparse matrix with random values
        SparseMatrix* matrix = new SparseMatrix(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                matrix->insert(i, j, rand() % 1000);
            }
        }
        return matrix;
    }

    SparseMatrix* add(SparseMatrix* B) {
        // Sum of two sparse matrices A and B
        SparseMatrix* sumC = new SparseMatrix(size);

        if (!B->head || !head) {
            return sumC;
        }

        Node* a = head, * b = B->head;
        // Both lists are non-empty and sorted

        // Traverse the circular linked lists of A and B to calculate the sum
        // A and B are sorted by row and column indices
        do {
            if (a->row < b->row || (a->row == b->row && a->col < b->col)) {
                sumC->insert(a->row, a->col, a->value);
                a = a->next;
            } else if (a->row > b->row || (a->row == b->row && a->col > b->col)) {
                sumC->insert(b->row, b->col, b->value);
                b = b->next;
            } else {
                // same (row, col)
                sumC->insert(a->row, a->col, a->value + b->value);
                a = a->next;
                b = b->next;
            }
        } while (a != head && b != B->head);

        return sumC;
    }

    SparseMatrix* scalarMultiply(int scalar) {
        SparseMatrix* C = new SparseMatrix(size);
        if (scalar == 0 || !head) {
            return C;
        }

        Node* a = head;
        do {
            C->insert(a->row, a->col, a->value * scalar);
            a = a->next;
        } while (a != head);

        return C;
    }

    SparseMatrix* transpose() {
        SparseMatrix* T = new SparseMatrix(size);

        if (!head) {
            return T;
        }

        Node* a = head;
        do {
            // Transpose the matrix by swapping the row and column indices
            T->insert(a->col, a->row, a->value);
            a = a->next;
        } while (a != head);

        return T;
    }

    void print(const string& name) {
        cout << name << ":\n";
        if (!head) {
            // Matrix is empty, so everything is 0
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    cout << "0 ";
                }
                cout << "\n";
            }
            return;
        }

        Node* a = head;
        // Traverse the circular linked list to print the matrix
        // If a node with the same row and column indices is found, print its value
        // Otherwise, print 0
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (a->row == i && a->col == j) {
                    cout << a->value << " ";
                    a = a->next;
                } else {
                    cout << "0 ";
                }
            }
            cout << "\n";
        }
    }
};

int main(int argc, char* argv[]) {
    int n = 1, s = INT32_MAX;
    for (int i = 1; i < argc; i++) {
        if (string(argv[i]) == "-n" && i + 1 < argc) {
            n = atoi(argv[i + 1]);
        } else if (string(argv[i]) == "-s" && i + 1 < argc) {
            s = atoi(argv[i + 1]);
        }
    }
    if (n < 1 || s == INT32_MAX) {
        cout << "Usage: make\n./basic_matrix_operations -n <size> -s <scalar>\nmake clean\n";
        return 1;
    }

    cout << "# Addition of two " << n << "x" << n << " random matrices\n";
    SparseMatrix* A = SparseMatrix::createRandomMatrix(n);
    SparseMatrix* B = SparseMatrix::createRandomMatrix(n);
    SparseMatrix* C = A->add(B);
    A->print("A");
    B->print("B");
    C->print("C");

    cout << "\n# Scalar multiplication of " << n << "x" << n << " matrix by " << s << "\n" << "s: " << s << "\n";
    SparseMatrix* D = SparseMatrix::createRandomMatrix(n);
    SparseMatrix* E = D->scalarMultiply(s);
    D->print("A");
    E->print("C");

    cout << "\n# Transposition of a random matrix\n";
    SparseMatrix* F = SparseMatrix::createRandomMatrix(n);
    SparseMatrix* G = F->transpose();
    F->print("A");
    G->print("A.T");

    delete A;
    delete B;
    delete C;
    delete D;
    delete E;
    delete F;
    delete G;

    return 0;
}
