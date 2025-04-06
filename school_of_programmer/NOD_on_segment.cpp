#include <iostream>
#include <vector>

using namespace std;

struct Node {
    int value;
    int l;
    int r;
    Node *left;
    Node *right;
} typedef Node;

int NOD(int A, int B) {
    while (B != 0) {
        int c = B;
        B = A % B;
        A = c;
    }
    return A;
}

Node *build(const vector<int> &array, int l, int r) {
    if (l == r) {
        Node *leaf = new Node();
        leaf->value = array[l];
        leaf->r = r;
        leaf->l = l;
        return leaf;
    }
    Node *node = new Node();
    node->l = l;
    node->r = r;
    int middle = (r + l) / 2;
    node->left = build(array, l, middle);
    node->right = build(array, middle + 1, r);
    node->value = NOD(node->left->value, node->right->value);
    return node;
}

int find(const vector<int> &array, Node *node, int l, int r) {
    if (l > r) {
        return 0;
    }
    if (l == node->l && r == node->r) {
        return node->value;
    }
    int middle = (node->l + node->r) / 2;
    if (r <= middle) {
        return find(array, node->left, l, r);
    } else if (l > middle) {
        return find(array, node->right, l, r);
    } else {
        int left_gcd = find(array, node->left, l, middle);
        int right_gcd = find(array, node->right, middle + 1, r);
        return NOD(left_gcd, right_gcd);
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int N, M;
    cin >> N;
    vector<int> array(N);
    for (int i = 0; i < N; i++) {
        cin >> array[i];
    }
    cin >> M;
    Node *node = build(array, 0, N - 1);
    for (int i = 0; i < M; i++) {
        int L, R;
        cin >> L >> R;
        cout << find(array, node, L - 1, R - 1) << " ";
    }
    return 0;
}
