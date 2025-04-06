#include <iostream>

using namespace std;

struct TreeNode {
    int value;
    TreeNode *left;
    TreeNode *right;

    TreeNode(int value) : value(value), left(nullptr), right(nullptr) {};
} typedef TreeNode;

string search(TreeNode* root, int value) {
    if (!root) {
        return "NO";
    }
    if (root->value == value) {
        return "YES";
    }
    if (value < root->value) {
        return search(root->left, value);
    } else {
        return search(root->right, value);
    }
}


string add(TreeNode *& root, int value, string answer) {
    if (answer == "YES") {
        return "ALREADY";
    } else {
        if (!root) {
            root = new TreeNode (value);
            return "DONE";
        }
        if (value < root->value) {
            return add(root->left, value, answer);
        } else {
            return add(root->right, value, answer);
        }
    }
}

void print_tree(TreeNode * root, int count) {
    if (!root) {
        return;
    }
    print_tree(root->left, count + 1);
    string dots;
    for (int i = 0; i < count; i++) {
        dots += '.';
    }
    cout << dots << root->value << endl;
    print_tree(root->right, count + 1);

}


int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    string command;
    string operand;
    int n;
    TreeNode *tree = nullptr;
    while (cin >> command) {
        if (command == "ADD") {
            cin >> operand;
            n = stoi(operand);
            string answer = search(tree, n);
            cout << add(tree, n, answer) << endl;
        } else if (command == "SEARCH") {
            cin >> operand;
            n = stoi(operand);
            cout << search(tree, n) << endl;
        } else {
            int count = 1;
            print_tree(tree->left, count);
            cout << tree->value << endl;
            print_tree(tree->right, count);
        }
    }
    return 0;
}
