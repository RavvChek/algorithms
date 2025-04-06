#include <iostream>
#include <vector>

using namespace std;

int calculate_size(const vector<vector<int>> &tree, int current, int parent, vector<int> &result) {
    int count = 1;
    if (parent < tree.size() || parent == -1) {
        for (int v: tree[current]) {
            if (parent != v) {
                count += calculate_size(tree, v, current, result);
            }
        }
    }
    result[current] = count;
    return count;
}

int main() {
    int V;
    cin >> V;
    vector<int> result (V);
    vector<vector<int>> tree(V);
    for (int i = 0; i < V - 1; i++) {
        int v1, v2;
        cin >> v1 >> v2;
        tree[v1 - 1].push_back(v2 - 1);
        tree[v2 - 1].push_back(v1 - 1);
    }
    calculate_size(tree, 0, -1, result);
    for (int value: result) {
        cout << value << " ";
    }
    return 0;
}
