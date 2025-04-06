#include <iostream>
#include <vector>

using namespace std;

int main() {
    int N;
    cin >> N;
    vector<vector<int>> tree (N);
    for (int i = 0; i < N; i++) {
        int v, u;
        cin >> v >> u;
        tree[v].push_back(u);
        tree[u].push_back(v);
    }
    
    return 0;
}
