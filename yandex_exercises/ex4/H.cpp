#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

vector<int> costs;
vector<int> result;
vector<vector<int>> tree;
vector<vector<long long>> dp;

void coating(int v, int parent) {
    if (tree[v].size() == 1) {
        dp[v][0] = -1;
        dp[v][1] = costs[v];
    } else {
        dp[v][0] = 0;
        dp[v][1] = costs[v];
    }
    for (int u: tree[v]) {
        if (u == parent) continue;
        coating(u, v);
        if (dp[v][0] == -1) {
            dp[v][0] = 0;
            dp[v][0] += dp[u][1];
        } else {
            dp[v][0] += dp[u][1];
        }
        if (dp[u][0] != -1) {
            dp[v][1] += min(dp[u][0], dp[u][1]);
        }
    }
}

void dfs(int v, int parent, bool is_covered_parent) {
    if (parent == -1) {
        if (dp[v][0] >= dp[v][1]) {
            result.push_back(v + 1);
            is_covered_parent = true;
        }
    } else {
        if (is_covered_parent) {
            if (dp[v][0] >= dp[v][1]) {
                result.push_back(v + 1);
                is_covered_parent = true;
            } else {
                is_covered_parent = false;
            }
        } else {
            is_covered_parent = true;
            result.push_back(v + 1);
        }
    }
    for (int u: tree[v]) {
        if (u == parent) continue;
        dfs(u, v, is_covered_parent);
    }
}

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int N;
    cin >> N;
    costs.resize(N);
    tree.resize(N);
    dp.assign(N, vector<long long>(2, 0));
    for (int i = 0; i < N - 1; ++i) {
        int u, v;
        cin >> u >> v;
        --u, --v;
        tree[u].push_back(v);
        tree[v].push_back(u);
    }
    for (int i = 0; i < N; ++i) {
        cin >> costs[i];
    }
    if (N == 1) {
        cout << costs[0] << " " << 1 << endl;
        cout << 1;
        return 0;
    }
    coating(0, -1);
    dfs(0, -1, false);
    cout << min(dp[0][0], dp[0][1]) << " " << result.size() << "\n";
    for (int v: result) {
        cout << v << " ";
    }
    return 0;
}
