#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

int main() {
    int n;
    cin >> n;
    vector<long long> execution_time(n + 1, 0);
    vector<vector<int>> dependencies(n + 1);
    vector<int> indegree(n + 1, 0);
    vector<long long> completion_time(n + 1, 0);
    for (int i = 1; i <= n; ++i) {
        long long ti;
        cin >> ti;
        execution_time[i] = ti;
        int dep;
        while (cin >> dep) {
            if (dep == -1) break;
            dependencies[dep].push_back(i);
            ++indegree[i];
        }
    }
    queue<int> q;
    for (int i = 1; i <= n; ++i) {
        if (indegree[i] == 0) {
            completion_time[i] = execution_time[i];
            q.push(i);
        }
    }
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int v : dependencies[u]) {
            completion_time[v] = max(completion_time[v], completion_time[u] + execution_time[v]);
            --indegree[v];
            if (indegree[v] == 0) {
                q.push(v);
            }
        }
    }
    long long result = *max_element(completion_time.begin(), completion_time.end());
    cout << result << endl;
    return 0;
}
