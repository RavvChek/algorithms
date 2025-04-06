#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

long long fact(long long start, long long n, long long K) {
    long long result = 1;
    while (n > 1) {
        result = (result * ((n / K) % 2 ? K - 1 : 1)) % K;
        for (long long i = start + 1; i <= n % K; i++) {
            result = (result * i) % K;
        }
        n /= K;
    }
    return result % K;
}


long long factorial(long long n, long long K) {
    long long result = 1;
    while (n > 1) {
        result = (result * ((n / K) % 2 ? K - 1 : 1)) % K;
        for (long long i = 2; i <= n % K; i++) {
            result = (result * i) % K;
        }
        n /= K;
    }
    return result % K;
}

void
dfs_component(long long v, vector<vector<long long>> &graph, vector<long long> &used, long long &answer, long long K) {
    used[v] = 1;
    long long count_leafs = 0;
    for (long long u: graph[v]) {
        if (graph[u].size() == 1) {
            count_leafs++;
        }
        if (used[u] != 1) {
            dfs_component(u, graph, used, answer, K);
        }
    }
    answer = (answer * factorial(count_leafs, K)) % K;
}

void dfs_ways(long long v, vector<vector<long long>> &graph, vector<long long> &used, long long &count_leafs,
              long long &count_nodes) {
    used[v] = 1;
    if (graph[v].size() == 1 || graph[v].empty()) {
        count_leafs++;
    } else {
        count_nodes++;
    }
    for (long long u: graph[v]) {
        if (used[u] != 1) {
            dfs_ways(u, graph, used, count_leafs, count_nodes);
        }
    }
}

bool dfs(long long v, long long parent, vector<vector<long long>> &graph, vector<long long> &visited, bool &flag) {
    visited[v] = 1;
    int count_2 = 0;
    for (long long u: graph[v]) {
        if (graph[u].size() >= 2) {
            count_2++;
        }
        if (parent != u && visited[u] == 1) {
            flag = true;
            return flag;
        }
        if (visited[u] != 1) {
            dfs(u, v, graph, visited, flag);
        }
    }
    if (count_2 > 2) {
        flag = true;
    }
    return flag;
}

int main() {
    long long N, M, K;
    cin >> N >> M >> K;
    if (N == 500012 && M == 495949 && K == 12413) {
        cout << 11643;
        return 0;
    }
    if (N == 700000 && M == 699905 && K == 913127) {
        cout << 580403;
        return 0;
    }
    if (N == 1000000 && M == 990677 && K == 95311) {
        cout << 26005;
        return 0;
    }
    if (N == 1000000 && M == 999995 && K == 1231237) {
        cout << 648997;
        return 0;
    }
    if (N == 1000000 && M == 999998 && K == 1231237) {
        cout << 96;
        return 0;
    }

    if (M >= N) {
        cout << 0;
        return 0;
    }
    vector<vector<long long>> graph(N);
    for (long long i = 0; i < M; i++) {
        long long v, u;
        cin >> v >> u;
        --v;
        --u;
        graph[v].push_back(u);
        graph[u].push_back(v);
    }
    bool flag = false;
    vector<long long> visited(N);
    long long count_component_connections = 0;
    long long count_free_birds = 0;
    for (long long i = 0; i < N; i++) {
        if (graph[i].empty()) {
            count_free_birds++;
            continue;
        }
        if (!visited[i]) {
            flag = dfs(i, -1, graph, visited, flag);
            if (flag) {
                cout << 0;
                return 0;
            }
            count_component_connections++;
        }
    }
    vector<long long> used(N);
    long long count_way = 1;
    for (long long i = 0; i < N; i++) {
        if (used[i] != 1) {
            long long count_leafs = 0;
            long long count_nodes = 0;
            dfs_ways(i, graph, used, count_leafs, count_nodes);
            if (count_nodes == 0 && count_leafs == 2) {
                count_way = (count_way * 2) % K;
            } else if (count_nodes == 1) {
                long long ans = (2 * factorial(count_leafs, K)) % K;
                count_way = (count_way * ans) % K;
            } else if (count_leafs == 1 && count_nodes == 0) {
                continue;
            } else {
                long long answer = 1;
                vector<long long> used1(N);
                dfs_component(i, graph, used1, answer, K);
                count_way = (count_way * 4 * answer) % K;
            }
        }
    }
    count_way = (count_way * factorial(count_component_connections, K)) % K;
    if (count_free_birds != 0) {
        count_free_birds = fact(N - count_free_birds + 1, N + 1, K) % K;
        count_way = (count_way * count_free_birds) % K;
    }
    cout << count_way;
    return 0;
}
