#include <iostream>
#include <map>
#include <vector>
#include <unordered_map>
#include <fstream>

using namespace std;

map<int, long long> result;
unordered_map<long long, vector<long long>> tree;

pair<long, long> solve(int root) {
    long long money = 1;
    long long count_children = tree[root].size();
    for (long child: tree[root]) {
        pair<long, long> money_and_count_children = solve(child);
        money += money_and_count_children.first;
        count_children += money_and_count_children.second;
    }
    result[root] = money + count_children;
    return make_pair(money + count_children, count_children);
}

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    fstream in("input.txt");
    int N;
    in >> N;
    for (int i = 2; i < N + 1; i++) {
        int value;
        in >> value;
        tree[value].push_back(i);
    }
    solve(1);
    fstream out("output.txt");
    for (pair<long, long> value: result) {
        out << value.second << " ";
    }
    return 0;
}
