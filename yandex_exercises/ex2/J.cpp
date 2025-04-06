#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int n;
    cin >> n;
    int value;
    vector<int> evidence(n);
    for (int i = 0; i < n; i++) {
        cin >> value;
        evidence[i] = value;
    }
    int m, k;
    cin >> m >> k;
    vector<int> experiments(m);
    for (int i = 0; i < m; i++) {
        cin >> value;
        experiments[i] = value;
    }
    vector<int> prefix_sum_count_rep(n);
    unordered_map<int, int> map_count_repeat;
    for (int i = 1; i < n; i++) {
        if (evidence[i] == evidence[i - 1]) {
            prefix_sum_count_rep[i] = prefix_sum_count_rep[i - 1] + 1;
            map_count_repeat.insert({prefix_sum_count_rep[i], i});
        } else {
            prefix_sum_count_rep[i] = prefix_sum_count_rep[i - 1];
        }
    }
    vector<int> arr_facts(n);
    unordered_map<int, int> map_facts;
    arr_facts[0] = 0;
    for (int i = 1; i < n; i++) {
        if (evidence[i] < evidence[i - 1]) {
            arr_facts[i] = arr_facts[i - 1] + 1;
            map_facts.insert({arr_facts[i], i});
        } else {
            arr_facts[i] = arr_facts[i - 1];
        }
    }
    int right;
    for (int i = 0; i < m; i++) {
        int left = 0;
        right = experiments[i] - 1;
        int border_repeat = map_count_repeat[prefix_sum_count_rep[right] - k];
        int border_increase;
        if (arr_facts[right] - 1 < 0) {
            border_increase = 0;
        } else {
            border_increase = map_facts[arr_facts[right]];
        }
        if (border_repeat > border_increase) {
            cout << border_repeat + 1 << " ";
        } else {
            cout << border_increase + 1 << " ";
        }
    }
    return 0;
}
