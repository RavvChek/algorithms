#include <iostream>
#include <vector>
#include <climits>

using namespace std;

int main() {
    int n;
    cin >> n;
    int value;
    vector<int> arr(n);
    vector<long long> prefix_sum(n + 1);
    for (int i = 0; i < n; i++) {
        cin >> value;
        arr[i] = value;
    }
    for (int i = 1; i < n + 1; i++) {
        prefix_sum[i] = prefix_sum[i - 1] + arr[i - 1];
    }
    long long count_transitions = 0;
    long long min_count_transitions = INT_MAX;
    for (int i = 1; i < n; i++) {
        count_transitions += (long long)i * arr[i];
    }
    min_count_transitions = count_transitions;
    for (int i = 1; i < n; i++) {
        count_transitions =
                count_transitions - arr[i] + (prefix_sum[i] - prefix_sum[0]) - (prefix_sum[n] - prefix_sum[i + 1]);
        min_count_transitions = min(count_transitions, min_count_transitions);
    }
    cout << min_count_transitions;
    return 0;
}