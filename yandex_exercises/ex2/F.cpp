#include <iostream>
#include <vector>

using namespace std;

int main() {
    int n;
    cin >> n;
    int value;
    vector<int> arr;
    vector<long long > prefix_sum(n + 1);
    for (int i = 0; i < n; i++) {
        cin >> value;
        arr.push_back(value);
        if (i != 0) {
            prefix_sum[i] = prefix_sum[i - 1] + arr[i - 1];
        }
    }
    prefix_sum[n] = prefix_sum[n - 1] + arr[n - 1];
    vector<long long> arr_pairs;
    for (int left = 1; left < n; left++) {
        arr_pairs.push_back(arr[left] * (prefix_sum[n] - prefix_sum[left + 1]));
    }
    vector<long long> prefix_pairs_sum(arr_pairs.size() + 1);
    for (int i = 1; i < arr_pairs.size() + 1; i++) {
        prefix_pairs_sum[i] = (prefix_pairs_sum[i - 1] + arr_pairs[i - 1]);
    }

    long long sum = 0;
    for (int i = 0; i < n; i++) {
        sum += (arr[i] * (prefix_pairs_sum[n - 1] - prefix_pairs_sum[i]));
    }
    cout << sum % 1000000007;
    return 0;
}

