#include <iostream>
#include <vector>

using namespace std;

int main() {
    int n, k;
    cin >> n >> k;
    int x;
    vector<int> prefix_sum(n + 1);
    for (int i = 1; i < n + 1; i++) {
        cin >> x;
        prefix_sum[i] = prefix_sum[i - 1] + x;
    }
    int count = 0;
    int right = 0;
    for (int left = 0; left < n + 1; left++) {
        while (right != n && prefix_sum[right] - prefix_sum[left] < k) {
            right++;
        }
        if (prefix_sum[right] - prefix_sum[left] == k) {
            count++;
        }
    }
    cout << count;
    return 0;
}