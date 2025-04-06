#include <iostream>
#include <vector>

using namespace std;

int main() {
    int n, r;
    cin >> n >> r;
    vector<int> prefix_sum(n + 1);
    int x;
    for (int i = 1; i < n + 1; i++) {
        cin >> x;
        prefix_sum[i] = x;
    }
    long count = 0;
    int right = 2;
    for (int left = 1; left < n + 1; left++) {
        while (right < n + 1 && prefix_sum[right] - prefix_sum[left] <= r) {
            right++;
        }
        count += n + 1 - right;
    }
    cout << count;
    return 0;
}