#include <iostream>
#include <vector>

using namespace std;

int main() {
    int n;
    cin >> n;
    int x;
    vector<int> arr;
    for (int i = 0; i < n; i++) {
        cin >> x;
        arr.push_back(x);
    }
    vector<int> prefix_sum(n + 1);
    prefix_sum[0] = 0;
    for (int i = 1; i < n + 1; i++) {
        prefix_sum[i] = prefix_sum[i - 1] + arr[i - 1];
    }
    for (int i = 1; i < n + 1; i++) {
        cout << prefix_sum[i] << " ";
    }
    return 0;
}