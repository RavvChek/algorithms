#include <iostream>
#include <vector>
#include <deque>

using namespace std;

int main() {
    int n, k;
    cin >> n >> k;
    vector<int> arr;
    vector<int> result;
    deque<int> d;
    for (int i = 0; i < n; i++) {
        int value;
        cin >> value;
        arr.push_back(value);
    }
    for (int i = 0; i < k; i++) {
        while (!d.empty() && d.back() > arr[i]) {
            d.pop_back();
        }
        d.push_back(arr[i]);
    }
    int min = d.front();
    result.push_back(min);
    for (int i = k; i < n; i++) {
        while (!d.empty() && d.back() > arr[i]) {
            d.pop_back();
        }
        if (!d.empty() && arr[i - k] == min) {
            d.pop_front();
        }
        d.push_back(arr[i]);
        min = d.front();
        result.push_back(min);
    }
    for (int res: result) {
        cout << res << endl;
    }
    return 0;
}
