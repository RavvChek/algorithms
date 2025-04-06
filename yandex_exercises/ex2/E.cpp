#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    int n;
    cin >> n;
    int value;
    vector<int> arr;
    for (int i = 0; i < n; i++) {
        cin >> value;
        arr.push_back(value);
    }
    sort(arr.begin(), arr.end());
    int l = 0;
    int r = arr.size() - 1;
    int middle = (l + r) / 2;
    while (!arr.empty()) {
        if (arr.size() / 2 != 0) {
            cout << arr[middle] << " ";
            arr.erase(arr.begin() + middle);
        } else {
            cout << arr[middle] << " ";
            arr.erase(arr.begin() + middle);
        }
        r--;
        middle = (l + r) / 2;
    }
    return 0;
}