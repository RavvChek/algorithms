#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>
#include <fstream>

using namespace std;

int main() {
    int n, k;
    cin >> n >> k;
    int value;
    vector<int> arr;
    for (int i = 0; i < n; i++) {
        cin >> value;
        arr.push_back(value);
    }
    sort(arr.begin(), arr.end());
    int right = 0;
    int left = 0;
    int day_count = 1;
    int count = 0;
    while (right < n) {
        if (arr[right] - arr[left] <= k) {
            right++;
            count++;
        } else {
            day_count = max(day_count, count);
            left++;
            count--;
        }
    }
    day_count = max(day_count, count);
    cout << day_count;
    return 0;
}