#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <deque>
#include <climits>

using namespace std;

struct Chair {
    int height;
    int width;
} typedef Chair;

bool comparator(Chair x, Chair y) {
    return x.height < y.height;
}

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int n, H;
    cin >> n >> H;
    vector<int> heights;
    vector<int> widths;
    vector<Chair> chairs;
    for (int i = 0; i < n; i++) {
        int value;
        cin >> value;
        heights.push_back(value);
    }
    for (int i = 0; i < n; i++) {
        int value;
        cin >> value;
        widths.push_back(value);
    }
    for (int i = 0; i < n; i++) {
        Chair c = {.height = heights[i], .width = widths[i]};
        chairs.push_back(c);
    }
    sort(chairs.begin(), chairs.end(), comparator);
    int right = 0;
    int left = 0;
    int w = chairs[0].width;
    deque<int> d;
    int min_convenience = INT_MAX;
    for (int i = 0; i < n; i++) {
        if (chairs[i].width >= H) {
            min_convenience = 0;
            cout << min_convenience;
            return 0;
        }
    }
    while (right < n) {
        if (w < H || left >= right) {
            right++;
            w += chairs[right].width;
            if (right > 0) {
                int dif = fabs(chairs[right].height - chairs[right - 1].height);
                while (!d.empty() && d.back() < dif) {
                    d.pop_back();
                }
                d.push_back(dif);
            }
        } else {
            min_convenience = min(d.front(), min_convenience);
            if (left < n - 1) {
                int dif = fabs(chairs[left + 1].height - chairs[left].height);
                if (dif == d.front()) {
                    d.pop_front();
                }
            }
            w -= chairs[left].width;
            left++;
        }
    }
    cout << min_convenience;
    return 0;
}
