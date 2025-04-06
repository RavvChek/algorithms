#include <iostream>

using namespace std;

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int n, b;
    cin >> n >> b;
    long long time = 0;
    long long count_visitors = 0;
    for (int i = 0; i < n; i++) {
        int value;
        cin >> value;
        count_visitors += value;
        time += count_visitors;
        if (count_visitors >= b) {
            count_visitors -= b;
        } else {
            count_visitors = 0;
        }
    }
    time += count_visitors;
    cout << time;
    return 0;
}
