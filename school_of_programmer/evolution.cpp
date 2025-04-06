#include <iostream>

using namespace std;

int main() {
    int N, x, y;
    cin >> N >> x >> y;
    while (x != y) {
        if (x > y) {
            x /= 2;
        } else {
            y /= 2;
        }
    }
    cout << x;
    return 0;
}
