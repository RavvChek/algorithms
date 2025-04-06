#include <iostream>

using namespace std;

int main() {
    int x, y, x1, y1, x2, y2;
    cin >> x1 >> y1 >> x2 >> y2 >> x >> y;
    if (x > x1 && x < x2) {
        if (y > y2) {
            cout << "N";
        }
        else {
            cout << "S";
        }
    }
    else if (x > x2) {
        if (y < y2 && y > y1) {
            cout << "E";
        }
        else if (y > y2) {
            cout << "NE";
        } else {
            cout << "SE";
        }
    }
    else {
        if (y < y2 && y > y1) {
            cout << "W";
        }
        else if (y > y2) {
            cout << "NW";
        } else {
            cout << "SW";
        }
    }
    return 0;
}
