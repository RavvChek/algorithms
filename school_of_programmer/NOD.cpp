#include <iostream>

using namespace std;

int main() {
    int A, B;
    cin >> A >> B;
    while (B != 0) {
        int c = B;
        B = A % B;
        A = c;
    }
    cout << A;
    return 0;
}
