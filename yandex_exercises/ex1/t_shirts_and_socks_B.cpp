#include <iostream>

using namespace std;

int main() {
    int A, B, C, D;
    cin >> A >> B >> C >> D;
    int max_AB = max(A, B);
    int max_CD = max(C, D);
    if (A + C < B + D && A + C < min(max_AB, max_CD) && (A != 0 && B != 0 && C != 0 && D != 0)) {
        cout << A + 1 << " " << C + 1;
    } else if (B + D < A + C && B + D < min(max_AB, max_CD) && (A != 0 && B != 0 && C != 0 && D != 0)) {
        cout << B + 1 << " " << D + 1;
    } else if (max_AB < max_CD && (A != 0 && B != 0 && C != 0 && D != 0)) {
        cout << max_AB + 1 << " " << 1;
    } else if (max_AB > max_CD && (A != 0 && B != 0 && C != 0 && D != 0)){
        cout << 1 << " " << max_CD + 1;
    } else if (A == 0) {
        cout << 1 << " " << C + 1;
    } else if (B == 0) {
        cout << 1 << " " << D + 1;
    } else if (C == 0) {
        cout << A + 1 << " " << 1;
    } else {
        cout << B + 1 << " " << 1;
    }
    return 0;
}