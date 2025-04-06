#include <iostream>

using namespace std;

int main() {
    long long l, r;
    cin >> l >> r;
    int count = 0;
    for (long long i = l; i <= r; i++) {
        int divisors = 0;
        for (long long j = 1; j * j <= i; j++) {
            if (i % j == 0) {
                divisors++;
                if (j != i / j) divisors++;
            }
        }
        if (divisors > 2) {
            bool isPrime = true;
            for (int k = 2; k * k <= divisors; k++) {
                if (divisors % k == 0) {
                    isPrime = false;
                    break;
                }
            }
            if (isPrime) {
                count++;
            }
        }
    }
    cout << count << std::endl;
    return 0;
}
