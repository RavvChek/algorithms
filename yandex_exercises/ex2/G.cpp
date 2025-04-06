#include <iostream>
#include <string>

using namespace std;

int main() {
    int n, c;
    cin >> n >> c;
    string s;
    cin >> s;
    long long count_rude = 0;
    int max_len = 0;
    int len = 0;
    int right = 0;
    int left = 0;
    int count_a = 0;
    int count_b = 0;
    while (right < n) {
        if (count_rude <= c) {
            if (s[right] == 'a') {
                count_a++;
            } else if (s[right] == 'b'){
                count_b++;
                count_rude += count_a;
            }
            max_len = max(max_len, len);
            len++;
            right++;
        } else {
            if (s[left] == 'a') {
                count_a--;
                count_rude -= count_b;
            } else if (s[left] == 'b'){
                count_b--;
            }
            len--;
            max_len = max(max_len, len);
            left++;
        }
    }
    if (count_rude <= c) {
        max_len = max(max_len, len);
    }
    cout << max_len;
    return 0;
}
