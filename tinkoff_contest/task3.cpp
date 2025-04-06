#include <iostream>
#include <map>
#include <unordered_map>
#include <string>
#include <algorithm>

using namespace std;

int main() {
    string str;
    cin >> str;
    string req;
    cin >> req;
    int k;
    cin >> k;
    unordered_map<char, int> right_ch;
    vector<string> strings;
    string l = "";
    for (int i = 0; i < str.size(); i++) {
        if (right_ch.count(str[i]) == 0) {
            if (!l.empty()) {
                strings.push_back(l);
            }
        }
        else {
            l += str[i];
        }
    }
    for (int i = (int)strings.size() - 1; i >= 0 ; i--) {
        str = strings[i];
        unordered_map<char, int> h;
        for (int j = str.size() - 1; j >= 0; j--) {
            h[str[j]]++;
            if (k + j <= str.size() - 1) {
                h[str[j]]--;
                if (h[str[j + k]] == 0) {
                    h.erase(str[j + k]);
                }
            }
            if ((int) h.size() == (int) right_ch.size()) {
                int min = (j + k < str.size()) ? j + k : str.size();
                for (int m = j; m < min; k--) {
                    cout << str[j];
                }
                return 0;
            }
        }
    }
    cout << -1;
    return 0;
}
