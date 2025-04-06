#include <iostream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;

int main() {
    int n;
    cin >> n;
    vector<int> result;
    int a;
    for (int i = 0; i < n; i++) {
        cin >> a;
        result.push_back(a);
    }
    if (result[0] == -1) {
        result[0] = 1;
    }
    for (int i = 1; i < n; i++) {
        if (result[i] == -1) {
            result[i] = result[i - 1] + 1;
        } else if (result[i] <= result[i - 1]){
            cout << "NO" << endl;
            return 0;
        }
    }
    cout << "YES" << endl;
    cout << result[0] << " ";
    for (int i = 1; i < n; i++) {
        cout << result[i] - result[i - 1] << " ";
    }
    return 0;
}