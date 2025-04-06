#include <iostream>
#include <stack>
#include <vector>

using namespace std;

struct cort {
    int index;
    int value;
} typedef cort;

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int N;
    cin >> N;
    stack<cort> s;
    vector<int> result(N);
    for (int i = 0; i < N; i++) {
        int value;
        cin >> value;
        if (!s.empty()) {
            while (!s.empty() && s.top().value > value) {
                int index = s.top().index;
                result[index] = i;
                s.pop();
            }
            cort c = {.index = i, .value = value};
            s.push(c);
        } else {
            cort c = {.index = i, .value = value};
            s.push(c);
        }
    }
    while (!s.empty()) {
        int index = s.top().index;
        result[index] = -1;
        s.pop();
    }
    for (int res: result) {
        cout << res << " ";
    }
    return 0;
}
