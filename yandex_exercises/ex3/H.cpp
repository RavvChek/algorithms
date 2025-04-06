#include <iostream>
#include <vector>
#include <stack>

using namespace std;

int main() {
    int n;
    string str;
    cin >> str;
    n = stoi(str);
    stack<int> number_stack;
    vector<long long> prefix_sum;
    vector<long long> result;
    prefix_sum.push_back(0);
    for (int i = 0; i < n; i++) {
        cin >> str;
        if (str[0] == '+') {
            str.erase(0, 1);
            int num1 = stoi(str);
            long long num2 = prefix_sum[prefix_sum.size() - 1];
            number_stack.push(num1);
            prefix_sum.push_back(num2 + num1);
        } else if (str[0] == '?') {
            str.erase(0, 1);
            int num = stoi(str);
            result.push_back(prefix_sum[prefix_sum.size() - 1] - prefix_sum[prefix_sum.size() - 1 - num]);
        } else if (str[0] == '-') {
            result.push_back(number_stack.top());
            prefix_sum.pop_back();
            number_stack.pop();
        }
    }
    for (long long val: result) {
        cout << val << endl;
    }
    return 0;
}
