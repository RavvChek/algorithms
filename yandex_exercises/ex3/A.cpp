#include <iostream>
#include <string>
#include <stack>

using namespace std;

int main() {
    string str;
    getline(cin, str);
    stack<char> s;
    for (char bracket: str) {
        if (!s.empty()) {
            if (s.top() == '{' && bracket == '}' || s.top() == '[' && bracket == ']' ||
                s.top() == '(' && bracket == ')') {
                s.pop();
            } else if (bracket == '{' || bracket == '(' || bracket == '[') {
                s.push(bracket);
            } else {
                cout << "no";
                return 0;
            }
        } else {
            if (bracket == ']' || bracket == ')' || bracket == '}') {
                cout << "no";
                return 0;
            } else {
                s.push(bracket);
            }
        }
    }
    if (!s.empty()) {
        cout << "no";
    } else {
        cout << "yes";
    }
    return 0;
}
