#include <iostream>
#include <string>
#include <stack>

using namespace std;

int main() {
    string expression;
    string str;
    getline(cin, expression);
    for (char s: expression) {
        if (s != ' ') {
            str += s;
        }
    }
    expression = str;
    stack<int> number_stack;
    for (char symbol: expression) {
        if (symbol == '+') {
            int num1 = number_stack.top();
            number_stack.pop();
            int num2 = number_stack.top();
            number_stack.pop();
            number_stack.push(num1 + num2);
        } else if (symbol == '*') {
            int num1 = number_stack.top();
            number_stack.pop();
            int num2 = number_stack.top();
            number_stack.pop();
            number_stack.push(num1 * num2);
        } else if (symbol == '-') {
            int num1 = number_stack.top();
            number_stack.pop();
            int num2 = number_stack.top();
            number_stack.pop();
            number_stack.push(num2 - num1);
        } else {
            int number = symbol - '0';
            number_stack.push(number);
        }
    }
    cout << number_stack.top();
    return 0;
}
