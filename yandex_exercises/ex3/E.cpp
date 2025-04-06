#include <iostream>
#include <stack>
#include <vector>
#include <string>
#include <regex>

using namespace std;

bool validateSpaces(const string &expression) {
    if (expression.size() == 1) {
        return true;
    }
    for (int i = 0; i < expression.size() - 3; i++) {
        if (isdigit(expression[i]) && expression[i + 1] == ' ' && isdigit(expression[i + 2])) {
            return false;
        }
        if (expression[i] == '+' && expression[i + 1] == ' ' && expression[i + 2] == '+') {
            return false;
        }
        if (expression[i] == '(' && (expression[i + 1] == '*' || expression[i + 1] == '+')) {
            return false;
        }
        if (expression[i + 1] == ')' && (expression[i] == '*' || expression[i] == '+' || expression[i] == '-')) {
            return false;
        }
    }
    return true;
}

vector<string> tokenizeExpression(const string &expression) {
    regex pattern(R"((\d+|\+|\-|\*|\/|\(|\)))");
    sregex_token_iterator it(expression.begin(), expression.end(), pattern);
    sregex_token_iterator end;
    vector<string> tokens;
    for (; it != end; ++it) {
        tokens.push_back(it->str());
    }
    return tokens;
}

bool validateSymbols(const string &expression) {
    for (char i: expression) {
        if (!(isdigit(i) || i == '+' || i == '-' || i == '*' ||
              i == '(' || i == ')' || i == ' ')) {
            return false;
        }
    }
    return true;
}

vector<string> infix_to_postfix_notation(const vector<string> &infix_notation) {
    vector<string> postfix_notation;
    stack<string> operand_stack;
    for (const string &symbol: infix_notation) {
        if (symbol == "+") {
            while (!operand_stack.empty() &&
                   (operand_stack.top() == "+" || operand_stack.top() == "-" || operand_stack.top() == "*")) {
                postfix_notation.push_back(operand_stack.top());
                operand_stack.pop();
            }
            operand_stack.push(symbol);
        } else if (symbol == "-") {
            while (!operand_stack.empty() &&
                   (operand_stack.top() == "+" || operand_stack.top() == "-" || operand_stack.top() == "*")) {
                postfix_notation.push_back(operand_stack.top());
                operand_stack.pop();
            }
            operand_stack.push(symbol);
        } else if (symbol == "*") {
            while (!operand_stack.empty() && operand_stack.top() == "*") {
                postfix_notation.push_back(operand_stack.top());
                operand_stack.pop();
            }
            operand_stack.push(symbol);
        } else if (symbol == "(") {
            operand_stack.push(symbol);
        } else if (symbol == ")") {
            int flag = 0;
            while (!operand_stack.empty() && operand_stack.top() != "(") {
                postfix_notation.push_back(operand_stack.top());
                operand_stack.pop();
                if (!operand_stack.empty() && operand_stack.top() == "(") {
                    flag = 1;
                }
            }
            if (!operand_stack.empty() && operand_stack.top() == "(") {
                flag = 1;
            }
            if (!operand_stack.empty()) {
                operand_stack.pop();
            }
            if (!flag) {
                postfix_notation.clear();
                return postfix_notation;
            }
        } else {
            if (symbol == " ") {
                continue;
            }
            postfix_notation.push_back(symbol);
        }
    }
    while (!operand_stack.empty()) {
        if (operand_stack.top() == "(") {
            postfix_notation.clear();
            return postfix_notation;
        }
        postfix_notation.push_back(operand_stack.top());
        operand_stack.pop();
    }
    return postfix_notation;
}

int solve_expression(const vector<string> &postfix_notation) {
    stack<int> number_stack;
    for (const string &symbol: postfix_notation) {
        if (symbol == "+") {
            int num1 = number_stack.top();
            number_stack.pop();
            int num2 = number_stack.top();
            number_stack.pop();
            number_stack.push(num1 + num2);
        } else if (symbol == "*") {
            int num1 = number_stack.top();
            number_stack.pop();
            int num2 = number_stack.top();
            number_stack.pop();
            number_stack.push(num1 * num2);
        } else if (symbol == "-") {
            int num1 = number_stack.top();
            number_stack.pop();
            int num2 = number_stack.top();
            number_stack.pop();
            number_stack.push(num2 - num1);
        } else {
            int number = stoi(symbol);
            number_stack.push(number);
        }
    }
    return number_stack.top();
}

int main() {
    string expression;
    getline(cin, expression);
    if (!validateSymbols(expression) || !validateSpaces(expression)) {
        cout << "WRONG";
        return 0;
    }
    string str;
    for (char symbol: expression) {
        if (symbol != ' ') {
            str += symbol;
        }
    }
    expression = str;
    str = "";
    if (expression[0] == '-' || expression[0] == '+') {
        str += '0' + expression[0];
    }
    for (int i = 0; i < expression.size() - 1; i++) {
        if (expression[i] == '(' && (expression[i + 1] == '-' || expression[i + 1] == '+')) {
            str += std::to_string(expression[i] + '0' + expression[i + 1]);
        } else {
            str += expression[i];
        }
    }
    if (expression[expression.size() - 1] == '-' || expression[expression.size() - 1] == '+' || expression[expression.size() - 1] == '*') {
        cout << "WRONG";
        return 0;
    }
    vector<string> postfix = infix_to_postfix_notation(tokenizeExpression(expression));
    if (postfix.empty()) {
        cout << "WRONG";
        return 0;
    }
    int result = solve_expression(postfix);
    cout << result;
    return 0;
}
