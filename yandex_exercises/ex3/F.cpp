#include <iostream>
#include <string>
#include <stack>
#include <unordered_map>

using namespace std;

unordered_map<char, int> createWeightMap(const string &priorities) {
    unordered_map<char, int> weight;
    for (int i = 0; i < priorities.size(); ++i) {
        weight[priorities[i]] = i + 1;
    }
    return weight;
}

bool compareBrackets(char bracket1, char bracket2, const unordered_map<char, int> &weight) {
    return weight.at(bracket1) <= weight.at(bracket2);
}

int main() {
    int n;
    string str;
    string w, s;
    getline(cin, str);
    n = stoi(str);
    getline(cin, w);
    getline(cin, s);
    unordered_map<char, int> weightMap = createWeightMap(w);
    string start_string = s;
    stack<char> bracket_stack;
    for (char i: s) {
        if (!bracket_stack.empty()) {
            if (bracket_stack.top() == '[' && i == ']' ||
                bracket_stack.top() == '(' && i == ')') {
                bracket_stack.pop();
            } else {
                bracket_stack.push(i);
            }
        } else {
            bracket_stack.push(i);
        }
    }
    string end_string;
    while (!bracket_stack.empty()) {
        if (bracket_stack.top() == '(') {
            end_string += ')';
            bracket_stack.pop();
        } else {
            end_string += ']';
            bracket_stack.pop();
        }
    }
    int brackets_missing_count = n - end_string.size() - start_string.size();
    if (brackets_missing_count == 0) {
        cout << start_string + end_string;
        return 0;
    }
    int flag = 0;
    if (w[0] == '(') {
        flag = 3;
    } else if (w[0] == '[') {
        flag = 1;
    } else if (w[1] == '[' && w[0] == ']') {
        flag = 2;
    } else if (w[1] == '(' && w[0] == ')') {
        flag = 4;
    } else if (w[1] == '(') {
        flag = 3;
    } else if (w[1] == '[') {
        flag = 1;
    } else if (w[2] == '[') {
        flag = 2;
    } else if (w[2] == '(') {
        flag = 4;
    }
    string add_string;
    if (flag == 2) {
        for (int i = 0; i < brackets_missing_count / 2; i++) {
            add_string += "[]";
        }
    } else if (flag == 1) {
        for (int i = 0; i < brackets_missing_count / 2; i++) {
            add_string = '[' + add_string + ']';
        }
    } else if (flag == 3) {
        for (int i = 0; i < brackets_missing_count / 2; i++) {
            add_string = '(' + add_string + ')';
        }
    } else {
        for (int i = 0; i < brackets_missing_count / 2; i++) {
            add_string += "()";
        }
    }
    int index = 0;
    for (int i = 0; i < end_string.size(); i++) {
        if (compareBrackets(end_string[i], add_string[0], weightMap)) {
            start_string += end_string[i];
        } else {
            index = i;
            flag = 1;
            break;
        }
    }
    start_string += add_string;
    if (flag == 1) {
        for (int i = index; i < end_string.size(); i++) {
            start_string += end_string[i];
        }
    }
    cout << start_string << endl;
    return 0;
}


