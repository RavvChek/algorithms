#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <stack>
#include <algorithm>

using namespace std;

struct data {
    string str;
    int hash_before;
    int dots;
    int hash_after;
};

int comparator(struct data x, struct data y) {
    if (x.str == y.str && x.hash_after == y.hash_after && x.dots == y.dots && x.hash_before == y.hash_before) {
        return 1;
    } else {
        return 0;
    }
}

int main() {
    string str;
    ifstream in("input.txt");
    ofstream out("input.txt");
    getline(in, str);
    int n = stoi(str);
    vector<string> matrix;
    int top = n, bottom = -1, left = n, right = -1;
    for (int i = 0; i < n; i++) {
        getline(in, str);
        for (int j = 0; j < n; j++) {
            if (str[j] == '#') {
                top = min(top, i);
                bottom = max(bottom, i);
                left = min(left, j);
                right = max(right, j);
            }
        }
        matrix.push_back(str);
    }
    vector<string> grid;
    for (int i = top; i <= bottom; ++i) {
        grid.push_back(matrix[i].substr(left, right - left + 1));
    }
    vector<string> new_grid;
    stack<string> s;
    for (const string &val: grid) {
        if (!s.empty()) {
            if (s.top() == val) {
                continue;
            } else {
                s.push(val);
                new_grid.push_back(val);
            }
        } else {
            s.push(val);
            new_grid.push_back(val);
        }
    }
    vector<struct data> new_grid2;
    for (auto &i: new_grid) {
        stack<char> stack;
        int hash_before = 0;
        int dot_count = 0;
        int hash_after = 0;
        bool dot_found = false;
        for (char j: i) {
            if (j == '.') {
                dot_found = true;
                dot_count++;
            } else if (j == '#') {
                if (!dot_found) {
                    hash_before++;
                } else {
                    hash_after++;
                }
            }
        }
        struct data d;
        d.dots = dot_count;
        d.hash_after = hash_after;
        d.hash_before = hash_before;
        if (dot_count == 0) {
            d.str = "#";
        } else if (hash_before > 0 && hash_after > 0) {
            d.str = "#.#";
        } else if (hash_before > 0 && hash_after == 0) {
            d.str = "#.";
        }
        new_grid2.push_back(d);
    }
    vector<string> new_grid3;

    for (const struct data &val: new_grid2) {
        stack<struct data> s;
        if (!s.empty()) {
            if (comparator(s.top(), val)) {
                continue;
            } else {
                s.push(val);
                new_grid3.push_back(val.str);
            }
        } else {
            s.push(val);
            new_grid3.push_back(val.str);
        }
    }

    if (new_grid3.empty()) {
        out << "X";
        return 0;
    }
    if (new_grid3[0] == "#" && new_grid2.size() == 1) {
        out << "I";
    } else if (new_grid3[0] == "#" && new_grid3[1] == "#.#" && new_grid3[2] == "#" && new_grid3.size() == 3) {
        out << "O";
    } else if (new_grid3[0] == "#" && new_grid3[1] == "#." && new_grid3[2] == "#" && new_grid3.size() == 3) {
        out << "C";
    } else if (new_grid3[0] == "#." && new_grid3[1] == "#" && new_grid3.size() == 2) {
        out << "L";
    } else if (new_grid3[0] == "#.#" && new_grid3[1] == "#" && new_grid3[2] == "#.#" && new_grid2.size() == 3 &&
            comparator(new_grid2[0], new_grid2[2])) {
        out << "H";
    } else if (new_grid3[0] == "#" && new_grid3[1] == "#.#" && new_grid3[2] == "#" && new_grid3[3] == "#." &&
               new_grid3.size() == 4 && new_grid2[1].hash_before == new_grid2[3].hash_before) {
        out << "P";
    } else {
        out << "X";
    }
    in.close();
    out.close();
    return 0;
}
