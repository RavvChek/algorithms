#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>

using namespace std;

struct task {
    int a_value;
    int b_value;
    int index;
};

bool comp_a(struct task x1, struct task x2) {
    if (x1.a_value == x2.a_value) {
        if (x1.b_value == x2.b_value) {
            return x1.index < x2.index;
        } else {
            return x1.b_value > x2.b_value;
        }
    } else {
        return x1.a_value > x2.a_value;
    }
}

bool comp_b(struct task x1, struct task x2) {
    if (x1.b_value == x2.b_value) {
        if (x1.a_value == x2.a_value) {
            return x1.index < x2.index;
        } else {
            return x1.a_value > x2.a_value;
        }
    } else {
        return x1.b_value > x2.b_value;
    }
}

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int n;
    cin >> n;
    int value;
    vector<int> a(n);
    vector<int> b(n);
    vector<int> mood(n);
    vector<struct task> a_cort(n);
    vector<struct task> b_cort(n);
    unordered_map<int, struct task> tasks;
    for (int i = 0; i < n; i++) {
        cin >> value;
        a[i] = value;
    }
    for (int i = 0; i < n; i++) {
        cin >> value;
        b[i] = value;
    }
    for (int i = 0; i < n; i++) {
        cin >> value;
        mood[i] = value;
        struct task t{.a_value = a[i], .b_value = b[i], .index = i};
        a_cort[i] = t;
        b_cort[i] = t;
        tasks.insert({i, t});
    }
    sort(a_cort.begin(), a_cort.end(), comp_a);
    sort(b_cort.begin(), b_cort.end(), comp_b);
    int a_id = 0;
    int b_id = 0;
        for (int ind: mood) {
        if (ind == 1) {
            if (tasks.count(b_cort[b_id].index) != 0) {
                cout << b_cort[b_id].index + 1 << " ";
                tasks.erase(b_cort[b_id].index);
                b_id++;
            } else {
                while (tasks.count(b_cort[b_id].index) == 0) {
                    b_id++;
                }
                cout << b_cort[b_id].index + 1 << " ";
                tasks.erase(b_cort[b_id].index);
                b_id++;
            }
        } else {
            if (tasks.count(a_cort[a_id].index) != 0) {
                cout << a_cort[a_id].index + 1 << " ";
                tasks.erase(a_cort[a_id].index);
                a_id++;
            } else {
                while (tasks.count(a_cort[a_id].index) == 0) {
                    a_id++;
                }
                cout << a_cort[a_id].index + 1 << " ";
                tasks.erase(a_cort[a_id].index);
                a_id++;
            }
        }
    }
    return 0;
}