#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

int calculate_count_children(const string &person, map<string, vector<string>> &pedigree,
                             map<string, int> &count_children) {
    if (count_children.find(person) != count_children.end()) {
        return count_children[person];
    }
    int count = 0;
    for (const string &child: pedigree[person]) {
        count += calculate_count_children(child, pedigree, count_children) + 1;
    }
    return count;
}

int main() {
    int N;
    cin >> N;
    vector<string> people;
    map<string, vector<string>> pedigree;
    map<string, int> count_children;
    for (int i = 0; i < N - 1; i++) {
        string child, parent;
        cin >> child >> parent;
        people.push_back(child);
        people.push_back(parent);
        pedigree[parent].push_back(child);
    }
    sort(people.begin(), people.end());
    people.erase(unique(people.begin(), people.end()), people.end());
    for (const string &person: people) {
        count_children[person] = calculate_count_children(person, pedigree, count_children);
    }
    for (pair<string, int> value: count_children) {
        cout << value.first << " " << value.second << endl;
    }
    return 0;
}
