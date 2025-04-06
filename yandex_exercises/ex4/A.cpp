#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

int calculate_height(const string &person, map<string, string> &pedigree, map<string, int> &heights) {
    if (heights.find(person) != heights.end()) {
        return heights[person];
    }
    if (pedigree.find(person) == pedigree.end())
        heights[person] = 0;
    return heights[person] = calculate_height(pedigree[person], pedigree, heights) + 1;
}

int main() {
    int N;
    cin >> N;
    string str;
    vector<string> people;
    map<string, string> pedigree;
    map<string, int> heights;
    for (int i = 0; i < N - 1; i++) {
        string child, parent;
        cin >> child >> parent;
        people.push_back(child);
        people.push_back(parent);
        pedigree[child] = parent;
    }
    sort(people.begin(), people.end());
    people.erase(unique(people.begin(), people.end()), people.end());
    for (const string &person: people) {
        calculate_height(person, pedigree, heights);
    }
    for (pair<string, int> value: heights) {
        if (!value.first.empty()) {
            cout << value.first << " " << value.second - 2 << endl;
        }
    }
    return 0;
}
