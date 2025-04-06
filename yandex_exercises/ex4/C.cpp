#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>

using namespace std;

int calculate_height(const string &person, unordered_map<string, string> &pedigree, map<string, int> &heights) {
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
    string child, parent;
    vector<string> people;
    unordered_map<string, string> pedigree;
    map<string, int> heights;
    for (int i = 0; i < N - 1; i++) {
        cin >> child >> parent;
        people.push_back(child);
        people.push_back(parent);
        pedigree[child] = parent;
    }
    for (string person : people) {
        calculate_height(person, pedigree, heights);
    }
    string person1;
    string person2;
    while (cin >> person1 >> person2) {
        while (heights[person1] > heights[person2]) {
            person1 = pedigree[person1];
        }
        while (heights[person1] < heights[person2]) {
            person2 = pedigree[person2];
        }
        while (person1 != person2) {
            person1 = pedigree[person1];
            person2 = pedigree[person2];
        }
        cout << person1 << endl;
    }
    return 0;
}
