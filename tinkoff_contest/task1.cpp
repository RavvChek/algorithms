#include <iostream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;

int main() {
    string input;
    cin >> input;
    vector<int> result;
    stringstream ss(input);
    string part;
    while (getline(ss, part, ',')) {
        size_t dash_pos = part.find('-');
        if (dash_pos != string::npos) {
            int start = stoi(part.substr(0, dash_pos));
            int end = stoi(part.substr(dash_pos + 1));
            for (int i = start; i <= end; ++i) {
                result.push_back(i);
            }
        } else {
            result.push_back(stoi(part));
        }
    }
    for (auto it = result.begin(); it != result.end(); ++it) {
        if (it != result.begin()) cout << " ";
        cout << *it;
    }
    cout << endl;
    return 0;
}
