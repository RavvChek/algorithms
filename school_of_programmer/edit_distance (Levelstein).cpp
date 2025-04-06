#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

int m(char a, char b) {
    return a == b ? 0 : 1;
}

int D(const string& word1, const string& word2, int M, int N, vector<vector<int>>& memo) {
    if (M == 0) return N;
    if (N == 0) return M;

    if (memo[M][N] != -1) {
        return memo[M][N];
    }

    int val1 = D(word1, word2, M, N - 1, memo) + 1;
    int val2 = D(word1, word2, M - 1, N, memo) + 1;
    int val3 = D(word1, word2, M - 1, N - 1, memo) + m(word1[M - 1], word2[N - 1]);

    int result = min(val1, val2);
    result = min(result, val3);

    memo[M][N] = result;
    return result;
}

class Solution {
public:
    int minDistance(string word1, string word2) {
        int M = word1.size();
        int N = word2.size();
        vector<vector<int>> memo(M + 1, vector<int>(N + 1, -1));
        return D(word1, word2, M, N, memo);
    }
};

int main() {
    Solution solution;
    string word1, word2;
    cin >> word1 >> word2;
    cout << solution.minDistance(word1, word2) << endl;
    return 0;
}
