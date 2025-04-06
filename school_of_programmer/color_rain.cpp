# include <iostream>
# include <vector>

using namespace std;

int main() {
    int N;
    cin >> N;
    vector<vector<int>> matrix(N, vector<int>(N));
    vector<int> colors(N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> matrix[i][j];
        }
    }
    for (int i = 0; i < N; i++) {
        cin >> colors[i];
    }
    int count = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (matrix[i][j] == 1 && colors[i] != colors[j]) {
                count++;
            }
        }
    }
    cout << count;
    return 0;
}