#include <iostream>
#include <vector>

using namespace std;

int main() {
    int N;
    cin >> N;
    vector<vector<int>> matrix(N, std::vector<int>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> matrix[i][j];
        }
    }
    int min_distance = 1000000;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            for (int k = j + 1; k < N; k++) {
                int distance = matrix[i][j] + matrix[j][k] + matrix[k][i];
                if (distance < min_distance) {
                    min_distance = distance;
                }
            }
        }
    }
    cout << min_distance << endl;
    return 0;
}
