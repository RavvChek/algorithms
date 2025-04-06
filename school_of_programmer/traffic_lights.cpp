# include <iostream>
# include <map>

using namespace std;

int main() {
    int N, M;
    cin >> N >> M;
    map<int, int> tunnels;
    for (int i = 1; i <= N; i++) {
        tunnels[i] = 0;
    }
    int i, j;
    for (int k = 0; k < M; k++) {
        cin >> i >> j;
        tunnels[i] = ++tunnels[i];
        tunnels[j] = ++tunnels[j];
    }
    for (int k = 1; k <= N; k++) {
        cout << tunnels[k] << " ";
    }
    return 0;
}
