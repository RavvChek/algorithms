#include<iostream>

using namespace std;
int a[100][100], b[100];

void dfs(int j, int n, int *k) {
    if (b[j]) return;
    b[j] = 1;
    (*k)++;
    for (int i = 0; i < n; ++i)
        if (a[j][i])dfs(i, n, k);
}

int main() {
    int s, n, k;
    cin >> n >> s;
    s--;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; j++) {
            cin >> a[i][j];
        }
    }
    b[s] = 1;
    for (int j = 0; j < n; ++j) {
        if (a[s][j]) {
            dfs(j, n, &k);
        }
    }
    cout << k;
    return 0;
}