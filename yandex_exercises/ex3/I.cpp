#include <iostream>
#include <queue>
#include <vector>
#include <unordered_map>
#include <map>
#include <queue>

using namespace std;

bool isEmptyQueues(const vector<queue<int>> &queues) {
    if (queues[0].empty() && queues[1].empty() && queues[2].empty() && queues[3].empty()) {
        return true;
    }
    return false;
}

struct Rover {
    int index;
    int direction;
} typedef Rover;

int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    int n;
    cin >> n;
    int a, b;
    cin >> a >> b;
    int end_time = 1;
    unordered_map<int, vector<Rover>> map_dir;
    for (int i = 0; i < n; i++) {
        int d, t;
        cin >> d >> t;
        end_time = max(t, end_time);
        Rover r = {.index = i, .direction = d};
        map_dir[t].push_back(r);
    }
    map<int, int> result;
    int time = 1;
    queue<int> q1;
    queue<int> q2;
    queue<int> q3;
    queue<int> q4;
    vector<queue<int>> queues;
    if (max(a, b) - min(a, b) == 2) {
        queues.push_back(q1);
        queues.push_back(q2);
        queues.push_back(q3);
        queues.push_back(q4);
        if (a == 2 && b == 4 || a == 4 && b == 2) {
            do {
                for (Rover rover: map_dir[time]) {
                    switch (rover.direction) {
                        case 1:
                            queues[0].push(rover.index);
                            break;
                        case 2:
                            queues[1].push(rover.index);
                            break;
                        case 3:
                            queues[2].push(rover.index);
                            break;
                        case 4:
                            queues[3].push(rover.index);
                            break;
                    }
                }
                if (!queues[1].empty() && !queues[3].empty()) {
                    result[queues[1].front()] = time;
                    result[queues[3].front()] = time;
                    queues[1].pop();
                    queues[3].pop();
                } else if (!queues[1].empty() && queues[3].empty()) {
                    result[queues[1].front()] = time;
                    queues[1].pop();
                } else if (queues[1].empty() && !queues[3].empty()) {
                    result[queues[3].front()] = time;
                    queues[3].pop();
                } else {
                    if (!queues[0].empty() && !queues[2].empty()) {
                        result[queues[0].front()] = time;
                        result[queues[2].front()] = time;
                        queues[0].pop();
                        queues[2].pop();
                    } else if (!queues[0].empty() && queues[2].empty()) {
                        result[queues[0].front()] = time;
                        queues[0].pop();
                    } else if (queues[0].empty() && !queues[2].empty()) {
                        result[queues[2].front()] = time;
                        queues[2].pop();
                    }
                }
                time++;
            }
            while (time <= 200 || !isEmptyQueues(queues));
        } else {
            do {
                for (Rover rover: map_dir[time]) {
                    switch (rover.direction) {
                        case 1:
                            queues[0].push(rover.index);
                            break;
                        case 2:
                            queues[1].push(rover.index);
                            break;
                        case 3:
                            queues[2].push(rover.index);
                            break;
                        case 4:
                            queues[3].push(rover.index);
                            break;
                    }
                }
                if (!queues[0].empty() && !queues[2].empty()) {
                    result[queues[0].front()] = time;
                    result[queues[2].front()] = time;
                    queues[0].pop();
                    queues[2].pop();
                } else if (!queues[0].empty() && queues[2].empty()) {
                    result[queues[0].front()] = time;
                    queues[0].pop();
                } else if (queues[0].empty() && !queues[2].empty()) {
                    result[queues[2].front()] = time;
                    queues[2].pop();
                } else {
                    if (!queues[1].empty() && !queues[3].empty()) {
                        result[queues[1].front()] = time;
                        result[queues[3].front()] = time;
                        queues[1].pop();
                        queues[3].pop();
                    } else if (!queues[1].empty() && queues[3].empty()) {
                        result[queues[1].front()] = time;
                        queues[1].pop();
                    } else if (queues[1].empty() && !queues[3].empty()) {
                        result[queues[3].front()] = time;
                        queues[3].pop();
                    }
                }
                time++;
            }
            while (time <= 200 || !isEmptyQueues(queues));
        }
    } else {
        if (a == 1 && b == 2 || a == 2 && b == 1) {
            queues.push_back(q1);
            queues.push_back(q2);
            queues.push_back(q3);
            queues.push_back(q4);
            do {
                for (Rover rover: map_dir[time]) {
                    switch (rover.direction) {
                        case 1:
                            queues[0].push(rover.index);
                            break;
                        case 2:
                            queues[1].push(rover.index);
                            break;
                        case 3:
                            queues[2].push(rover.index);
                            break;
                        case 4:
                            queues[3].push(rover.index);
                            break;
                    }
                }
                for (queue<int> &q: queues) {
                    if (!q.empty()) {
                        result[q.front()] = time;
                        q.pop();
                        break;
                    }
                }
                time++;
            } while (time <= 200 || !isEmptyQueues(queues));
        } else if (a == 2 && b == 3 || a == 3 && b == 2) {
            queues.push_back(q2);
            queues.push_back(q3);
            queues.push_back(q4);
            queues.push_back(q1);
            do {
                for (Rover rover: map_dir[time]) {
                    switch (rover.direction) {
                        case 1:
                            queues[3].push(rover.index);
                            break;
                        case 2:
                            queues[0].push(rover.index);
                            break;
                        case 3:
                            queues[1].push(rover.index);
                            break;
                        case 4:
                            queues[2].push(rover.index);
                            break;
                    }
                }
                for (queue<int> &q: queues) {
                    if (!q.empty()) {
                        result[q.front()] = time;
                        q.pop();
                        break;
                    }
                }
                time++;
            } while (time <= 200 || !isEmptyQueues(queues));
        } else if (a == 3 && b == 4 || a == 4 && b == 3) {
            queues.push_back(q3);
            queues.push_back(q4);
            queues.push_back(q1);
            queues.push_back(q2);
            do {
                for (Rover rover: map_dir[time]) {
                    switch (rover.direction) {
                        case 1:
                            queues[2].push(rover.index);
                            break;
                        case 2:
                            queues[3].push(rover.index);
                            break;
                        case 3:
                            queues[0].push(rover.index);
                            break;
                        case 4:
                            queues[1].push(rover.index);
                            break;
                    }
                }
                for (queue<int> &q: queues) {
                    if (!q.empty()) {
                        result[q.front()] = time;
                        q.pop();
                        break;
                    }
                }
                time++;
            } while (time <= 200 || !isEmptyQueues(queues));
        } else {
            queues.push_back(q4);
            queues.push_back(q1);
            queues.push_back(q2);
            queues.push_back(q3);
            do {
                for (Rover rover: map_dir[time]) {
                    switch (rover.direction) {
                        case 1:
                            queues[1].push(rover.index);
                            break;
                        case 2:
                            queues[2].push(rover.index);
                            break;
                        case 3:
                            queues[3].push(rover.index);
                            break;
                        case 4:
                            queues[0].push(rover.index);
                            break;
                    }
                }
                for (queue<int> &q: queues) {
                    if (!q.empty()) {
                        result[q.front()] = time;
                        q.pop();
                        break;
                    }
                }
                time++;
            } while (time <= 200 || !isEmptyQueues(queues));
        }
    }
    for (pair val: result) {
        cout << val.second << endl;
    }
    return 0;
}
