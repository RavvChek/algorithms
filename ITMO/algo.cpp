// яндекс-контест ИТМО 1
#include <iostream>
#include <vector>
using namespace std;

int main()
{
    int count;
    cin >> count;
    vector<int> flowers(count);
    for (int i = 0; i < count; i++) {
        cin >> flowers[i];
    }

    int prevNum = -1;
    int prevPrevNum = -1;
    int maxLength = 0;
    int curLength = 0;
    int index = 0;
    int startIndex = 1;
    int maxStartIndex = 1;
    for (int flower : flowers) {
        if (flower != prevNum || flower != prevPrevNum) {
            curLength++;
        }
        else {
            if (maxLength < curLength) {
                maxLength = curLength;
                maxStartIndex = startIndex;
            }
            curLength = 2;
            startIndex = index;
        }
        prevPrevNum = prevNum;
        prevNum = flower;
        index++;
    }
    if (maxLength < curLength) {
        maxLength = curLength;
        maxStartIndex = startIndex;
    }
    /*cout << maxLength << endl;*/
    cout << maxStartIndex << " ";
    cout << maxStartIndex + maxLength - 1;
}



// яндекс-контест ИТМО 2
#include <iostream>
#include <string>
#include <stack>
using namespace std;

int main()
{
	string input ;
	cin >> input;

	const int size = input.size();
	const int n = size / 2;
	stack<char> letters;
	stack<int> traps_ids;
	stack<int> animals_ids;
	string result = "Possible\n";
	int count = -1;
	int* result_ids = new int[n];

	for (int i = 0; i < size; i++) {
		if (islower(input[i])) {
			animals_ids.push(i - count + 1);
		}
		else {
			count++;
			traps_ids.push(count + 1);
		}
		if (letters.empty() || input[i] == letters.top()) {
			letters.push(input[i]);
		}
		else if (toupper(input[i]) == toupper(letters.top())) {
			result_ids[traps_ids.top() - 1] = animals_ids.top()  - 1;
			letters.pop();
			traps_ids.pop();
			animals_ids.pop();
		}
		else {
			letters.push(input[i]);
		}
	}

	if (letters.empty()) {
		cout << result;
		for (int i = 0; i < n; i++) {
			cout << result_ids[i];
			if (i + 1!= n) {
				cout << " ";
			}
		}
	}
	else {
		cout << "Impossible";
	}
	return 0;
}



//яндекс-контест ИТМО 3
#include <iostream>
#include <string>
#include <unordered_map>
#include <stack>
#include <fstream>
#include <list>
using namespace std;

int main()
{
	unordered_map<string, stack<int>> variables;
	stack<list<string>> variablesStack;
	string str;
	list<string> l;
	variablesStack.push(l);
	ifstream input("input.txt");
	if (!input.is_open()) {
		return 1;
	}
	while (getline(input, str)) {
		if (str == "{") {
			list<string> arr;
			variablesStack.push(arr);
		} 
		else if (str == "}") {
			for (const auto& variable : variablesStack.top()) {
				variables[variable].pop();
			}
			variablesStack.pop();
		}
		else {
			int symbol = str.find('=');
			string var1 = str.substr(0, symbol);
			string var2 = str.substr(symbol + 1);

			if (isdigit(var2[0]) || var2[0] == '-') {
				variables[var1].push(stoi(var2));
			}
			else {
				if (variables.find(var2) == variables.end() || variables[var2].empty()) {
					variables[var2].push(0);
				}
				variables[var1].push(variables[var2].top());
				cout << variables[var1].top() << endl;
			}
			variablesStack.top().push_front(var1);
		}
	}
	input.close();
	return 0;
}


// яндекс-контест ИТМО 4
#include <iostream>
using namespace std;

int experiment(int a, int b, int c, int d, int k) {
    int l = a;
    if (k > 0) {
        k--;
        if (a * b <= c) {
            return 0;
        }
        a = max(min(a * b - c, d), 0);
        if (a == l) {
            return a;
        }
        return experiment(a, b, c, d, k);
    }
    else {
        return a;
    }
}

int main()
{
    int a, b, c, d, k;
    cin >> a >> b >> c >> d >> k;

    cout << experiment(a, b, c, d, k);
    return 0;
}



// яндекс-контест ИТМО 5
#include <iostream>
#include <vector>

int check_distance(int x, int count_cows, std::vector<int> stalls) {
	int cows = 1;
	int last_cow = stalls[0];
	int flag = 0;
	for (int i = 1; i < stalls.size(); i++) {
		if (stalls[i] - last_cow >= x) {
			last_cow = stalls[i];
			cows++;
		}
		if (cows == count_cows) {
			flag = 1;
			break;
		}
		
	}
	return flag;
}


int main()
{
	int count_stalls;
	int count_cows;
	std::cin >> count_stalls;
	std::cin >> count_cows;
	std::vector<int> stalls(count_stalls);
	for (int i = 0; i < count_stalls; i++) {
		std::cin >> stalls[i];
	}
	int l = 0;
	int r = stalls[stalls.size() - 1];
	while (r - l > 1) {
		int middle = (l + r) / 2;
		if (check_distance(middle, count_cows, stalls)) {
			l = middle;
		}
		else {
			r = middle;
		}
	}
	std::cout << l;
	return 0;
}



// яндекс-контест ИТМО 6
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <list>

bool comparator(std::string str1, std::string str2) {
	return str1 + str2 > str2 + str1;
}

int main()
{
	std::list<std::string> numbers;
	std::string str;
	while (std::cin >> str) {
		numbers.push_back(str);
	}
	numbers.sort(comparator);
	for (std::string num : numbers) {
		std::cout << num;
	}
	return 0;
}



// яндекс-контест ИТМО 7
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>

int main()
{
	std::cin.tie(0);
	std::ios_base::sync_with_stdio(0);
	std::string str;
	std::cin >> str;
	int count_of_letters = 26;
	int size_str = str.size();
	std::vector<int> weights(count_of_letters);
	std::unordered_map<char, int> letters_with_weights;
	std::vector<std::pair<char, int>> letters_in_string_sorted;
	std::string result1 = "";
	std::string result2 = "";
	int weight;
	for (int i = 0; i < count_of_letters; i++) {
		std::cin >> weight;
		weights[i] = weight;
	}
	for (int i = 0; i < count_of_letters; i++) {
		letters_with_weights['a' + i] = weights[i];
	}
	for (int i = 0; i < size_str; i++) {
		letters_in_string_sorted.push_back({str[i], letters_with_weights[str[i]]});
	}
	std::sort(letters_in_string_sorted.begin(), letters_in_string_sorted.end(), [](const std::pair<char, int>& a, const std::pair<char, int>& b) {
		if (a.second == b.second) {
			return a.first > b.first;
		}
		return a.second > b.second;
		});
	for (int i = 0; i < size_str - 1; i++) {
		if (letters_in_string_sorted[i].first == letters_in_string_sorted[i + 1].first && ((!result1.empty() && result1[result1.size() - 1] != letters_in_string_sorted[i].first) || result1.empty())) {
			result1 = result1 + letters_in_string_sorted[i].first;
			i++;
		}
		else {
			result2 = result2 + letters_in_string_sorted[i].first;
		}
	}
	if (size_str > result2.size() + 2 * result1.size()) {
		result2.push_back(letters_in_string_sorted[size_str - 1].first);
	}	
	std::cout << result1 + result2;
	for (int i = result1.size() - 1; i > -1; i--) {
		std::cout << result1[i];
	}
	return 0;
}




//яндекс-контест ИТМО 8
#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>

int main()
{
	int n;
	int k;
	std::cin >> n >> k;
	std::vector<int> prices(n);
	for (int i = 0; i < n; i++) {
		std::cin >> prices[i];
	}
	std::sort(prices.begin(), prices.end(), std::greater<>());
	for (int i = 1; i < n + 1; i++) {
		if (i % k == 0) {
			prices[i - 1] = 0;
		}
	}
	int sum = 0;
	for (int price : prices) {
		sum += price;
	}
	std::cout << sum;
	return 0;
}



// яндекс-контест ИТМО 9
#include <iostream>
#include <vector>
#include <unordered_set>
#include <queue>
#include <list>
#include <limits.h>

struct priority_car {
	int priority;
	int car;
	bool operator<(const priority_car& other) const {
		if (priority != other.priority) {
			return priority < other.priority;
		}
		else {
			return car > other.car;
		}
	}
};

int main()
{
	int N, K, P;
	std::cin >> N >> K >> P;
	std::unordered_set<int> cashe;
	std::priority_queue<struct priority_car> cars;
	std::vector<std::list<int>> priority(N + 1);
	std::vector<int> sequence(P);
	int count = 0;
	for (int i = 0; i < P; i++) {
		std::cin >> sequence[i];
		priority[sequence[i]].push_back(i);
	}
	for (int i = 0; i < P; i++) {
		int current_car = sequence[i];
		priority[current_car].pop_front();
		if (cashe.find(current_car) == cashe.end()) {
			if (cashe.size() >= K) {
				cashe.erase(cars.top().car);
				cars.pop();
			}
			count++;
			cashe.insert(current_car);

		}
		struct priority_car a;
		if (priority[current_car].empty()) {
			a.car = current_car;
			a.priority = INT_MAX;
			cars.push(a);
		}
		else {
			a.car = current_car;
			a.priority = priority[current_car].front();
			cars.push(a);
		}
	}
	std::cout << count;
	return 0;
}	



// яндекс-контест ИТМО 10
#include <iostream>
#include <list>
#include <string>
#include <vector>

int main()
{
	int N;
	std::cin >> N;
	std::list<int> queue1;
	std::list<int> queue2;
	std::vector<int> result;
	for (int i = 0; i < N; i++) {
		char symbol;
		int number;
		std::cin >> symbol;
		if (symbol == '-') {
			result.push_back(queue1.front());
			queue1.pop_front();
		}	
		else if (symbol == '+') {
			std::cin >> number;
			queue2.push_back(number);
		}
		else {
			std::cin >> number;
			queue1.push_back(number);
		}
		if (queue1.size() > queue2.size() + 1) {
			queue2.push_front(queue1.back());
			queue1.pop_back();
		}
		else if (queue1.size() < queue2.size()){
			queue1.push_back(queue2.front());
			queue2.pop_front();
		}
	}
	for (int i = 0; i < result.size(); i++) {
		std::cout << result[i] << std::endl;
	}
	return 0;
}



// яндекс-контест ИТМО 11
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>


int main() {
	int N, M, index_x,  size_x;
	std::cin >> N >> M;
	std::unordered_map<int, std::pair<int, int>> requests;
	std::multimap<int, int> free_blocks_by_size;
	std::map<int, int> free_blocks;
	free_blocks_by_size.insert({N, 1});
	free_blocks.insert({ 1, N });
	int req;
	std::vector<int> results;
	auto it_d = free_blocks.begin();
	for (int i = 1; i <= M; i++) {
		std::cin >> req;
		if (req <= 0) {
			requests.insert({ i, {req, 0} });
			std::pair<int, int> block = requests.at(abs(req));
			int index = block.second;
			int size = block.first;
			if (index == -1) {
				continue;
			}
			auto it_right = free_blocks.lower_bound(index);
			auto it_left = (it_right != free_blocks.begin()) ? std::prev(it_right) : free_blocks.end();
			
			if (it_right != free_blocks.end() && it_right->first == index + size) {
				if (it_left != free_blocks.end() && it_left->first + it_left->second == index) {
					index_x = it_left->first;
					size_x = it_left->second + it_right->second;
					it_d = free_blocks_by_size.find(it_left->second);
					while (it_d->second != it_left->first) it_d++;
					free_blocks_by_size.erase(it_d);
					free_blocks.erase(it_left);
					it_d = free_blocks_by_size.find(it_right->second);
					while (it_d->second != it_right->first) it_d++;
					free_blocks_by_size.erase(it_d);
					free_blocks.erase(it_right);
					free_blocks.insert({ index_x, size + size_x });
					free_blocks_by_size.insert({ size + size_x, index_x });
				}
				else {
					size_x = it_right->second;
					it_d = free_blocks_by_size.find(it_right->second);
					while (it_d->second != it_right->first) it_d++;
					free_blocks_by_size.erase(it_d);
					free_blocks.erase(it_right);
					free_blocks.insert({ index, size + size_x });
					free_blocks_by_size.insert({ size + size_x, index });
				}
			}
			else {
				if (it_left != free_blocks.end() && it_left->first + it_left->second == index) {
					index_x = it_left->first;
					size_x = it_left->second;
					it_d = free_blocks_by_size.find(it_left->second);
					while (it_d->second != it_left->first) it_d++;
					free_blocks_by_size.erase(it_d);
					free_blocks.erase(it_left);
					free_blocks.insert({index_x, size + size_x});
					free_blocks_by_size.insert({size + size_x, index_x});
				}
				else {
					free_blocks.insert({index, size});
					free_blocks_by_size.insert({size, index});
				}
			}
		}
		else {
			
			auto it = free_blocks_by_size.lower_bound(req);
			if (it == free_blocks_by_size.end()) {
				requests.insert({i, {req, -1}});
				results.push_back(-1);
			}
			else {
				int new_block_size = it->first - req;
				int index = it->second;
				requests.insert({ i, { req, index } });
				results.push_back(index);
				free_blocks.erase(it->second);
				free_blocks_by_size.erase(it);
				if (new_block_size > 0) {
					free_blocks.insert({index + req, new_block_size});
					free_blocks_by_size.insert({new_block_size, index + req});
				}
			}
		}
	}
	for (int val : results) {
		std::cout << val << std::endl;
	}
	return 0;
}




//яндекс-контест ИТМО 12
#include <iostream>
#include <vector>

int main()
{
	int N, K;
	std::cin >> N >> K;
	std::vector<int> numbers;
	for (int i = 0; i < N; i++) {
		int x;
		std::cin >> x;
		numbers.push_back(x);
	}
	int min;
	for (int i = 0; i < N - K + 1; i++) {
		min = 100000;
		int k = 0;
		int j = i;
		while (k < K) {
			if (numbers[j+k] <= min) {
				min = numbers[j+k];
			}
			k++;
		}
		std::cout << min << " ";
	}
}



// яндекс-контест ИТМО 13
#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <limits>
#include <algorithm>

enum Moves {
	DOWN = 'S',
	RIGHT = 'E',
	LEFT = 'W',
	UP = 'N'
};

struct CompareSecond {
	bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) {
		return a.second > b.second;
	}
};

int main()
{
	std::cin.tie(0);
	std::ios_base::sync_with_stdio(0);
	int INF = std::numeric_limits<int>::max();
	int N, M, x1, y1, x2, y2;
	std::cin >> N >> M >> x1 >> y1 >> x2 >> y2;
	std::string moves = "";
	std::string str;
	std::vector<std::string> grid(N);
	int start = (x1 - 1) * M + (y1 - 1);
	int goal = (x2 - 1) * M + (y2 - 1);
	std::vector<int> visited(N * M, INF);
	visited[start] = 0;
	for (int i = 0; i < N; i++) {
		std::cin >> grid[i];
	}
	std::vector<int> p (N * M);
	std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, CompareSecond> q;
	q.push({start, 0});
	std::vector<std::pair<int, int>> vertices;
	while (!q.empty()) {
		std::pair<int, int> v = q.top();
		x1 = v.first / M;
		y1 = v.first % M;
		q.pop();
		if (x1 == x2 && y1 == y2) {
			break;
		}
		if (grid[x1][y1] == '#') {
			continue;
		}
		if (x1 - 1 >= 0 && grid[x1 - 1][y1] != '#') {
			vertices.push_back({ v.first - M, (grid[x1 - 1][y1] == 'W') ? 2 : 1 });
		}	
		if (x1 + 1 < N && grid[x1 + 1][y1] != '#') {
			vertices.push_back({ v.first + M, (grid[x1 + 1][y1] == 'W') ? 2 : 1 });
		}
		if (y1 - 1 >= 0 && grid[x1][y1 - 1] != '#') {
			vertices.push_back({ v.first - 1, (grid[x1][y1 - 1] == 'W') ? 2 : 1 });
		}
		if (y1 + 1 < M && grid[x1][y1 + 1] != '#') {
			vertices.push_back({ v.first + 1, (grid[x1][y1 + 1] == 'W') ? 2 : 1 });
		}
		for (std::pair<int, int> val : vertices) {
			int new_len = val.second + v.second;
			if (visited[val.first] == INF || visited[val.first] > new_len) {
				q.push({val.first, new_len});
				visited[val.first] = new_len;
				p[val.first] = v.first;
			}
		}
		vertices.clear();
	}
	if (visited[goal] == INF) {
		std::cout << -1;
		return 0;
	}
	if (visited[goal] == 0) {
		std::cout << visited[goal];
	}
	else {
		std::cout << visited[goal] << std::endl;
	}
	for (int v = goal; v != start; v = p[v]) {
		p.push_back(v);
	}
	p.push_back(start);
	std::reverse(p.begin(), p.end());
	for (int i = 0; i < N * M; i++) {
		if (p[i] == goal) {
			break;
		}
		int diff = p[i + 1] - p[i];
		if (diff == -1) {
			moves += LEFT;
		}
		else if (diff == 1) {
			moves += RIGHT;
		}
		else if (diff == M) {
			moves += DOWN;
		}
		else if (diff == -M) {
			moves += UP;	
		}
	}
	std::cout << moves;
	return 0;
}




// яндекс-контест ИТМО 14
#include <iostream>
#include <vector>
#include <queue>

enum Color {
	WHITE = 0,
	BLACK = 1
};

void dfs( std::vector<std::vector<int>>& graph, std::vector<int>& visited, int v) {
	visited[v] = BLACK;
	for (int u : graph[v]) {
		if (visited[u] == WHITE) {
			dfs(graph, visited, u);
		}
	}
}

int main()
{
	int n;
	std::cin >> n;
	std::vector<std::vector<int>> graph(n + 1);
	std::vector<int> visited(n + 1, WHITE);
	int e;
	for (int v = 1; v < n + 1; ++v) {
		std::cin >> e;
		graph[v].push_back(e);
		graph[e].push_back(v);
	}
	int count = 0;
	for (int i = 1; i < graph.size(); i++) {
		if (visited[i] == WHITE) {
			count++;
			dfs(graph, visited, i);
		}
	}
	std::cout << count;
	return 0;
}



// яндекс-контест ИТМО 15
#include <iostream>
#include <vector>
#include <queue>

enum Color {
	WHITE,
	BLUE,
	RED
};

#define WHITE -1
#define BLUE 0
#define RED 1

bool is_bipartite(std::vector<std::vector<int>>& graph) {
	std::queue<int> q;
	std::vector<int> visited(graph.size(), WHITE);
	for (int s = 1; s < graph.size(); s++) {
		if (visited[s] == WHITE) {
			visited[s] = BLUE;
			q.push(s);
			while (!q.empty()) {
				int u = q.front();
				q.pop();
				for (int v : graph[u]) {
					if (visited[v] == -1) {
						visited[v] = 1 - visited[u];
						q.push(v);
					}
					else if (visited[u] == visited[v]) {
						return false;
					}
				}
			}
		}
	}
	return true;
}

int main()
{
	int N, M;
	std::cin >> N >> M;
	std::vector<std::vector<int>> graph(N + 1);
	int v, u;
	for (int i = 0; i < M; i++) {
		std::cin >> v >> u;
		graph[u].push_back(v);
		graph[v].push_back(u);
	}
	if (is_bipartite(graph)) {
		std::cout << "YES";
	}
	else {
		std::cout << "NO";
	}
	return 0;
}




//яндекс-контест ИТМО 16
#include <iostream>
#include <vector>
#include <limits.h>

struct vertex {
	int id;
	int weight;
};


void dfs1(std::vector<std::vector<vertex>>& graph, int val, int v, std::vector<int>& visited) {
	visited[v] = 1;
	for (int i = 0; i < graph.size(); ++i) {
		if (!visited[i] && graph[i][v].weight <= val) {
			dfs1(graph, val, i, visited);
		}
	}
}

void dfs2(std::vector<std::vector<vertex>>& graph, int val, int v, std::vector<int>& visited) {
	visited[v] = 1;
	for (int i = 0; i < graph.size(); ++i) {
		if (!visited[i] && graph[v][i].weight <= val) {
			dfs2(graph, val, i, visited);
		}
	}
}

bool check(std::vector<std::vector<vertex>>& graph, int val) {
	int size = graph.size();
	std::vector<int> visited (size, 0);
	dfs1(graph, val, 0, visited);
	for (int i = 0; i < size; ++i) {
		if (!visited[i]) {
			return false;
		}
	}
	std::fill(visited.begin(), visited.end(), 0);
	dfs2(graph, val, 0, visited);
	for (int i = 0; i < size; ++i) {
		if (!visited[i]) {
			return false;
		}
	}
	return true;
}

int bin_search(int max, std::vector<std::vector<vertex>>& graph) {
	int l = 0;
	int r = max;
	while (r - l > 1) {
		int m = (l + r) / 2;
		if (!check(graph, m)) {
			l = m;
		}
		else {
			r = m;
		}
	}
	return r;
}

int main()
{
	int n;
	std::cin >> n;
	vertex v;
	std::vector<std::vector<vertex>> graph(n);
	int oil;
	int max = INT_MIN;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cin >> oil;
			max = std::max(oil, max);
			vertex v;
			v.id = j;
			v.weight = oil;
			graph[i].push_back(v	);
		}
	}
	std::cout << bin_search(max, graph);
	return 0;
}




// тимус 1005
#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>
using namespace std;

int main()
{
	int size;
	cin >> size;

	vector<int> stones(size);
	for (int i = 0; i < size; i++) {
		cin >> stones[i];
	}

	int minDiff = INT_MAX;
	for (int code = 0; code < (1 << (size - 1)); code++) {
		int sum0 = 0;
		int sum1 = 0;
		for (int i = 0; i < size; i++) {
			if (((code >> i) & 1) == 0) {
				sum0 += stones[i];
			}
			else {
				sum1 += stones[i];
			}
		}
		minDiff = min(minDiff, abs(sum0 - sum1));
	}
	cout << minDiff;
	return 0;
}




// тимус 1155
#include <iostream>
#include <unordered_map>
using namespace std;

int main()
{
	unordered_map<char, int> cameras;
	for (int i = 0; i < 8; i++) {
		cin >> cameras['A' + i];
	}
	if (cameras['A'] + cameras['C'] + cameras['F'] + cameras['H'] != cameras['B'] + cameras['D'] + cameras['E'] + cameras['G']) {
		cout << "IMPOSSIBLE";
	}
	else {
		int flag = 0;
		while (true) {
			if (cameras['B'] > 0 && cameras['C'] > 0) {
				cout << "BC-" << endl;
				cameras['B']--;
				cameras['C']--;
			}
			if (cameras['B'] > 0 && cameras['A'] > 0) {
				cout << "BA-" << endl;
				cameras['B']--;
				cameras['A']--;
			}
			if (cameras['B'] > 0 && cameras['F'] > 0) {
				cout << "BF-" << endl;
				cameras['B']--;
				cameras['F']--;
			}
			if (cameras['B'] > 0 && cameras['H'] > 0) {
				cout << "CG+" << endl;
				cout << "HG-" << endl;
				cout << "BC-" << endl;
				cameras['B']--;
				cameras['H']--;
			}
			if (cameras['A'] > 0 && cameras['E'] > 0) {
				cout << "AE-" << endl;
				cameras['A']--;
				cameras['E']--;
			}
			if (cameras['A'] > 0 && cameras['D'] > 0) {
				cout << "AD-" << endl;
				cameras['A']--;
				cameras['D']--;
			}
			if (cameras['D'] > 0 && cameras['C'] > 0) {
				cout << "DC-" << endl;
				cameras['D']--;
				cameras['C']--;
			}
			if (cameras['D'] > 0 && cameras['H'] > 0) {
				cout << "DH-" << endl;
				cameras['D']--;
				cameras['H']--;
			}
			if (cameras['G'] > 0 && cameras['C'] > 0) {
				cout << "GC-" << endl;
				cameras['G']--;
				cameras['C']--;
			}
			if (cameras['G'] > 0 && cameras['H'] > 0) {
				cout << "GH-" << endl;
				cameras['G']--;
				cameras['H']--;
			}
			if (cameras['G'] > 0 && cameras['F'] > 0) {
				cout << "GF-" << endl;
				cameras['G']--;
				cameras['F']--;
			}
			if (cameras['E'] > 0 && cameras['F'] > 0) {
				cout << "EF-" << endl;
				cameras['E']--;
				cameras['F']--;
			}
			if (cameras['E'] > 0 && cameras['H'] > 0) {
				cout << "EH-" << endl;
				cameras['E']--;
				cameras['H']--;
			}
			if (cameras['D'] > 0 && cameras['F'] > 0) {
				cout << "AE+" << endl;
				cout << "DA-" << endl;
				cout << "FE-" << endl;
				cameras['D']--;
				cameras['F']--;
			}
			if (cameras['C'] > 0 && cameras['E'] > 0) {
				cout << "DH+" << endl;
				cout << "CD-" << endl;
				cout << "EH-" << endl;
				cameras['C']--;
				cameras['E']--;
			}
			if (cameras['A'] > 0 && cameras['G'] > 0) {
				cout << "BF+" << endl;
				cout << "AB-" << endl;
				cout << "GF-" << endl;
				cameras['A']--;
				cameras['G']--;
			}
			flag = 0;
			for (int i = 0; i < 8; i++) {
				if (cameras['A' + i] != 0) {
					flag = 1;
					break;
				}
			}
			if (flag == 0) {
				break;
			}
		}
	}
	return 0;
}




// тимус 1162
#include <iostream>
#include <list>
#include <vector>
using namespace std;

struct edge {
    int from, to;
    double rate, commission;
};

int main() {
    int N, M, S, A, B;
    double V, RAB, CAB, RBA, CBA;
    cin >> N >> M >> S >> V;
    list<edge> edges;
    S--;
    for (int i = 0; i < M; i++) {
        cin >> A >> B >> RAB >> CAB >> RBA >> CBA;
        A--; B--;
        edges.push_back({ A, B, RAB, CAB });
        edges.push_back({ B, A, RBA, CBA });
    }

    std::vector<double> currency(N, 0);
    currency[S] = V;

    for (int i = 0; i < N - 1; i++)
        for (edge x : edges)
            if ((currency[x.from] - x.commission) * x.rate > currency[x.to])
                currency[x.to] = (currency[x.from] - x.commission) * x.rate;

    for (edge x : edges)
        if ((currency[x.from] - x.commission) * x.rate > currency[x.to]) {
            cout << "YES";
            return 0;
        }
    cout << "NO";
    return 0;
}




// тимус 1444
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <cmath>

struct pumpkin {
	int x, y;
	int i;
};

int comp_angle(struct pumpkin p1, struct pumpkin p2) {
	return p1.x * p2.y - p2.x * p1.y;
}

int distance(const struct pumpkin p) {
	return p.x * p.x + p.y * p.y;
}

bool comparator(const struct pumpkin p1, const struct pumpkin p2) {
	return comp_angle(p1, p2) > 0 || (comp_angle(p1, p2) == 0 && distance(p1) < distance(p2));
}
int main()
{
	int N;
	std::cin >> N;
	std::vector<struct pumpkin> pumpkins1;
	std::vector<struct pumpkin> pumpkins2;
	int x;
	int y;
	struct pumpkin first_pumpkin;
	struct pumpkin pk;
	int start_x;
	int start_y;
	std::cin >> start_x;
	std::cin >> start_y;
	first_pumpkin.i = 1;
	first_pumpkin.x = start_x;
	first_pumpkin.y = start_y;
	for (int i = 1; i < N; i++) {
		std::cin >> x;
		std::cin >> y;
		pk.i = i + 1;
		pk.x = x;
		pk.y = y;
		pk.x -= start_x;
		pk.y -= start_y;
		if (pk.y > 0 || (pk.y == 0 && pk.x > 0)) {
			pumpkins1.push_back(pk);
		}
		else {
			pumpkins2.push_back(pk);
		}
	}
	std::sort(pumpkins1.begin(), pumpkins1.end(), comparator);
	std::sort(pumpkins2.begin(), pumpkins2.end(), comparator);
	std::cout << N << std::endl << first_pumpkin.i << std::endl;
 	if ((pumpkins1.empty() || pumpkins2.empty()) || comp_angle(pumpkins1.back(), pumpkins2.front()) > 0) {
		for (int i = 0; i < pumpkins1.size(); i++) {
			std::cout << pumpkins1[i].i << std::endl;
		}
		for (int i = 0; i < pumpkins2.size(); i++) {
			std::cout << pumpkins2[i].i << std::endl;
		}
	}
	else {
		for (int i = 0; i < pumpkins2.size(); i++) {
			std::cout << pumpkins2[i].i << std::endl;
		}
		for (int i = 0; i < pumpkins1.size(); i++) {
			std::cout << pumpkins1[i].i << std::endl;
		}
	}
	return 0;
}




// тимус 1450 
#include <iostream>
#include <vector>

using namespace std;


int main() {

    ios::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);

    int n, m, start, finish;
    cin >> n >> m;

    vector<vector<int> > edges(m, vector<int>(3));
    for (int i = 0; i < m; ++i) {
        cin >> edges[i][0] >> edges[i][1] >> edges[i][2];
        edges[i][0]--;
        edges[i][1]--;
    }
    cin >> start >> finish;
    vector<int> d(n, -1);

    d[--start] = 0;
    for (int i = 0; i < n - 1; ++i){
        for (int j = 0; j < m; ++j){
            if (d[edges[j][0]] > -1) {
                d[edges[j][1]] = max(d[edges[j][1]], d[edges[j][0]] + edges[j][2]);
            }
        }
    }
    if(d[--finish] == 1) {
        cout << "No solution";
        return 0;
    }
    cout << d[finish];

    return 0;
}



// тимус 1628
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <unordered_map>

int main()
{
	int m, n, k, count = 0;
	std::cin >> m >> n >> k;
	int x, y;
	std::set<std::pair<int, int>> single_elements;
	std::unordered_map<int, std::set<int>> rows;
	std::unordered_map<int, std::set<int>> cols;
	for (int i = 0; i < k; i++) {
		std::cin >> x >> y;
		rows[x].insert(y);
		cols[y].insert(x);
	}
	for (int i = 1; i <= m; i++) {
		int prev = 0;
		for (int cur : rows[i]) {
			if (cur - prev > 2) {
				count++;
			}
			else if (cur - prev == 2){
				single_elements.insert({ i, cur - 1});
			}
			prev = cur;
		}
		if (n - prev > 1) {
			count++;
		}
		else if (n - prev == 1){
			single_elements.insert({ i, n } );
		}
	}
	for (int i = 1; i <= n; i++ ) {
		int prev = 0;
		for (int cur : cols[i]) {
			if (cur - prev > 2) {
				count++;
			}
			else if (cur - prev == 2){
				if (single_elements.find({cur - 1, i}) != single_elements.end()) {
					count++;
				}
			}
			prev = cur;
		}
		if (m - prev > 1) {
			count++;
		}
		else if (m - prev == 1){
			if (single_elements.find({m, i}) != single_elements.end()) {
				count++;
			}
		}
	}
	std::cout << count;
	return 0;
}




// тимус 1650 
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <set>


int main()
{
	std::cin.tie(0);
	std::ios_base::sync_with_stdio(0);
	int n;
	std::cin >> n;
	std::unordered_map<std::string, std::pair<std::string, long long>> people;
	std::map<std::string, int> cities_days;
	std::map<long long, std::set<std::string>> rich_city;
	std::unordered_map<std::string, long long> cities_sums;
	std::string name, city;
	long long money;
	for (int i = 0; i < n; i++) {
		std::cin >> name >> city >> money;
		people[name] = { city, money };
		if (cities_sums.find(city) != cities_sums.end()) {
			long long sum = cities_sums[city];
			rich_city[sum].erase(city);
			if (rich_city[sum].size() == 0) {
				rich_city.erase(sum);
			}
		}
		cities_sums[city] += money;
		rich_city[cities_sums[city]].insert(city);
	}
	int m, k, day;
	std::cin >> m >> k;
	int cur_day, prev_day = 0;
	for (int i = 0; i <= k; i++) {
		if (i == k) {
			day = m;
		}
		else {
			std::cin >> day >> name >> city;
		}
		cur_day = day;
		std::map<long long, std::set<std::string>>::reverse_iterator it = rich_city.rbegin();
		if (cur_day != prev_day && it->second.size() == 1) cities_days[*(it->second.begin())] += cur_day - prev_day;
		if (i < k) {
			std::string from_city = people[name].first;
			long long old_money = cities_sums[from_city];
			rich_city[old_money].erase(from_city);
			if (rich_city[old_money].size() == 0) {
				rich_city.erase(old_money);
			}
			money = people[name].second;
			cities_sums[from_city] -= money;
			rich_city[cities_sums[from_city]].insert(from_city);

			long long new_money = cities_sums[city];
			rich_city[new_money].erase(city);
			if (rich_city[new_money].size() == 0) {
				rich_city.erase(new_money);
			}
			cities_sums[city] += money;
			rich_city[cities_sums[city]].insert(city);
			std::move(people[name].first) = city;			
			prev_day = cur_day;
		}
	}
	for (auto it = cities_days.begin(); it != cities_days.end(); it++) {
		std::cout << it->first << " " << it->second << std::endl;
	}
	return 0;
}



// тимус 1726
#include <iostream>
#include <vector>
#include <algorithm>

int main()
{
	long long n;
	std::cin >> n;
	std::vector<int> X;
	std::vector<int> Y;
	int x;
	int y;
	for (int i = 0; i < n; i++) {
		std::cin >> x;
		std::cin >> y;
		X.push_back(x);
		Y.push_back(y);
	}
	std::sort(X.begin(), X.end());
	std::sort(Y.begin(), Y.end());
	long long s = 0;
	for (long long i = 1; i < n; i++) {
		s += ((X[i] - X[i - 1]) + (Y[i] - Y[i - 1])) * 2 * i * (n - i);
	}
	std::cout << s / (n * (n - 1)) << std::endl;
	return 0;
}
