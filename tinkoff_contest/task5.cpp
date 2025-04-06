#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <sstream>

long long timeToSeconds(const std::string &time_str) {
    int h, m, s;
    char colon;
    std::istringstream iss(time_str);
    iss >> h >> colon >> m >> colon >> s;
    return h * 3600 + m * 60 + s;
}

int main() {
    std::string start_time;
    std::cin >> start_time;
    long long start_seconds = timeToSeconds(start_time);
    int n;
    std::cin >> n;
    std::map<std::string, std::map<char, std::pair<bool, int>>> teams;
    for (int i = 0; i < n; ++i) {
        std::string team_name, query_time_str, server_id, result;
        std::cin.ignore();
        std::getline(std::cin, team_name, '\"');
        std::cin.ignore();
        std::cin >> query_time_str >> server_id >> result;
        long long query_seconds = timeToSeconds(query_time_str) - start_seconds;
        auto &team_info = teams[team_name];
        if (result == "ACCESSED") {
            team_info[server_id] = {true, query_seconds};
        } else if (result == "DENIED" || result == "FORBIDEN") {
            if (team_info[server_id].first == false) {
                team_info[server_id].second += 20;
            }
        }
    }
    std::vector<std::tuple<int, std::string, int, int>> results;
    for (const auto &team_pair : teams) {
        const std::string &team_name = team_pair.first;
        const auto &server_info = team_pair.second;
        int hacked_servers = 0;
        int penalty_time = 0;
        for (const auto &server_pair : server_info) {
            if (server_pair.second.first) {
                ++hacked_servers;
                penalty_time += server_pair.second.second;
            }
        }
        results.emplace_back(hacked_servers, team_name, penalty_time, team_name);
    }
    std::sort(results.begin(), results.end(), [](const auto &a, const auto &b) {
        if (std::get<0>(a) != std::get<0>(b)) {
            return std::get<0>(a) > std::get<0>(b);
        }
        if (std::get<2>(a) != std::get<2>(b)) {
            return std::get<2>(a) < std::get<2>(b);
        }
        return std::get<3>(a) < std::get<3>(b);
    });

    int rank = 1;
    for (const auto &result : results) {
        std::cout << rank++ << " \"" << std::get<1>(result) << "\" " << std::get<0>(result) << " " << std::get<2>(result) << std::endl;
    }

    return 0;
}
