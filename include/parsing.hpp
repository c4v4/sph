#ifndef AC_CFT_INCLUDE_PARSING_HPP_
#define AC_CFT_INCLUDE_PARSING_HPP_

#include <algorithm>
#include <cassert>
#include <fstream>
#include <sstream>

struct InstanceData {
    idx_t nrows{};
    std::vector<real_t> costs{};
    std::vector<idx_t> matbeg{};
    std::vector<idx_t> matval{};
    std::vector<idx_t> warmstart;
};

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
      return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
      return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

std::vector<std::string> split(std::string &s, char delim) {
    trim(s);
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {

        elems.push_back(item);
    }
    return elems;
}

InstanceData parse_scp_instance(const std::string& path) {

    auto in = std::ifstream(path);

    // rows columns
    auto line = std::string();
    std::getline(in, line);
    auto tokens = split(line, ' ');
    const auto nrows = std::stoul(tokens[0]);
    const auto ncols = std::stoul(tokens[1]);

    // cost for each column
    auto costs = std::vector<real_t>(ncols);
    auto j = 0UL;
    while(j < costs.size()) {
        std::getline(in, line);
        tokens = split(line, ' ');
        for(const auto& token : tokens) {
            assert(j < costs.size());
            costs[j] = std::stof(token);
            j++;
        }
    }
    assert(j == costs.size());

    // for each row, the number of columns which cover row i followed by a list of the columns which cover row i
    auto cols = std::vector<std::vector<unsigned long>>(ncols);
    for(auto i = 0UL; i < nrows; i++) {
        std::getline(in, line);
        const auto icols = std::stoul(line);

        auto n = 0UL;
        while(n < icols) {
            std::getline(in, line);
            tokens = split(line, ' ');
            for(const auto& token : tokens) {
                const auto cidx = std::stoul(token) - 1;
                assert(cidx < cols.size());
                cols[cidx].emplace_back(i);
                n++;
            }
        }

        assert(n == icols);

    }

    auto matbeg = std::vector<idx_t>();
    auto matval = std::vector<idx_t>();

    for(const auto& col : cols) {
        matbeg.emplace_back(matval.size());
        for (const auto i : col) {
            matval.emplace_back(i);
        }
    }

    matbeg.emplace_back(matval.size());

    return {static_cast<idx_t>(nrows), costs, matbeg, matval, {}};

}

InstanceData parse_rail_instance(const std::string& path) {

    auto in = std::ifstream(path);

    // rows columns
    auto line = std::string();
    std::getline(in, line);
    auto tokens = split(line, ' ');
    const auto nrows = std::stoul(tokens[0]);
    const auto ncols = std::stoul(tokens[1]);

    auto costs = std::vector<real_t>(ncols);
    auto matbeg = std::vector<idx_t>();
    auto matval = std::vector<idx_t>();

    for(auto j = 0UL; j < ncols; j++) {

        std::getline(in, line);
        tokens = split(line, ' ');
        costs[j] = std::stof(tokens[0]);
        const auto jrows = std::stoul(tokens[1]);

        matbeg.emplace_back(matval.size());

        for(auto n = 0UL; n < jrows; n++) {
            const auto i = std::stoul(tokens[n+2]) - 1;
            matval.emplace_back(i);
        }

    }

    matbeg.emplace_back(matval.size());

    return {static_cast<idx_t>(nrows), costs, matbeg, matval, {}};

}

InstanceData parse_cvrp_instance(const std::string& path) {

    auto in = std::ifstream(path);

    // rows columns
    auto line = std::string();
    std::getline(in, line);
    auto tokens = split(line, ' ');
    const auto nrows = std::stoul(tokens[0]);
    const auto ncols = std::stoul(tokens[1]);

    auto costs = std::vector<real_t>(ncols);
    auto matbeg = std::vector<idx_t>();
    auto matval = std::vector<idx_t>();

    for(auto j = 0UL; j < ncols; j++) {

        std::getline(in, line);
        tokens = split(line, ' ');
        costs[j] = std::stof(tokens[0]);

        matbeg.emplace_back(matval.size());

        for(auto n = 1UL; n < tokens.size(); n++) {
            const auto i = std::stoul(tokens[n]);
            matval.emplace_back(i);
        }

    }

    matbeg.emplace_back(matval.size());

    auto warmstart = std::vector<idx_t>();
    std::getline(in, line);
    tokens = split(line, ' ');
    for(const auto& v : tokens) {
        warmstart.emplace_back(std::stoul(v));
    }

    return {static_cast<idx_t>(nrows), costs, matbeg, matval, warmstart};

}

#endif  // AC_CFT_INCLUDE_PARSING_HPP_
