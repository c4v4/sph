// Copyright (c) 2022 Francesco Cavaliere
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef SPH_INCLUDE_PARSING_HPP_
#define SPH_INCLUDE_PARSING_HPP_

#include <algorithm>
#include <cassert>
#include <fstream>
#include <sstream>

#include "cft.hpp"

struct InstanceData {
    unsigned nrows{};
    std::vector<double> costs{};
    std::vector<double> solcosts{};
    std::vector<sph::idx_t> matbeg{};
    std::vector<sph::idx_t> matval{};
    std::vector<sph::idx_t> warmstart;
};

// trim from start (in place)
static inline void ltrim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
}

// trim from end (in place)
static inline void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string& s) {
    ltrim(s);
    rtrim(s);
}

std::vector<std::string> split(std::string& s, char delim) {
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

    InstanceData inst;
    inst.nrows = std::stoul(tokens[0]);
    const auto ncols = std::stoul(tokens[1]);
    inst.costs = std::vector<double>(ncols);
    inst.solcosts = std::vector<double>(ncols, 10e10);

    auto j = 0UL;
    while (j < inst.costs.size()) {
        std::getline(in, line);
        tokens = split(line, ' ');
        for (const auto& token : tokens) {
            assert(j < inst.costs.size());
            inst.costs[j] = std::stof(token);
            j++;
        }
    }
    assert(j == inst.costs.size());

    // for each row, the number of columns which cover row i followed by a list of the columns which cover row i
    auto cols = std::vector<std::vector<unsigned long>>(ncols);
    for (auto i = 0UL; i < inst.nrows; i++) {
        std::getline(in, line);
        const auto icols = std::stoul(line);

        auto n = 0UL;
        while (n < icols) {
            std::getline(in, line);
            tokens = split(line, ' ');
            for (const auto& token : tokens) {
                const auto cidx = std::stoul(token) - 1;
                assert(cidx < cols.size());
                cols[cidx].emplace_back(i);
                n++;
            }
        }

        assert(n == icols);
    }

    for (const auto& col : cols) {
        inst.matbeg.emplace_back(inst.matval.size());
        for (const auto i : col) {
            inst.matval.emplace_back(i);
        }
    }

    inst.matbeg.emplace_back(inst.matval.size());

    return inst;
}

InstanceData parse_rail_instance(const std::string& path) {

    auto in = std::ifstream(path);

    // rows columns
    auto line = std::string();
    std::getline(in, line);
    auto tokens = split(line, ' ');

    InstanceData inst;
    inst.nrows = std::stoul(tokens[0]);
    const auto ncols = std::stoul(tokens[1]);
    inst.costs = std::vector<double>(ncols);
    inst.solcosts = std::vector<double>(ncols, 10e10);

    for (auto j = 0UL; j < ncols; j++) {

        std::getline(in, line);
        tokens = split(line, ' ');
        inst.costs[j] = std::stof(tokens[0]);
        const auto jrows = std::stoul(tokens[1]);

        inst.matbeg.emplace_back(inst.matval.size());

        for (auto n = 0UL; n < jrows; n++) {
            const auto i = std::stoul(tokens[n + 2]) - 1;
            inst.matval.emplace_back(i);
        }
    }

    inst.matbeg.emplace_back(inst.matval.size());

    return inst;
}

InstanceData parse_cvrp_instance(const std::string& path) {

    auto in = std::ifstream(path);

    // rows columns
    auto line = std::string();
    std::getline(in, line);
    auto tokens = split(line, ' ');

    InstanceData inst;
    inst.nrows = std::stoul(tokens[0]);
    const auto ncols = std::stoul(tokens[1]);
    inst.costs = std::vector<double>(ncols);
    inst.solcosts = std::vector<double>(ncols);

    for (auto j = 0UL; j < ncols; j++) {

        std::getline(in, line);
        tokens = split(line, ' ');
        inst.costs[j] = std::stof(tokens[0]);
        inst.solcosts[j] = std::stof(tokens[1]);

        inst.matbeg.emplace_back(inst.matval.size());

        for (auto n = 2UL; n < tokens.size(); n++) {
            const auto i = std::stoul(tokens[n]);
            inst.matval.emplace_back(i);
        }
    }

    inst.matbeg.emplace_back(inst.matval.size());

    std::getline(in, line);
    tokens = split(line, ' ');
    for (const auto& v : tokens) {
        inst.warmstart.emplace_back(std::stoul(v));
    }

    return inst;
}

#endif  // SPH_INCLUDE_PARSING_HPP_
