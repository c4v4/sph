#ifndef AC_CFT_INCLUDE_RENDUNDACYREMOVAL_HPP
#define AC_CFT_INCLUDE_RENDUNDACYREMOVAL_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

#include "Instance.hpp"
#include "LowerBound.hpp"
#include "MStar.hpp"
#include "Scores.hpp"
#include "Solution.hpp"
#include "cft.hpp"

#define ENUM_THRESH 0

template <typename T, size_t N>
constexpr std::array<T, N> make_array(T value) {
    std::array<T, N> a;
    for (auto& x : a) { x = value; }
    return a;
}


class RendundacyRemoval {
public:
    RendundacyRemoval(SubInstance& subinst_, MStar& M_star_) : subinst(subinst_), M_star(M_star_) { }

    void reset(const LocalSolution& S) {
        redundant_cols.clear();
        for (auto j : S) {
            if (M_star.is_redundant(subinst.get_col(j))) { redundant_cols.emplace_back(j); }
        }
        // assert(M_star.get_uncovered() == 0);
    }

    inline auto size() { return redundant_cols.size(); }

    auto _enumeration_removal(std::array<bool, ENUM_THRESH>& vars, idx_t end, real_t& UB, real_t partial_cost, const LocalMultipliers& u_k) {

        if (end == 0) { return std::make_pair(partial_cost, std::array<bool, ENUM_THRESH>(vars)); }

        idx_t last = end - 1;
        auto& col = subinst.get_col(redundant_cols[last]);

        auto sol_pair0 = std::make_pair(UB, std::array<bool, ENUM_THRESH>());
        if (M_star.is_redundant(col)) {

            vars[last] = false;
            for (auto i : col) { M_star.uncover(i); }
            assert(M_star.get_uncovered() == 0);

            sol_pair0 = _enumeration_removal(vars, end - 1, UB, partial_cost, u_k);
            if (UB > sol_pair0.first) { UB = sol_pair0.first; }

            vars[last] = true;
            for (auto i : col) { M_star.cover(i); }
        }

        if (sol_pair0.first == partial_cost) { return sol_pair0; }  // removed all from here, avoid computing LB

        // Compute a LB
        auto next_p_cost = partial_cost + col.get_cost();
        // auto LB1 = best_col_LB(subinst, M_star, redundant_cols.begin(), redundant_cols.begin() + end);
        if (next_p_cost + 1.0 > UB) { return sol_pair0; }

        auto sol_pair1 = _enumeration_removal(vars, end - 1, UB, next_p_cost, u_k);

        if (UB > sol_pair1.first) { UB = sol_pair1.first; }

        if (sol_pair0.first < sol_pair1.first) { return sol_pair0; }

        return sol_pair1;
    }

    auto heur_enum_removal(const LocalMultipliers& u_k) {
        // assert(M_star.get_uncovered() == 0);

        cols_to_remove.clear();
        if (redundant_cols.empty()) { return cols_to_remove; }

        std::sort(redundant_cols.begin(), redundant_cols.end(), [&](auto a, auto b) { return subinst.get_col(a).get_cost() < subinst.get_col(b).get_cost(); });

        //IF_VERBOSE { fmt::print("Found {} redundant cols {{{}}}, ", redundant_cols.size(), fmt::join(redundant_cols, ", ")); }

        if (redundant_cols.size() > ENUM_THRESH) {
            // the first one is free
            cols_to_remove.emplace_back(redundant_cols.back());
            for (auto i : subinst.get_col(redundant_cols.back())) { M_star.uncover(i); }
            redundant_cols.pop_back();

            auto& cols = subinst.get_cols();
            while (redundant_cols.size() > ENUM_THRESH) {
                auto j = redundant_cols.back();
                if (M_star.is_redundant(cols[j])) {
                    cols_to_remove.emplace_back(j);
                    for (auto i : subinst.get_col(j)) { M_star.uncover(i); }
                }
                redundant_cols.pop_back();
            }
        }

        if (!redundant_cols.empty()) {
            // start enumeration
            auto vars = make_array<bool, ENUM_THRESH>(true);

            real_t UB = std::numeric_limits<real_t>::max();  // enum first tries to remove all the column in order
            auto sol_pair = _enumeration_removal(vars, redundant_cols.size(), UB, 0.0, u_k);
            assert(sol_pair.first == UB);

            for (idx_t j = 0; j < redundant_cols.size(); ++j) {
                if (!sol_pair.second[j]) { cols_to_remove.emplace_back(redundant_cols[j]); }
            }
        }

        // IF_VERBOSE {
        //    real_t cost = std::accumulate(cols_to_remove.begin(), cols_to_remove.end(), 0.0,
        //                                  [&](auto sum, auto idx) { return sum + subinst.get_col(idx).get_cost(); });
        //    fmt::print("removed {} cols, cost improvement of {}: {{{}}}\n", cols_to_remove.size(), cost, fmt::join(cols_to_remove, ", "));
        //}

        return cols_to_remove;
    }

    void operator()(LocalSolution& S, const LocalMultipliers& u_k) {
        reset(S);
        cols_to_remove = heur_enum_removal(u_k);
        S.remove(cols_to_remove);
    }

private:
    SubInstance& subinst;
    MStar& M_star;
    std::vector<idx_t> redundant_cols;
    std::vector<idx_t> cols_to_remove;
};

#endif