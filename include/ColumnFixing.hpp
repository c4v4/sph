#ifndef AC_CFT_INCLUDE_COLUMNFIXING_HPP_
#define AC_CFT_INCLUDE_COLUMNFIXING_HPP_

#include "Greedy.hpp"
#include "Instance.hpp"
#include "cft.hpp"

#define Q_THRESHOLD (-0.001)

class ColumnFixing {
public:
    ColumnFixing(SubInstance& subinst_, Greedy& greedy_, MStar& covered_rows_) : subinst(subinst_), greedy(greedy_), covered_rows(covered_rows_) { }

    idx_t operator()(LocalMultipliers& u_star, GlobalSolution& S_star) {
        auto& cols = subinst.get_cols();
        auto nrows = subinst.get_nrows();

        std::sort(S_star.begin(), S_star.end());

        Q.clear();
        for (idx_t j = 0; j < cols.size(); ++j) {
            if (cols[j].compute_lagr_cost(u_star) < Q_THRESHOLD && std::binary_search(S_star.begin(), S_star.end(), subinst.get_global_col_idx(j))) {
                Q.emplace_back(j);
            }
        }

        covered_rows.reset_covered(cols, Q, nrows);

        columns_to_fix.clear();
        for (auto j : Q) {
            if (!covered_rows.is_redundant(cols[j])) { columns_to_fix.emplace_back(j); }
        }

        LocalSolution S(columns_to_fix);
        greedy(u_star, std::max<idx_t>(nrows / 200, 1UL), S);

        idx_t remaining_rows = subinst.fix_columns(S, u_star);

        IF_VERBOSE { fmt::print("Fixed {} columns.\n", S.size()); }

        IF_DEBUG {
            assert(!subinst.is_corrupted());
            for ([[maybe_unused]] auto& col : cols) { assert(!col.empty()); }
        }

        return remaining_rows;
    }

private:
    SubInstance& subinst;
    Greedy& greedy;
    MStar& covered_rows;

    std::vector<idx_t> Q;
    std::vector<idx_t> columns_to_fix;
};

#endif