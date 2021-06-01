#ifndef AC_CFT_INCLUDE_INSTANCE_HPP_
#define AC_CFT_INCLUDE_INSTANCE_HPP_

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cassert>
#include <numeric>
#include <unordered_map>

#include "IndexList.hpp"
#include "MStar.hpp"
#include "TrivialHeap.hpp"
#include "cft.hpp"
#include "queue"

#define MIN_COV 5U

class Instance {
public:
    explicit Instance(const idx_t nrows_) : nrows(nrows_), active_rows(nrows) {
        std::iota(active_rows.begin(), active_rows.end(), 0);
    }

    void add_columns(const std::vector<real_t>& costs, const std::vector<idx_t>& matbeg, const std::vector<idx_t>& matval) {
        const auto ncols = costs.size();
        IF_VERBOSE { fmt::print("Adding {} columns.\n", ncols); }
        for (idx_t j = 0; j < ncols; ++j) {
            cols.emplace_back(costs[j], matval.begin() + matbeg[j], matval.begin() + matbeg[j + 1]);
            active_cols.emplace_back(j);
        }
        IF_VERBOSE { fmt::print("Problem size = {}x{}.\n", nrows, cols.size()); }
    }

    [[nodiscard]] inline auto get_ncols() const { return cols.size(); }
    [[nodiscard]] inline idx_t get_nrows() const { return nrows; }
    [[nodiscard]] inline const auto& get_active_rows() const { return active_rows; }


    [[nodiscard]] inline auto& get_active_cols() { return active_cols; }
    [[nodiscard]] inline auto& get_fixed_cols() { return fixed_cols; }
    [[nodiscard]] inline auto& get_cols() { return cols; }
    [[nodiscard]] inline auto& get_col(idx_t idx) { return cols[idx]; }
    [[nodiscard]] inline const auto& get_col(idx_t idx) const { return cols[idx]; }

    void reset_fixing() {

        assert(std::is_sorted(active_cols.begin(), active_cols.end()));
        assert(std::is_sorted(fixed_cols.begin(), fixed_cols.end()));

        // restore previously removed rows
        for (const auto& [i, col_idxs] : cols_of_rem_rows) {
            for (const auto j : col_idxs) { cols[j].emplace_back(i); }
        }

        active_rows.resize(nrows);
        std::iota(active_rows.begin(), active_rows.end(), 0);


        cols_of_rem_rows.clear();

        // merge active and fixed columns
        idx_t actlen = active_cols.size();
        active_cols.insert(active_cols.end(), fixed_cols.begin(), fixed_cols.end());
        fixed_cols.clear();

        std::inplace_merge(active_cols.begin(), active_cols.begin() + actlen, active_cols.end());
        assert(std::is_sorted(active_cols.begin(), active_cols.end()));
    }

    void fix_columns(const std::vector<idx_t>& idxs, const MStar& M_star) {

        assert(fixed_cols.empty());
        assert(active_rows.size() == nrows);
        assert(std::is_sorted(active_cols.begin(), active_cols.end()));
        assert(std::is_sorted(fixed_cols.begin(), fixed_cols.end()));

        idx_t actv_i = 0;
        for (idx_t i = 0; i < M_star.size(); ++i) {
            if (M_star[i] == 0) { active_rows[actv_i++] = i; }
        }
        active_rows.resize(actv_i);

        // remove rows on all columns
        for (idx_t j = 0; j < cols.size(); ++j) {
            auto& col = cols[j];
            for (idx_t n = 0; n < col.size();) {
                idx_t i = col[n];
                if (M_star[i] > 0) {
                    col[n] = col.back();
                    col.pop_back();
                    cols_of_rem_rows[i].emplace_back(j);
                } else {
                    ++n;
                }
            }
        }

        IF_VERBOSE { fmt::print("Instance: fixed {} columns and removed {} rows.\n", idxs.size(), cols_of_rem_rows.size()); };

        idx_t empty_cols = 0;

        // Prologue: untouched part
        idx_t iok = 0;
        while (idxs[0] != active_cols[iok] && !cols[active_cols[iok]].empty()) { ++iok; }

        // Body: remove idxs cols
        idx_t islide = iok;
        for (idx_t irem = 0; irem < idxs.size(); ++islide) {
            assert(islide < active_cols.size());

            idx_t j1 = idxs[irem];
            idx_t j2 = active_cols[islide];

            if (j1 != j2 && !cols[j2].empty()) {
                // This column is ok
                active_cols[iok++] = j2;
            } else {
                // Skip this column
                if (j1 == j2) {
                    ++irem;
                } else
                    IF_DEBUG { ++empty_cols; }
            }
        }

        // Epilogue: finish copying, excluding empty columns
        for (; islide < active_cols.size(); ++islide) {

            idx_t j = active_cols[islide];
            if (!cols[j].empty()) {
                active_cols[iok++] = j;
            } else
                IF_DEBUG { ++empty_cols; }
        }
        assert(iok == active_cols.size() - idxs.size() - empty_cols);

        active_cols.resize(iok);
        fixed_cols = idxs;
    }

    auto compute_fixed_cost() {
        real_t fixed_cost = 0.0;
        for (idx_t j : fixed_cols) { fixed_cost += cols[j].get_cost(); }
        return fixed_cost;
    }

private:
    const idx_t nrows;
    std::vector<Column> cols;
    std::vector<idx_t> active_cols;
    std::vector<idx_t> fixed_cols;

    std::unordered_map<idx_t, std::vector<idx_t>> cols_of_rem_rows;
    std::vector<idx_t> active_rows;
};

class SubInstance {

public:
    explicit SubInstance(Instance& inst_) : inst(inst_) { }

    [[nodiscard]] inline auto get_ncols() const { return cols.size(); }
    [[nodiscard]] inline auto get_nrows() const { return rows.size(); }

    [[nodiscard]] inline auto get_global_col_idx(idx_t local_j) const { return local_to_global_col_idxs[local_j]; }
    [[nodiscard]] inline auto get_global_row_idx(idx_t local_i) const { return local_to_global_row_idxs[local_i]; }
    [[nodiscard]] inline auto get_local_row_idx(idx_t global_i) const { return global_to_local_row_idxs[global_i]; }

    [[nodiscard]] inline auto& get_cols() { return cols; }
    [[nodiscard]] inline auto& get_rows() { return rows; }

    [[nodiscard]] inline const auto& get_cols() const { return cols; }
    [[nodiscard]] inline const auto& get_rows() const { return rows; }

    [[nodiscard]] inline auto& get_col(idx_t idx) { return cols[idx]; }
    [[nodiscard]] inline auto& get_row(idx_t idx) { return rows[idx]; }

    [[nodiscard]] inline const auto& get_col(idx_t idx) const { return cols[idx]; }
    [[nodiscard]] inline const auto& get_row(idx_t idx) const { return rows[idx]; }

    [[nodiscard]] inline auto& get_fixed_cols() { return fixed_cols_global_idxs; }

    [[nodiscard]] inline auto& get_instance() { return inst; }

    [[nodiscard]] inline auto compute_fixed_cost() const {
        real_t fixed_cost = inst.compute_fixed_cost();
        for (idx_t j : fixed_cols_global_idxs) {
            assert(std::find(inst.get_active_cols().begin(), inst.get_active_cols().end(), j) != inst.get_active_cols().end());
            fixed_cost += inst.get_col(j).get_cost();
        }
        return fixed_cost;
    }

    [[nodiscard]] inline auto is_corrupted() const {
        for (idx_t j = 0; j < cols.size(); ++j) {
            assert(!cols[j].empty());
            for (idx_t i : cols[j]) {
                if (std::find(rows[i].begin(), rows[i].end(), j) == rows[i].end()) {
                    IF_DEBUG { fmt::print("Col {} not found in row {}. \n Row: ", j, i, fmt::join(rows[i], ", ")); }
                    return true;
                }
            }
        }

        for (idx_t i = 0; i < rows.size(); ++i) {
            assert(!rows[i].empty());
            for (idx_t j : rows[i]) {
                if (std::find(cols[j].begin(), cols[j].end(), i) == cols[j].end()) {
                    IF_DEBUG { fmt::print("Row {} not found in col {}. \n Col: ", i, j, fmt::join(cols[j], ", ")); }
                    return true;
                }
            }
        }

        return false;
    }

    void reset() {
        idx_t active_nrows = inst.get_active_rows().size();

        // setup rows
        local_to_global_row_idxs.resize(active_nrows);
        global_to_local_row_idxs.resize(inst.get_nrows());

        idx_t p = 0;
        for (idx_t i : inst.get_active_rows()) {
            local_to_global_row_idxs[p] = i;
            global_to_local_row_idxs[i] = p;
            p++;
        }

        fixed_cols_global_idxs.clear();
        tot_removed_rows = 0;

        auto priced_cols = _cost_active_cols();
        covering_times.reset_uncovered(inst.get_nrows());
        auto global_col_idxs = _select_C2_cols(priced_cols, covering_times);

        // fill cols and rows
        replace_columns(global_col_idxs);

        IF_VERBOSE { fmt::print("Sub-instance size = {}x{}.\n", rows.size(), cols.size()); }

        assert(!is_corrupted());
    }

    real_t price(const std::vector<real_t>& u_k) {

        auto priced_cols = std::vector<PricedCol>();
        real_t global_LB = _price_active_cols(u_k, priced_cols);

        auto global_C1_idxs = _select_C1_cols(priced_cols, covering_times);
        auto global_C2_idxs = _select_C2_cols(priced_cols, covering_times);

        global_C1_idxs.insert(global_C1_idxs.end(), global_C2_idxs.begin(), global_C2_idxs.end());
        replace_columns(global_C1_idxs);

        IF_DEBUG {
            covering_times.reset_covered(cols, rows.size());
            assert(covering_times.get_uncovered() == 0);
            assert(!is_corrupted());
        }

        return global_LB;
    }

    idx_t fix_columns(const std::vector<idx_t>& local_idxs_to_fix, std::vector<real_t>& u_star) {

        if (local_idxs_to_fix.empty()) { return rows.size(); }
        idx_t removed = 0;

        // mark rows to remove
        for (idx_t lj : local_idxs_to_fix) {
            idx_t gj = local_to_global_col_idxs[lj];
            local_to_global_col_idxs[lj] = REMOVED_INDEX;
            assert(gj != REMOVED_INDEX);
            assert(std::find(fixed_cols_global_idxs.begin(), fixed_cols_global_idxs.end(), gj) == fixed_cols_global_idxs.end());

            fixed_cols_global_idxs.emplace_back(gj);
            const auto& col = cols[lj];
            for (idx_t li : col) {
                if (!_is_local_row_active(li)) { continue; }
                removed++;
                tot_removed_rows++;
                idx_t gi = local_to_global_row_idxs[li];
                local_to_global_row_idxs[li] = REMOVED_INDEX;
                global_to_local_row_idxs[gi] = REMOVED_INDEX;
            }
        }

        // compact rows
        idx_t li = 0, rows_left = 0;
        while (_is_local_row_active(li)) { ++li, ++rows_left; }

        for (; li < rows.size(); ++li) {
            if (_is_local_row_active(li)) {
                idx_t gs = local_to_global_row_idxs[li];

                assert(_is_global_row_active(gs));
                assert(!rows[rows_left].empty());

                global_to_local_row_idxs[gs] = rows_left;
                local_to_global_row_idxs[rows_left] = gs;
                u_star[rows_left] = u_star[li];
                ++rows_left;
            }
        }
        assert(rows_left == local_to_global_row_idxs.size() - removed);

        local_to_global_row_idxs.resize(rows_left);
        u_star.resize(rows_left);

        IF_DEBUG {
            assert(rows_left == inst.get_active_rows().size() - tot_removed_rows);
            for ([[maybe_unused]] idx_t gi : local_to_global_row_idxs) { assert(_is_global_row_active(gi)); }
        }

        if (rows_left > 0) {
            auto global_idxs = std::vector<idx_t>();
            global_idxs.reserve(cols.size() - local_idxs_to_fix.size());
            for (idx_t gj : local_to_global_col_idxs) {
                if (gj != REMOVED_INDEX) {
                    for (auto gi : inst.get_col(gj)) {
                        if (_is_global_row_active(gi)) {
                            global_idxs.emplace_back(gj);
                            break;
                        }
                    }
                }
            }

            replace_columns(global_idxs);
        } else {
            cols.clear();
            rows.clear();
        }

        return rows_left;
    }

    void replace_columns(const std::vector<idx_t>& glob_cols_idxs) {
        assert(!glob_cols_idxs.empty());

        rows.resize(local_to_global_row_idxs.size());
        for (auto& row : rows) { row.clear(); }

        idx_t ncols = glob_cols_idxs.size();
        cols.resize(ncols);
        local_to_global_col_idxs.resize(ncols);

        idx_t lj = 0;
        for (idx_t gj : glob_cols_idxs) {
            assert(std::find(fixed_cols_global_idxs.begin(), fixed_cols_global_idxs.end(), gj) == fixed_cols_global_idxs.end());

            const auto& col = inst.get_col(gj);

            cols[lj].clear();
            cols[lj].set_cost(col.get_cost());
            for (const auto gi : col) {
                if (_is_global_row_active(gi)) {
                    const auto li = global_to_local_row_idxs[gi];
                    cols[lj].emplace_back(li);
                    rows[li].emplace_back(lj);
                }
            }
            local_to_global_col_idxs[lj] = gj;
            assert(!cols[lj].empty());
            ++lj;
        }

        assert(!is_corrupted());
    }

    inline auto get_global_LB(const std::vector<real_t>& u_k) {
        real_t global_LB = std::reduce(u_k.begin(), u_k.end(), static_cast<real_t>(0.0));

        // price all active columns and add their contribution to the LB
        for (idx_t gj : inst.get_active_cols()) {
            const auto& col = inst.get_col(gj);
            real_t c_u = col.get_cost();

            bool is_empty = true;
            for (idx_t gi : col) {
                if (_is_global_row_active(gi)) {
                    is_empty = false;
                    idx_t li = global_to_local_row_idxs[gi];  // retrieve the mapped row index
                    c_u -= u_k[li];
                }
            }

            if (!is_empty && c_u < 0.0) { global_LB += c_u; }
        }

        return global_LB;
    }

private:
    struct PricedCol {
        PricedCol(idx_t idx, real_t cost) : gj(idx), c_u(cost), selected(false) { }
        idx_t gj;
        real_t c_u;
        bool selected;
    };

    std::vector<PricedCol> _cost_active_cols() {
        std::vector<PricedCol> priced_cols;
        priced_cols.reserve(rows.size());

        for (idx_t gj : inst.get_active_cols()) {
            assert(gj < inst.get_ncols());
            const auto& col = inst.get_col(gj);

            for (idx_t gi : col) {
                if (_is_global_row_active(gi)) {
                    priced_cols.emplace_back(gj, col.get_cost());
                    break;
                }
            }
        }

        return priced_cols;
    }

    real_t _price_active_cols(const std::vector<real_t>& u_k, std::vector<PricedCol>& priced_cols) {

        priced_cols.clear();

        // price all active columns and add their contribution to the LB
        real_t global_LB = std::reduce(u_k.begin(), u_k.end(), static_cast<real_t>(0.0));
        for (idx_t gj : inst.get_active_cols()) {
            const auto& col = inst.get_col(gj);
            real_t c_u = col.get_cost();

            bool is_empty = true;
            for (idx_t gi : col) {
                if (_is_global_row_active(gi)) {
                    is_empty = false;
                    idx_t li = global_to_local_row_idxs[gi];  // retrieve the mapped row index
                    c_u -= u_k[li];
                }
            }

            if (!is_empty) {  // check for empty columns
                if (c_u < 0.0) { global_LB += c_u; }
                priced_cols.emplace_back(gj, c_u);

                assert(gj < inst.get_ncols());
            }
        }

        return global_LB;
    }

    std::vector<idx_t> _select_C1_cols(std::vector<PricedCol>& priced_cols, MStar& _covering_times) {

        idx_t fivem = std::min<idx_t>(MIN_COV * inst.get_active_rows().size(), priced_cols.size());
        auto global_column_idxs = std::vector<idx_t>();
        global_column_idxs.reserve(fivem);

        std::nth_element(priced_cols.begin(), priced_cols.begin() + fivem, priced_cols.end(), [](const auto& c1, const auto& c2) { return c1.c_u < c2.c_u; });

        _covering_times.reset_uncovered(inst.get_nrows());

        for (idx_t n = 0; n < fivem; n++) {
            assert(n < priced_cols.size());
            auto& prcol = priced_cols[n];

            if (prcol.c_u >= 0.1) { continue; }

            global_column_idxs.emplace_back(prcol.gj);
            prcol.selected = true;

            for (auto gi : inst.get_col(prcol.gj)) {  // update covering info
                // here we don't care whether `gi` is removed or not
                _covering_times.cover(gi);
            }
        }

        assert(std::unique(global_column_idxs.begin(), global_column_idxs.end()) == global_column_idxs.end());

        return global_column_idxs;
    }

    std::vector<idx_t> _select_C2_cols(std::vector<PricedCol>& priced_cols, MStar& _covering_times) {
        auto global_column_idxs = std::vector<idx_t>();

        // check for still-uncovered rows
        auto still_uncovered_rows = std::vector<idx_t>();
        for (idx_t gi : inst.get_active_rows()) {
            if (_is_global_row_active(gi)) {
                _covering_times[gi] = MIN_COV - std::min<idx_t>(MIN_COV, _covering_times[gi]);
                if (_covering_times[gi] > 0) { still_uncovered_rows.emplace_back(gi); }
            }
        }

        if (!still_uncovered_rows.empty()) {

            // iterate over the remaining columns searching for the best covering
            auto best_cols = std::vector<TrivialHeap<idx_t, MIN_COV>>(inst.get_nrows());

            const auto comp_lambda = [&](const auto n1, const auto n2) {  // keep largest at the end
                assert(n1 < priced_cols.size() && n2 < priced_cols.size());
                return priced_cols[n1].c_u > priced_cols[n2].c_u;
            };

            for (idx_t n = global_column_idxs.size(); n < priced_cols.size(); ++n) {

                const auto& prcol = priced_cols[n];
                if (prcol.selected) { continue; }

                const auto& col = inst.get_col(prcol.gj);
                for (idx_t gi : col) {
                    if (!_is_global_row_active(gi) || _covering_times[gi] <= 0) { continue; }
                    assert(gi < best_cols.size());

                    if (best_cols[gi].size() < _covering_times[gi]) {
                        best_cols[gi].insert(n, comp_lambda);
                    } else {
                        if (comp_lambda(best_cols[gi].back(), n)) {
                            best_cols[gi].pop_back();
                            best_cols[gi].insert(n, comp_lambda);
                        }
                    }
                }
            }

            for (idx_t gi : still_uncovered_rows) {
                for (idx_t n : best_cols[gi]) {
                    if (!priced_cols[n].selected) {
                        priced_cols[n].selected = true;
                        global_column_idxs.emplace_back(priced_cols[n].gj);
                    }
                }
            }
        }

        assert(std::unique(global_column_idxs.begin(), global_column_idxs.end()) == global_column_idxs.end());

        return global_column_idxs;
    }


    inline bool _is_global_row_active(const idx_t gbl_idx) {
        assert(gbl_idx < global_to_local_row_idxs.size());
        return global_to_local_row_idxs[gbl_idx] != REMOVED_INDEX;
    }

    inline bool _is_local_row_active(const idx_t lcl_idx) {
        assert(lcl_idx < local_to_global_row_idxs.size());
        return local_to_global_row_idxs[lcl_idx] != REMOVED_INDEX;
    }

    Instance& inst;
    MStar covering_times;
    std::vector<Column> cols;
    std::vector<Row> rows;
    std::vector<idx_t> local_to_global_col_idxs;  // map local to original indexes
    std::vector<idx_t> fixed_cols_global_idxs;    // original indexes of locally fixed cols

    std::vector<idx_t> global_to_local_row_idxs;
    std::vector<idx_t> local_to_global_row_idxs;

    idx_t tot_removed_rows = 0;
};


#endif  // AC_CFT_INCLUDE_INSTANCE_HPP_
