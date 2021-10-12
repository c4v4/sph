#ifndef SPH_INCLUDE_SUBINSTANCE_HPP_
#define SPH_INCLUDE_SUBINSTANCE_HPP_


#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cassert>
#include <numeric>

#include "CollectionOf.hpp"
#include "IndexList.hpp"
#include "Instance.hpp"
#include "MStar.hpp"
#include "Stopwatch.hpp"
#include "SubInstance.hpp"
#include "TrivialHeap.hpp"
#include "cft.hpp"

class SubInstance {

public:
    explicit SubInstance(Instance &inst_) : inst(inst_), best_cols(inst_.get_nrows()), fixed_cost(inst_.get_fixed_cost()) { }

    [[nodiscard]] inline auto get_ncols() const { return cols.size(); }
    [[nodiscard]] inline auto get_nrows() const { return rows.size(); }

    [[nodiscard]] inline auto get_global_col_idx(idx_t local_j) const { return local_to_global_col_idxs[local_j]; }
    [[nodiscard]] inline auto get_global_row_idx(idx_t local_i) const { return local_to_global_row_idxs[local_i]; }
    [[nodiscard]] inline auto get_local_row_idx(idx_t global_i) const { return global_to_local_row_idxs[global_i]; }

    [[nodiscard]] inline auto &get_cols() { return cols; }
    [[nodiscard]] inline auto &get_rows() { return rows; }

    [[nodiscard]] inline const auto &get_cols() const { return cols; }
    [[nodiscard]] inline const auto &get_rows() const { return rows; }

    [[nodiscard]] inline auto &get_col(idx_t idx) { return cols[idx]; }
    [[nodiscard]] inline auto &get_row(idx_t idx) { return rows[idx]; }

    [[nodiscard]] inline const auto &get_col(idx_t idx) const { return cols[idx]; }
    [[nodiscard]] inline const auto &get_row(idx_t idx) const { return rows[idx]; }

    [[nodiscard]] inline auto &get_fixed_cols() { return fixed_cols_global_idxs; }
    [[nodiscard]] inline auto get_fixed_cost() const { return fixed_cost; }
    [[nodiscard]] inline auto &get_instance() { return inst; }

    [[nodiscard]] Timer &get_timelimit() { return inst.get_timelimit(); }

    [[nodiscard]] bool is_corrupted() const {
        idx_t j_counter = 0;
        for (auto &col : cols) {
            if (std::addressof(col) != std::addressof(cols[j_counter])) {
                fmt::print("Subinstance cols iterator corrupted at {}: {} != {}\n", j_counter, (void *)std::addressof(col),
                           (void *)std::addressof(cols[j_counter]));
                return true;
            }
            ++j_counter;
        }

        for (idx_t j = 0; j < cols.size(); ++j) {
            if (cols[j].empty()) {
                IF_DEBUG { fmt::print("Col {} is empty.\n ", j); }
                return true;
            }
            if (j > get_ncols()) {
                IF_DEBUG { fmt::print("Col {} does not exist. \n Col: ", j, fmt::join(cols[j], ", ")); }
                return true;
            }

            for (idx_t i : cols[j]) {
                if (std::find(rows[i].begin(), rows[i].end(), j) == rows[i].end()) {
                    IF_DEBUG { fmt::print("Col {} not found in row {}. \n Row: ", j, i, fmt::join(rows[i], ", ")); }
                    return true;
                }
            }
        }

        for (idx_t i = 0; i < rows.size(); ++i) {
            if (rows[i].empty()) {
                IF_DEBUG { fmt::print("Row {} is empty.\n ", i); }
                return true;
            }
            if (i > get_nrows()) {
                IF_DEBUG { fmt::print("Row {} does not exist. \n Row: ", i, fmt::join(rows[i], ", ")); }
                return true;
            }

            for (idx_t j : rows[i]) {
                if (std::find(cols[j].begin(), cols[j].end(), i) == cols[j].end()) {
                    IF_DEBUG { fmt::print("Row {} not found in col {}. \n Col: ", i, j, fmt::join(cols[j], ", ")); }
                    return true;
                }
            }
        }

        return false;
    }

    [[nodiscard]] real_t get_global_LB(const std::vector<real_t> &u_k) {
        real_t global_LB = std::reduce(u_k.begin(), u_k.end(), static_cast<real_t>(0.0));

        // price all active columns and add their contribution to the LB
        for (idx_t gj : inst.get_active_cols()) {
            const auto &col = inst.get_col(gj);
            real_t c_u = col.get_cost();

            for (idx_t gi : col) {
                if (_is_global_row_active(gi)) {
                    idx_t li = global_to_local_row_idxs[gi];  // retrieve the mapped row index
                    c_u -= u_k[li];
                }
            }

            if (c_u < 0.0) { global_LB += c_u; }
        }

        return global_LB;
    }
    [[nodiscard]] idx_t find_local_col_idx(idx_t gj) {

        for (idx_t gi : inst.get_col(gj)) {
            if (_is_global_row_active(gi)) {
                idx_t active_li = global_to_local_row_idxs[gi];
                for (idx_t lj : rows[active_li]) {
                    if (local_to_global_col_idxs[lj] == gj) { return lj; }
                }
                break;
            }
        }
        return REMOVED_INDEX;
    }
    [[nodiscard]] std::vector<idx_t> get_localized_solution(const std::vector<idx_t> &glob_sol) {
        assert(glob_sol.size() >= fixed_cols_global_idxs.size());

        std::vector<idx_t> local_sol;
        local_sol.reserve(glob_sol.size() - fixed_cols_global_idxs.size() - inst.get_fixed_cols().size());
        for (idx_t gj : glob_sol) {
            idx_t lj = find_local_col_idx(gj);
            if (lj != REMOVED_INDEX) { local_sol.emplace_back(lj); }
        }

        IF_DEBUG {
            auto &sifc = fixed_cols_global_idxs;
            auto &ifc = inst.get_fixed_cols();
            for ([[maybe_unused]] idx_t lj : local_sol) { assert(std::find(glob_sol.begin(), glob_sol.end(), local_to_global_col_idxs[lj]) != glob_sol.end()); }
            for ([[maybe_unused]] idx_t gj : sifc) { assert(std::find(glob_sol.begin(), glob_sol.end(), gj) != glob_sol.end()); }
            for ([[maybe_unused]] idx_t gj : ifc) { assert(std::find(glob_sol.begin(), glob_sol.end(), gj) != glob_sol.end()); }
            for ([[maybe_unused]] idx_t gj : glob_sol) {

                [[maybe_unused]] bool check1 = [&]() {
                    for (idx_t lj : local_sol)
                        if (local_to_global_col_idxs[lj] == gj) return true;
                    return false;
                }();
                [[maybe_unused]] bool check2 = std::find(sifc.begin(), sifc.end(), gj) != sifc.end();
                [[maybe_unused]] bool check3 = std::find(ifc.begin(), ifc.end(), gj) != ifc.end();

                assert(check1 || check2 || check3);
            }
            assert(sifc.size() + ifc.size() + local_sol.size() == glob_sol.size());
        }

        return local_sol;
    }

    std::vector<idx_t> add_columns(const Cols &new_cols) {

        std::vector<idx_t> inserted_cols_idxs = inst.add_columns(new_cols);
        assert(!inserted_cols_idxs.empty());

        idx_t new_ncols = cols.size() + inserted_cols_idxs.size();
        idx_t old_ncols = cols.size();

        cols.reserve(new_ncols);
        local_to_global_col_idxs.resize(new_ncols);

        idx_t lj = old_ncols;
        for (idx_t &gj : inserted_cols_idxs) {

            local_to_global_col_idxs[lj] = gj;

            InstCol &gcol = inst.get_col(gj);
            cols.new_col_create(gcol.get_cost());
            for (idx_t gi : gcol) {
                assert(_is_global_row_active(gi));
                idx_t li = global_to_local_row_idxs[gi];
                cols.new_col_push_back(li);
                rows[li].emplace_back(lj);
            }

            gj = lj;  // convert global to local indices
            ++lj;
        }
        assert(lj == new_ncols);
        assert(!is_corrupted());

        return inserted_cols_idxs;
    }

    inline void update_sol_cost(const std::vector<idx_t> &local_sol) {
        real_t sol_cost = std::accumulate(local_sol.begin(), local_sol.end(), 0.0, [&](real_t sum, idx_t lj) { return sum + cols[lj].get_cost(); });
        sol_cost += get_fixed_cost();
        update_sol_costs(local_sol, sol_cost);
    }

    void update_sol_costs(const std::vector<idx_t> &local_sol, real_t sol_cost) {

        for (idx_t lj : local_sol) {
            InstCol &gcol = inst.get_col(local_to_global_col_idxs[lj]);
            if (gcol.get_solcost() > sol_cost) { gcol.set_solcost(sol_cost); }
        }

        for (idx_t gj : fixed_cols_global_idxs) {
            if (inst.get_col(gj).get_solcost() > sol_cost) { inst.get_col(gj).set_solcost(sol_cost); }
        }

        for (idx_t gj : inst.get_fixed_cols()) {
            if (inst.get_col(gj).get_solcost() > sol_cost) { inst.get_col(gj).set_solcost(sol_cost); }
        }
    }

    void reset() {
        local_to_global_row_idxs.resize(inst.get_nrows());
        global_to_local_row_idxs.resize(inst.get_nrows());

        idx_t li = 0;
        for (idx_t gi = 0; gi < inst.get_nrows(); ++gi) {
            if (inst.is_row_active(gi)) {
                local_to_global_row_idxs[li] = gi;
                global_to_local_row_idxs[gi] = li;
                ++li;
            } else {
                global_to_local_row_idxs[gi] = REMOVED_INDEX;
            }
        }
        local_to_global_row_idxs.resize(li);

        fixed_cols_global_idxs.clear();
        fixed_cost = inst.get_fixed_cost();

        _init_priced_cols(priced_cols);
        local_to_global_col_idxs.clear();
        covering_times.reset_uncovered(inst.get_nrows());
        std::sort(priced_cols.begin(), priced_cols.end(), [](const Priced_Col &c1, const Priced_Col &c2) { return c1.c_u < c2.c_u; });

        _select_C2_cols(priced_cols, covering_times, local_to_global_col_idxs);
        _select_C0_cols(priced_cols, local_to_global_col_idxs);

        replace_columns(local_to_global_col_idxs);

        IF_VERBOSE { fmt::print("Sub-instance size = {}x{}.\n", rows.size(), cols.size()); }

        assert(!is_corrupted());
    }
    NO_INLINE real_t price(const std::vector<real_t> &u_k) {

        real_t global_LB = _price_active_cols(u_k, priced_cols);

        local_to_global_col_idxs.clear();
        covering_times.reset_uncovered(inst.get_nrows());  // reset convered_rows to consider only reduced costs covering for C2
        std::sort(priced_cols.begin(), priced_cols.end(), [](const Priced_Col &c1, const Priced_Col &c2) { return c1.c_u < c2.c_u; });

        _select_C1_cols(priced_cols, covering_times, local_to_global_col_idxs);
        _select_C2_cols(priced_cols, covering_times, local_to_global_col_idxs);
        _select_C0_cols(priced_cols, local_to_global_col_idxs);

        replace_columns(local_to_global_col_idxs);

        IF_DEBUG {
            covering_times.reset_covered(cols, rows.size());
            assert(covering_times.get_uncovered() == 0);
        }

        return global_LB;
    }

    idx_t fix_columns(std::vector<idx_t> &local_idxs_to_fix, std::vector<real_t> &u_star) {

        if (local_idxs_to_fix.empty()) { return rows.size(); }

        // mark rows to remove
        for (idx_t lj : local_idxs_to_fix) {
            idx_t gj = local_to_global_col_idxs[lj];
            local_to_global_col_idxs[lj] = REMOVED_INDEX;

            assert(gj != REMOVED_INDEX);
            assert(std::find(fixed_cols_global_idxs.begin(), fixed_cols_global_idxs.end(), gj) == fixed_cols_global_idxs.end());

            fixed_cols_global_idxs.emplace_back(gj);
            const auto &col = cols[lj];
            for (idx_t li : col) {
                if (local_to_global_row_idxs[li] == REMOVED_INDEX) { continue; }

                idx_t gi = local_to_global_row_idxs[li];
                local_to_global_row_idxs[li] = REMOVED_INDEX;
                global_to_local_row_idxs[gi] = REMOVED_INDEX;
            }
        }

        fixed_cost = inst.get_fixed_cost();
        for (idx_t gj : fixed_cols_global_idxs) { fixed_cost += inst.get_col(gj).get_cost(); }

        // compact rows
        idx_t li = 0;
        while (local_to_global_row_idxs[li] != REMOVED_INDEX) { ++li; }

        idx_t rows_left = li;
        for (; li < rows.size(); ++li) {
            if (local_to_global_row_idxs[li] == REMOVED_INDEX) { continue; }

            idx_t gi = local_to_global_row_idxs[li];

            assert(_is_global_row_active(gi));
            assert(!rows[rows_left].empty());

            global_to_local_row_idxs[gi] = rows_left;
            local_to_global_row_idxs[rows_left] = gi;
            u_star[rows_left] = u_star[li];
            ++rows_left;
        }
        local_to_global_row_idxs.resize(rows_left);
        u_star.resize(rows_left);

        IF_DEBUG {
            for ([[maybe_unused]] idx_t gi : local_to_global_row_idxs) { assert(_is_global_row_active(gi)); }
        }

        if (rows_left == 0) {
            cols.clear();
            rows.clear();
            return 0;
        }

        idx_t lj = 0;
        for (idx_t gj : local_to_global_col_idxs) {
            if (gj == REMOVED_INDEX) { continue; }

            for (auto gi : inst.get_col(gj)) {
                if (_is_global_row_active(gi)) {
                    local_to_global_col_idxs[lj++] = gj;
                    break;
                }
            }
        }
        local_to_global_col_idxs.resize(lj);

        replace_columns(local_to_global_col_idxs);
        return rows_left;
    }

    void replace_columns(const std::vector<idx_t> &glob_cols_idxs) {
        assert(!glob_cols_idxs.empty());

        rows.resize(local_to_global_row_idxs.size());
        for (auto &row : rows) { row.clear(); }

        idx_t ncols = glob_cols_idxs.size();

        cols.clear();
        cols.reserve(ncols);

        for (idx_t lj = 0; lj < glob_cols_idxs.size(); ++lj) {
            idx_t gj = glob_cols_idxs[lj];

            assert(std::find(fixed_cols_global_idxs.begin(), fixed_cols_global_idxs.end(), gj) == fixed_cols_global_idxs.end());
            assert(!inst.get_col(gj).empty());

            const auto &gcol = inst.get_col(gj);

            cols.new_col_create(gcol.get_cost());
            for (idx_t gi : gcol) {
                if (_is_global_row_active(gi)) {
                    idx_t li = global_to_local_row_idxs[gi];
                    cols.new_col_push_back(li);
                    rows[li].emplace_back(lj);
                }
            }
            assert(!cols[lj].empty());
        }
        assert(!is_corrupted());
    }


private:
    struct Priced_Col {
        idx_t j;
        real_t c_u;
        real_t sol_cost;
    };

    class Priced_Columns : public std::vector<Priced_Col> {
    public:
        Priced_Columns() { }

        void reset(idx_t ncols) {
            assert(ncols > 0);
            resize(ncols);
        }

        inline void select(idx_t n) {
            (*this)[n].j = REMOVED_INDEX;
            (*this)[n].c_u = (*this)[n].sol_cost = REAL_MAX;
        }
        inline bool is_selected(idx_t n) const { return (*this)[n].j == REMOVED_INDEX; }
    };

    struct Col_Comp {
        bool operator()(std::pair<idx_t, real_t> &p1, std::pair<idx_t, real_t> &p2) { return p1.second > p2.second; }
    };  // keep largest at the end

    void _init_priced_cols(Priced_Columns &_priced_cols) {
        _priced_cols.reset(inst.get_active_cols().size());

        idx_t p_idx = 0;
        for (idx_t gj : inst.get_active_cols()) {
            assert(gj < inst.get_ncols());

            InstCol &col = inst.get_col(gj);
            for (idx_t gi : col) {
                if (_is_global_row_active(gi)) {
                    _priced_cols[p_idx++] = {gj, col.get_cost(), col.get_solcost()};
                    break;
                }
            }
        }
        _priced_cols.resize(p_idx);
    }

    real_t _price_active_cols(const std::vector<real_t> &u_k, Priced_Columns &_priced_cols) {

        _priced_cols.reset(inst.get_active_cols().size());

        // price all active columns and add their contribution to the LB
        real_t global_LB = std::reduce(u_k.begin(), u_k.end(), static_cast<real_t>(0.0));

        idx_t p_idx = 0;
        for (idx_t gj : inst.get_active_cols()) {
            assert(gj < inst.get_ncols());

            InstCol &col = inst.get_col(gj);
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

                _priced_cols[p_idx++] = {gj, c_u, col.get_solcost()};
            }
        }

        _priced_cols.resize(p_idx);

        return global_LB;
    }

    NO_INLINE void _select_C0_cols(Priced_Columns &_priced_cols, std::vector<idx_t> &global_col_idxs) {
        idx_t fivem = std::min<idx_t>(HARD_CAP, std::min<idx_t>(MIN_SOLCOST_COV * inst.get_active_rows_size(), _priced_cols.size()));
        global_col_idxs.reserve(fivem);

        std::nth_element(_priced_cols.begin(), _priced_cols.begin() + fivem, _priced_cols.end(),
                         [](const Priced_Col &c1, const Priced_Col &c2) { return c1.sol_cost < c2.sol_cost; });

        if (_priced_cols[0].sol_cost == REAL_MAX) { return; }

        for (idx_t n = 0; n < fivem; ++n) {
            assert(n < _priced_cols.size());

            if (_priced_cols.is_selected(n) || _priced_cols[n].sol_cost == REAL_MAX) { continue; }

            idx_t gj = _priced_cols[n].j;
            assert(gj < inst.get_ncols());

            global_col_idxs.emplace_back(gj);
            _priced_cols.select(n);
        }

        IF_DEBUG {
            [[maybe_unused]] auto old_end = global_col_idxs.end();
            assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
        }
    }
    NO_INLINE void _select_C1_cols(Priced_Columns &_priced_cols, MStar &_covering_times, std::vector<idx_t> &global_col_idxs) {

        idx_t fivem = std::min<idx_t>(HARD_CAP, std::min<idx_t>(MIN_COV * inst.get_active_rows_size(), _priced_cols.size()));
        global_col_idxs.reserve(fivem);

        for (idx_t n = 0; n < fivem; n++) {
            assert(n < _priced_cols.size());

            if (_priced_cols.is_selected(n) || _priced_cols[n].c_u >= 0.1) { continue; }

            idx_t gj = _priced_cols[n].j;
            assert(gj < inst.get_ncols());

            global_col_idxs.emplace_back(gj);
            _covering_times.cover_rows(inst.get_col(gj));
            _priced_cols.select(n);
        }

        IF_DEBUG {
            [[maybe_unused]] auto old_end = global_col_idxs.end();
            assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
        }
    }
    NO_INLINE void _select_C2_cols(Priced_Columns &_priced_cols, MStar &_covering_times, std::vector<idx_t> &global_col_idxs) {

        assert(std::is_sorted(_priced_cols.begin() + global_col_idxs.size(), _priced_cols.end(),
                              [](const Priced_Col &c1, const Priced_Col &c2) { return c1.c_u < c2.c_u; }));

        idx_t nrows = get_nrows() == 0 ? inst.get_nrows() : get_nrows();
        idx_t min_cov = std::min<idx_t>(MIN_COV, HARD_CAP / nrows);
        idx_t fivem = std::min<idx_t>(min_cov * inst.get_active_rows_size(), _priced_cols.size());
        global_col_idxs.reserve(fivem);

        // check for still-uncovered rows
        idx_t rows_to_cover = 0;
        for (idx_t gi = 0; gi < inst.get_nrows(); ++gi) {
            if (_is_global_row_active(gi)) {
                _covering_times[gi] = min_cov - std::min<idx_t>(min_cov, _covering_times[gi]);
                rows_to_cover += static_cast<idx_t>(_covering_times[gi] > 0);
            } else {
                _covering_times[gi] = 0;
            }
        }

        for (auto &heap : best_cols) { heap.clear(); }

        for (idx_t n = global_col_idxs.size(); n < _priced_cols.size(); ++n) {
            assert(!_priced_cols.is_selected(n));

            InstCol &col = inst.get_col(_priced_cols[n].j);
            auto pair = std::make_pair(n, _priced_cols[n].c_u);
            for (idx_t gi : col) {

                if (_covering_times[gi] == 0) { continue; }

                assert(gi < best_cols.size());
                --_covering_times[gi];
                best_cols[gi].push_back(pair);

                rows_to_cover -= static_cast<idx_t>(_covering_times[gi] == 0);
                if (rows_to_cover == 0) {
                    assert(std::count(_covering_times.begin(), _covering_times.end(), 0) == _covering_times.size());
                    break;
                }
            }
        }

        for (auto &heap : best_cols) {
            for (auto [n, c_u] : heap) {
                if (_priced_cols.is_selected(n)) { continue; }

                idx_t gj = _priced_cols[n].j;
                assert(gj < inst.get_ncols());

                global_col_idxs.emplace_back(gj);

                _priced_cols.select(n);
            }
        }

        IF_DEBUG {
            [[maybe_unused]] auto old_end = global_col_idxs.end();
            assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
        }
    }


    [[nodiscard]] inline bool _is_global_row_active(const idx_t gbl_idx) {
        assert(gbl_idx < global_to_local_row_idxs.size());
        return global_to_local_row_idxs[gbl_idx] != REMOVED_INDEX;
    }

private:
    Instance &inst;
    MStar covering_times;
    SubInstCols cols;
    std::vector<Row> rows;
    std::vector<idx_t> local_to_global_col_idxs;  // map local to original indexes
    std::vector<idx_t> fixed_cols_global_idxs;    // original indexes of locally fixed cols

    std::vector<idx_t> global_to_local_row_idxs;
    std::vector<idx_t> local_to_global_row_idxs;

    std::vector<TrivialHeap<std::pair<idx_t, real_t>, MIN_COV, Col_Comp>> best_cols;
    Priced_Columns priced_cols;

    real_t fixed_cost;
};

#endif  // SPH_INCLUDE_SUBINSTANCE_HPP_