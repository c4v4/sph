#ifndef AC_CFT_INCLUDE_INSTANCE_HPP_
#define AC_CFT_INCLUDE_INSTANCE_HPP_

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cassert>
#include <numeric>

#include "IndexList.hpp"
#include "MStar.hpp"
#include "TrivialHeap.hpp"
#include "cft.hpp"
#include "queue"

#define MIN_COV 10U

class Instance {
public:
    explicit Instance(const idx_t nrows_) : nrows(nrows_), active_rows(nrows, true), nactive_rows(nrows) { }

    [[nodiscard]] inline auto get_ncols() const { return cols.size(); }
    [[nodiscard]] inline idx_t get_nrows() const { return nrows; }
    [[nodiscard]] inline idx_t get_active_rows_size() const { return nactive_rows; }

    [[nodiscard]] inline auto &get_active_cols() { return active_cols; }
    [[nodiscard]] inline auto &get_fixed_cols() { return fixed_cols; }
    [[nodiscard]] inline auto &get_cols() { return cols; }
    [[nodiscard]] inline auto &get_col(idx_t idx) { return cols[idx]; }
    [[nodiscard]] inline const auto &get_col(idx_t idx) const { return cols[idx]; }

    void reset_fixing() {

        assert(std::is_sorted(active_cols.begin(), active_cols.end()));
        assert(std::is_sorted(fixed_cols.begin(), fixed_cols.end()));

        active_rows.assign(nrows, true);

        // merge active and fixed columns
        active_cols.resize(cols.size());
        std::iota(active_cols.begin(), active_cols.end(), 0);
        fixed_cols.clear();
    }
    void fix_columns(const std::vector<idx_t> &idxs) {
        for (idx_t j : idxs) {
            for (idx_t i : cols[j]) { active_rows[i] = false; }
        }

        _fix_columns(idxs);
    }
    void fix_columns(const std::vector<idx_t> &idxs, const MStar &M_star) {

        assert(fixed_cols.empty());
        assert(active_rows.size() == nrows);
        assert(M_star.size() == nrows);

        for (idx_t i = 0; i < nrows; ++i) { active_rows[i] = !M_star[i]; }
        nactive_rows = M_star.get_uncovered();

        _fix_columns(idxs);
    }

    auto compute_fixed_cost() {
        real_t fixed_cost = 0.0;
        for (idx_t j : fixed_cols) { fixed_cost += cols[j].get_cost(); }
        return fixed_cost;
    }

    inline bool is_row_active(idx_t gi) {
        assert(gi < nrows);
        return active_rows[gi];
    }

    std::vector<idx_t> add_columns(const std::vector<Column> &new_cols) {
        idx_t old_ncols = cols.size();
        idx_t ncols = new_cols.size();

        std::vector<idx_t> inserted_cols_idxs;
        inserted_cols_idxs.reserve(new_cols.size());

        for (idx_t j = 0; j < ncols; ++j) {
            cols.emplace_back(new_cols[j]);
            active_cols.emplace_back(old_ncols + j);
            inserted_cols_idxs.emplace_back(old_ncols + j);
        }

        return inserted_cols_idxs;
    }

    std::vector<idx_t> add_columns(const std::vector<real_t> &costs, const std::vector<real_t> &sol_costs, const std::vector<idx_t> &matbeg,
                                   const std::vector<idx_t> &matval) {
        idx_t old_ncols = cols.size();
        idx_t ncols = costs.size();

        std::vector<idx_t> inserted_cols_idxs;
        inserted_cols_idxs.reserve(costs.size());

        for (idx_t j = 0; j < ncols; ++j) {
            cols.emplace_back(costs[j], sol_costs[j], matval.begin() + matbeg[j], matval.begin() + matbeg[j + 1]);
            active_cols.emplace_back(old_ncols + j);
            inserted_cols_idxs.emplace_back(old_ncols + j);
        }

        return inserted_cols_idxs;
    }

private:
    void _fix_columns(const std::vector<idx_t> &idxs) {

        idx_t iok = 0;
        for (idx_t j = 0; j < cols.size(); ++j) {
            auto &col = cols[j];
            for (idx_t i : col) {
                if (active_rows[i]) {
                    active_cols[iok++] = j;
                    break;
                }
            }
        }

        active_cols.resize(iok);
        fixed_cols = idxs;
    }

    const idx_t nrows;
    std::vector<Column> cols;
    std::vector<idx_t> active_cols;
    std::vector<idx_t> fixed_cols;

    std::vector<bool> active_rows;
    idx_t nactive_rows;
};

class SubInstance {

public:
    explicit SubInstance(Instance &inst_) : inst(inst_), best_cols(inst_.get_nrows()), priced_cols(inst_.get_cols()) { }

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

    [[nodiscard]] inline auto &get_instance() { return inst; }

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

    std::vector<idx_t> add_columns(const std::vector<Column> &new_cols) {

        std::vector<idx_t> inserted_cols_idxs = inst.add_columns(new_cols);
        assert(!inserted_cols_idxs.empty());

        idx_t new_ncols = cols.size() + inserted_cols_idxs.size();
        idx_t old_ncols = cols.size();

        cols.resize(new_ncols);
        local_to_global_col_idxs.resize(new_ncols);

        idx_t lj = old_ncols;
        for (idx_t &gj : inserted_cols_idxs) {

            local_to_global_col_idxs[lj] = gj;
            cols[lj] = inst.get_col(gj);
            for (idx_t gi : cols[lj]) {
                assert(_is_global_row_active(gi));
                idx_t li = global_to_local_row_idxs[gi];
                rows[li].emplace_back(lj);
            }

            gj = lj;  // convert global to local indices
            ++lj;
        }
        assert(lj == new_ncols);
        assert(!is_corrupted());

        return inserted_cols_idxs;
    }

    void add_cols_if_changed(const std::vector<idx_t> &local_sol) {

        real_t sol_cost = std::accumulate(local_sol.begin(), local_sol.end(), 0.0, [&](real_t sum, idx_t lj) { return sum + cols[lj].get_cost(); });
        sol_cost += compute_fixed_cost();

        // fmt::print("solcost {}\n", sol_cost);

        std::vector<Column> cols_to_insert;
        cols_to_insert.reserve(local_sol.size());

        std::vector<idx_t> locl_idxs;
        locl_idxs.reserve(local_sol.size());

        for (idx_t lj : local_sol) {
            if (cols[lj].get_solcost() > sol_cost) { cols[lj].set_solcost(sol_cost); }

            idx_t gj = local_to_global_col_idxs[lj];
            if (inst.get_col(gj).size() != cols[lj].size()) {
                locl_idxs.emplace_back(lj);

                Column &gcol = cols_to_insert.emplace_back(Column(cols[lj].get_cost(), cols[lj].get_solcost(), cols[lj].size()));
                std::transform(cols[lj].begin(), cols[lj].end(), gcol.begin(), [&](idx_t n) { return local_to_global_row_idxs[n]; });
            } else if (inst.get_col(gj).get_solcost() > sol_cost) {
                // fmt::print("[{}] : {} --> {}\n", gj, inst.get_col(gj).get_solcost(), sol_cost);
                inst.get_col(gj).set_solcost(sol_cost);
            }
        }

        for (idx_t gj : fixed_cols_global_idxs) {
            if (inst.get_col(gj).get_solcost() > sol_cost) { inst.get_col(gj).set_solcost(sol_cost); }
        }

        for (idx_t gj : inst.get_fixed_cols()) {
            if (inst.get_col(gj).get_solcost() > sol_cost) { inst.get_col(gj).set_solcost(sol_cost); }
        }

        std::vector<idx_t> glob_idxs = inst.add_columns(cols_to_insert);
        assert(locl_idxs.size() == glob_idxs.size());

        for (idx_t n = 0; n < locl_idxs.size(); ++n) {
            local_to_global_col_idxs[locl_idxs[n]] = glob_idxs[n];
            assert(cols[locl_idxs[n]].size() == inst.get_col(glob_idxs[n]).size());
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
                li++;
            } else {
                global_to_local_row_idxs[gi] = REMOVED_INDEX;
            }
        }
        local_to_global_row_idxs.resize(li);

        fixed_cols_global_idxs.clear();

        _init_priced_cols(priced_cols);
        covering_times.reset_uncovered(inst.get_nrows());
        local_to_global_col_idxs.clear();
        _select_C0_cols(priced_cols, covering_times, local_to_global_col_idxs);
        _select_C2_cols(priced_cols, covering_times, local_to_global_col_idxs);

        replace_columns(local_to_global_col_idxs);

        IF_VERBOSE { fmt::print("Sub-instance size = {}x{}.\n", rows.size(), cols.size()); }

        assert(!is_corrupted());
    }

    real_t price(const std::vector<real_t> &u_k) {

        real_t global_LB = _price_active_cols(u_k, priced_cols);

        local_to_global_col_idxs.clear();
        _select_C0_cols(priced_cols, covering_times, local_to_global_col_idxs);
        _select_C1_cols(priced_cols, covering_times, local_to_global_col_idxs);
        _select_C2_cols(priced_cols, covering_times, local_to_global_col_idxs);

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

            // fmt::print("lj: {}, gj: {} \n", lj, gj);

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
        cols.resize(ncols);

        for (idx_t lj = 0; lj < glob_cols_idxs.size(); ++lj) {
            idx_t gj = glob_cols_idxs[lj];

            assert(std::find(fixed_cols_global_idxs.begin(), fixed_cols_global_idxs.end(), gj) == fixed_cols_global_idxs.end());
            assert(!inst.get_col(gj).empty());

            const auto &gcol = inst.get_col(gj);
            const idx_t csize = gcol.size();

            cols[lj].clear();
            idx_t gi = gcol[0];
            for (idx_t next_idx = 1; next_idx <= csize; ++next_idx) {
                idx_t next_gi = next_idx == csize ? REMOVED_INDEX : gcol[next_idx];
                if (_is_global_row_active(gi)) {
                    idx_t li = global_to_local_row_idxs[gi];
                    cols[lj].emplace_back(li);
                    rows[li].emplace_back(lj);
                }
                gi = next_gi;
            }
            assert(!cols[lj].empty());

            cols[lj].set_cost(gcol.get_cost());
            cols[lj].set_solcost(gcol.get_solcost());
        }

        assert(!is_corrupted());
    }

    [[nodiscard]] inline auto get_global_LB(const std::vector<real_t> &u_k) {
        real_t global_LB = std::reduce(u_k.begin(), u_k.end(), static_cast<real_t>(0.0));

        // price all active columns and add their contribution to the LB
        for (idx_t gj : inst.get_active_cols()) {
            const auto &col = inst.get_col(gj);
            real_t c_u = col.get_cost();

            // bool is_empty = true;
            for (idx_t gi : col) {
                if (_is_global_row_active(gi)) {
                    // is_empty = false;
                    idx_t li = global_to_local_row_idxs[gi];  // retrieve the mapped row index
                    c_u -= u_k[li];
                }
            }

            if (/* !is_empty && */ c_u < 0.0) { global_LB += c_u; }
        }

        return global_LB;
    }

    [[nodiscard]] inline idx_t find_local_col_idx(idx_t gj) {

        idx_t active_li = [&]() {
            for (idx_t gi : inst.get_col(gj)) {
                if (_is_global_row_active(gi)) { return global_to_local_row_idxs[gi]; }
            }
            return REMOVED_INDEX;
        }();

        if (active_li == REMOVED_INDEX) { return REMOVED_INDEX; }

        for (idx_t lj : rows[active_li]) {
            if (local_to_global_col_idxs[lj] == gj) { return lj; }
        }

        return REMOVED_INDEX;
    }

    [[nodiscard]] inline auto get_localized_solution(const std::vector<idx_t> &glob_sol) {
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

private:
    class Priced_Columns : public std::vector<Column *> {
    public:
        Priced_Columns(std::vector<Column> &cols) : col0_address(cols.data()) { }

        void reset(std::vector<Column> &cols) {
            assert(!cols.empty());
            resize(cols.size());
            col0_address = cols.data();
        }

        inline idx_t get_global_index(idx_t n) const { return get_global_index((*this)[n]); }
        inline idx_t get_global_index(Column *ptr) const {
            assert(ptr >= col0_address && col0_address < col0_address + size());
            assert(ptr != std::addressof(selected_col));
            return ptr - col0_address;
        }

        inline void select(idx_t n) { return select((*this)[n]); }
        inline void select(Column *&ptr) {
            assert(ptr >= col0_address && col0_address < col0_address + size());
            ptr = std::addressof(selected_col);
        }

        inline bool is_selected(idx_t n) const { return is_selected((*this)[n]); }
        inline bool is_selected(Column *ptr) const {
            assert(ptr >= col0_address && col0_address < col0_address + size());
            return ptr == std::addressof(selected_col);
        }

    private:
        Column *col0_address;
        Column selected_col = Column(std::numeric_limits<real_t>::max(), std::numeric_limits<real_t>::max());
    };

    struct Col_Comp {
        bool operator()(std::pair<idx_t, real_t> &p1, std::pair<idx_t, real_t> &p2) { return p1.second > p2.second; }
    };  // keep largest at the end

    void _init_priced_cols(Priced_Columns &_priced_cols) {
        _priced_cols.reset(inst.get_cols());

        idx_t p_idx = 0;
        for (idx_t gj : inst.get_active_cols()) {
            assert(gj < inst.get_ncols());

            Column *col = std::addressof(inst.get_col(gj));
            for (idx_t gi : *col) {
                if (_is_global_row_active(gi)) {
                    _priced_cols[p_idx++] = col;
                    break;
                }
            }
        }
        _priced_cols.resize(p_idx);
    }

    real_t _price_active_cols(const std::vector<real_t> &u_k, Priced_Columns &_priced_cols) {

        _priced_cols.reset(inst.get_cols());

        // price all active columns and add their contribution to the LB
        real_t global_LB = std::reduce(u_k.begin(), u_k.end(), static_cast<real_t>(0.0));
        // fmt::print("global LB base: {}\n", global_LB);

        idx_t p_idx = 0;
        for (idx_t gj : inst.get_active_cols()) {
            Column *col = std::addressof(inst.get_col(gj));
            real_t c_u = col->get_cost();

            bool is_empty = true;
            for (idx_t gi : *col) {
                if (_is_global_row_active(gi)) {
                    is_empty = false;
                    idx_t li = global_to_local_row_idxs[gi];  // retrieve the mapped row index
                    c_u -= u_k[li];
                }
            }

            if (!is_empty) {  // check for empty columns
                if (c_u < 0.0) {
                    // fmt::print("{}\n", gj);
                    global_LB += c_u;
                }
                _priced_cols[p_idx++] = col;

                assert(gj < inst.get_ncols());
            }
        }
        // fmt::print("global LB total: {}\n", global_LB);

        _priced_cols.resize(p_idx);

        return global_LB;
    }

    void _select_C0_cols(Priced_Columns &_priced_cols, MStar &_covering_times, std::vector<idx_t> &global_col_idxs) {

        idx_t fivem = std::min<idx_t>(MIN_COV * inst.get_active_rows_size(), _priced_cols.size());
        global_col_idxs.reserve(fivem);

        // std::sort(priced_cols.begin(), priced_cols.end(), [](PricedCol& c1, PricedCol& c2) { return c1.sol_c < c2.sol_c; });
        std::nth_element(_priced_cols.begin(), _priced_cols.begin() + fivem, _priced_cols.end(),
                         [](Column *c1, Column *c2) { return c1->get_solcost() < c2->get_solcost(); });

        _covering_times.reset_uncovered(inst.get_nrows());

        for (idx_t n = 0; n < fivem; n++) {
            assert(n < _priced_cols.size());

            if (_priced_cols.is_selected(n)) { continue; }

            idx_t gj = _priced_cols.get_global_index(n);
            global_col_idxs.emplace_back(gj);
            _covering_times.cover_rows(*_priced_cols[n]);
            _priced_cols.select(n);
        }

        IF_DEBUG {
            [[maybe_unused]] auto old_end = global_col_idxs.end();
            assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
        }
    }

    void _select_C1_cols(Priced_Columns &_priced_cols, MStar &_covering_times, std::vector<idx_t> &global_col_idxs) {

        idx_t fivem = std::min<idx_t>(MIN_COV * inst.get_active_rows_size(), _priced_cols.size());
        global_col_idxs.reserve(fivem);

        std::nth_element(_priced_cols.begin(), _priced_cols.begin() + fivem, _priced_cols.end(),
                         [](Column *c1, Column *c2) { return c1->get_cu() < c2->get_cu(); });

        _covering_times.reset_uncovered(inst.get_nrows());

        for (idx_t n = 0; n < fivem; n++) {
            assert(n < _priced_cols.size());

            if (_priced_cols.is_selected(n) || _priced_cols[n]->get_cu() >= 0.1) { continue; }

            idx_t gj = _priced_cols.get_global_index(n);
            global_col_idxs.emplace_back(gj);
            _covering_times.cover_rows(*_priced_cols[n]);
            _priced_cols.select(n);
        }

        IF_DEBUG {
            [[maybe_unused]] auto old_end = global_col_idxs.end();
            assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
        }
    }

    void _select_C2_cols(Priced_Columns &_priced_cols, MStar &_covering_times, std::vector<idx_t> &global_col_idxs) {

        idx_t fivem = std::min<idx_t>(MIN_COV * inst.get_active_rows_size(), _priced_cols.size());
        global_col_idxs.reserve(fivem);

        // check for still-uncovered rows
        auto still_uncovered_rows = std::vector<idx_t>();
        for (idx_t gi = 0; gi < inst.get_nrows(); ++gi) {
            if (_is_global_row_active(gi)) {
                _covering_times[gi] = MIN_COV - std::min<idx_t>(MIN_COV, _covering_times[gi]);
                if (_covering_times[gi] > 0) { still_uncovered_rows.emplace_back(gi); }
            }
        }
        if (still_uncovered_rows.empty()) { return; }

        // iterate over the remaining columns searching for the best covering
        for (auto &heap : best_cols) { heap.clear(); }

        for (idx_t n = global_col_idxs.size(); n < _priced_cols.size(); ++n) {
            if (_priced_cols.is_selected(n)) { continue; }

            Column *col = _priced_cols[n];
            auto pair = std::make_pair(n, _priced_cols[n]->get_cu());
            for (idx_t gi : *col) {
                if (!_is_global_row_active(gi) || _covering_times[gi] <= 0) { continue; }
                assert(gi < best_cols.size());

                if (best_cols[gi].size() < _covering_times[gi]) {
                    best_cols[gi].insert(pair);
                } else if (Col_Comp()(best_cols[gi].back(), pair)) {
                    best_cols[gi].pop_back();
                    best_cols[gi].insert(pair);
                }
            }
        }

        for (idx_t gi : still_uncovered_rows) {
            for (auto [n, c_u] : best_cols[gi]) {
                if (_priced_cols.is_selected(n)) { continue; }

                idx_t gj = _priced_cols.get_global_index(n);
                global_col_idxs.emplace_back(gj);
                _covering_times.cover_rows(*_priced_cols[n]);

                _priced_cols.select(n);
            }
        }

        IF_DEBUG {
            [[maybe_unused]] auto old_end = global_col_idxs.end();
            assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
        }
    }

    inline bool _is_global_row_active(const idx_t gbl_idx) {
        assert(gbl_idx < global_to_local_row_idxs.size());
        return global_to_local_row_idxs[gbl_idx] != REMOVED_INDEX;
    }

    Instance &inst;
    MStar covering_times;
    std::vector<Column> cols;
    std::vector<Row> rows;
    std::vector<idx_t> local_to_global_col_idxs;  // map local to original indexes
    std::vector<idx_t> fixed_cols_global_idxs;    // original indexes of locally fixed cols

    std::vector<idx_t> global_to_local_row_idxs;
    std::vector<idx_t> local_to_global_row_idxs;

    std::vector<TrivialHeap<std::pair<idx_t, real_t>, MIN_COV, Col_Comp>> best_cols;
    Priced_Columns priced_cols;
};

#endif  // AC_CFT_INCLUDE_INSTANCE_HPP_
