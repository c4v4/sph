#ifndef SPH_INCLUDE_SUBINSTANCE_HPP_
#define SPH_INCLUDE_SUBINSTANCE_HPP_


#include <algorithm>
#include <cassert>
#include <numeric>

#include "CollectionOf.hpp"
#include "IndexList.hpp"
#include "Instance.hpp"
#include "MStar.hpp"
#include "Stopwatch.hpp"
#include "cft.hpp"
#include "fmt/core.h"
#include "fmt/ranges.h"

namespace sph {

    class SubInstance {

    public:
        explicit SubInstance(Instance &inst_) : inst(inst_), fixed_cost(inst_.get_fixed_cost()) { }

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
                    SPH_DEBUG { fmt::print("Col {} is empty.\n ", j); }
                    return true;
                }
                if (j > get_ncols()) {
                    SPH_DEBUG { fmt::print("Col {} does not exist. \n Col: {}", j, fmt::join(cols[j], ", ")); }
                    return true;
                }

                for (idx_t i : cols[j]) {
                    if (std::find(rows[i].begin(), rows[i].end(), j) == rows[i].end()) {
                        SPH_DEBUG { fmt::print("Col {} not found in row {}. \n Row: {}", j, i, fmt::join(rows[i], ", ")); }
                        return true;
                    }
                }
            }

            for (idx_t i = 0; i < rows.size(); ++i) {
                if (rows[i].empty()) {
                    SPH_DEBUG { fmt::print("Row {} is empty.\n ", i); }
                    return true;
                }
                if (i > get_nrows()) {
                    SPH_DEBUG { fmt::print("Row {} does not exist. \n Row: {}", i, fmt::join(rows[i], ", ")); }
                    return true;
                }

                for (idx_t j : rows[i]) {
                    if (std::find(cols[j].begin(), cols[j].end(), i) == cols[j].end()) {
                        SPH_DEBUG { fmt::print("Row {} not found in col {}. \n Col: {}", i, j, fmt::join(cols[j], ", ")); }
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
            return NOT_AN_INDEX;
        }
        [[nodiscard]] std::vector<idx_t> get_localized_solution(const std::vector<idx_t> &glob_sol) {
            assert(glob_sol.size() >= fixed_cols_global_idxs.size());

            std::vector<idx_t> local_sol;
            local_sol.reserve(glob_sol.size() - fixed_cols_global_idxs.size() - inst.get_fixed_cols().size());
            for (idx_t gj : glob_sol) {
                idx_t lj = find_local_col_idx(gj);
                if (lj != NOT_AN_INDEX) { local_sol.emplace_back(lj); }
            }

            SPH_DEBUG {
                auto &sifc = fixed_cols_global_idxs;
                auto &ifc = inst.get_fixed_cols();
                for ([[maybe_unused]] idx_t lj : local_sol) {
                    assert(std::find(glob_sol.begin(), glob_sol.end(), local_to_global_col_idxs[lj]) != glob_sol.end());
                }
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

        template <typename ColContainer>
        std::vector<idx_t> add_columns(const ColContainer &new_cols) {

            std::vector<idx_t> inserted_cols_idxs = inst.add_columns(new_cols);
            assert(!inserted_cols_idxs.empty());

            idx_t new_ncols = cols.size() + inserted_cols_idxs.size();
            idx_t old_ncols = cols.size();

            cols.reserve(new_ncols);
            local_to_global_col_idxs.resize(new_ncols);

            idx_t lj = old_ncols;
            for (idx_t &gj : inserted_cols_idxs) {

                local_to_global_col_idxs[lj] = gj;

                Column &gcol = inst.get_col(gj);
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
            real_t sol_cost = std::accumulate(local_sol.begin(), local_sol.end(), 0.0,
                                              [&](real_t sum, idx_t lj) { return sum + cols[lj].get_cost(); });
            sol_cost += get_fixed_cost();
            update_sol_costs(local_sol, sol_cost);
        }

        void update_sol_costs(const std::vector<idx_t> &local_sol, real_t sol_cost) {

            for (idx_t lj : local_sol) {
                Column &gcol = inst.get_col(local_to_global_col_idxs[lj]);
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
                    global_to_local_row_idxs[gi] = NOT_AN_INDEX;
                }
            }
            local_to_global_row_idxs.resize(li);

            local_to_global_col_idxs.clear();
            fixed_cols_global_idxs.clear();
            fixed_cost = inst.get_fixed_cost();

            inst.fill_with_best_columns(local_to_global_col_idxs);
            replace_columns(local_to_global_col_idxs);

            SPH_VERBOSE(2) { fmt::print("Sub-instance size = {}x{}.\n", rows.size(), cols.size()); }

            assert(!is_corrupted());
        }

        real_t price(const std::vector<real_t> &u_k) {

            global_u_k.assign(inst.get_nrows(), 0);
            for (idx_t i = 0; i < u_k.size(); ++i) { global_u_k[local_to_global_row_idxs[i]] = u_k[i]; }

            local_to_global_col_idxs.clear();
            real_t global_LB = inst.fill_with_best_columns(local_to_global_col_idxs, global_u_k);

            replace_columns(local_to_global_col_idxs);

            return global_LB;
        }

        idx_t fix_columns(std::vector<idx_t> &local_idxs_to_fix, std::vector<real_t> &u_star) {

            if (local_idxs_to_fix.empty()) { return rows.size(); }

            // mark rows to remove
            for (idx_t lj : local_idxs_to_fix) {
                idx_t gj = local_to_global_col_idxs[lj];
                local_to_global_col_idxs[lj] = NOT_AN_INDEX;

                assert(gj != NOT_AN_INDEX);
                assert(std::find(fixed_cols_global_idxs.begin(), fixed_cols_global_idxs.end(), gj) == fixed_cols_global_idxs.end());

                fixed_cols_global_idxs.emplace_back(gj);
                const auto &col = cols[lj];
                for (idx_t li : col) {
                    if (local_to_global_row_idxs[li] == NOT_AN_INDEX) { continue; }

                    idx_t gi = local_to_global_row_idxs[li];
                    local_to_global_row_idxs[li] = NOT_AN_INDEX;
                    global_to_local_row_idxs[gi] = NOT_AN_INDEX;
                }
            }

            fixed_cost = inst.get_fixed_cost();
            for (idx_t gj : fixed_cols_global_idxs) { fixed_cost += inst.get_col(gj).get_cost(); }

            // compact rows
            idx_t li = 0;
            while (local_to_global_row_idxs[li] != NOT_AN_INDEX) { ++li; }

            idx_t rows_left = li;
            for (; li < rows.size(); ++li) {
                if (local_to_global_row_idxs[li] == NOT_AN_INDEX) { continue; }

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

            SPH_DEBUG {
                for ([[maybe_unused]] idx_t gi : local_to_global_row_idxs) { assert(_is_global_row_active(gi)); }
            }

            if (rows_left == 0) {
                cols.clear();
                rows.clear();
                return 0;
            }

            idx_t lj = 0;
            for (idx_t gj : local_to_global_col_idxs) {
                if (gj == NOT_AN_INDEX) { continue; }

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
        [[nodiscard]] inline bool _is_global_row_active(const idx_t gbl_idx) {
            assert(gbl_idx < global_to_local_row_idxs.size());
            return global_to_local_row_idxs[gbl_idx] != NOT_AN_INDEX;
        }

    private:
        Instance &inst;
        SubInstCols cols;
        std::vector<Row> rows;
        std::vector<idx_t> local_to_global_col_idxs;  // map local to original indexes
        std::vector<idx_t> fixed_cols_global_idxs;    // original indexes of locally fixed cols

        std::vector<idx_t> global_to_local_row_idxs;
        std::vector<idx_t> local_to_global_row_idxs;

        std::vector<real_t> global_u_k;
        real_t fixed_cost;
    };

}  // namespace sph

#endif  // SPH_INCLUDE_SUBINSTANCE_HPP_
