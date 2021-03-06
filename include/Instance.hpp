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

#ifndef SPH_INCLUDE_INSTANCE_HPP_
#define SPH_INCLUDE_INSTANCE_HPP_

#include <algorithm>
#include <cassert>
#include <numeric>

#include "CollectionOf.hpp"
#include "MStar.hpp"
#include "Stopwatch.hpp"
#include "UniqueColSet.hpp"
#include "cft.hpp"

namespace sph {

    constexpr unsigned SUBINST_MIN_COV = 4U;
    constexpr unsigned SUBINST_MIN_SOLCOST_COV = 4U;
    constexpr unsigned SUBINST_HARD_CAP = 15'000U;
    constexpr unsigned INST_HARD_CAP = 500'000U;

    /**
     * @brief Represents a complete instance of a Set Partitioning problem.
     *
     */
    class Instance {
    public:
        explicit Instance(const idx_t nrows_)
            : nrows(nrows_), active_rows(nrows, true), nactive_rows(nrows), fixed_cost(0.0), ncols_constr(0) { }

        [[nodiscard]] inline idx_t get_ncols() const { return cols.size(); }
        [[nodiscard]] inline idx_t get_nrows() const { return nrows; }
        [[nodiscard]] inline idx_t get_active_rows_size() const { return nactive_rows; }

        [[nodiscard]] inline std::vector<idx_t> &get_active_cols() { return active_cols; }
        [[nodiscard]] inline const std::vector<idx_t> &get_active_cols() const { return active_cols; }
        [[nodiscard]] inline std::vector<idx_t> &get_fixed_cols() { return fixed_cols; }
        [[nodiscard]] inline const std::vector<idx_t> &get_fixed_cols() const { return fixed_cols; }
        [[nodiscard]] inline UniqueColSet &get_cols() { return cols; }
        [[nodiscard]] inline const UniqueColSet &get_cols() const { return cols; }
        [[nodiscard]] inline Column &get_col(idx_t idx) { return cols[idx]; }
        [[nodiscard]] inline const Column &get_col(idx_t idx) const { return cols[idx]; }
        [[nodiscard]] inline real_t get_fixed_cost() const { return fixed_cost; }
        [[nodiscard]] inline idx_t get_ncols_constr() const {
            return ncols_constr > fixed_cols.size() ? ncols_constr - fixed_cols.size() : 0;
        }
        [[nodiscard]] inline KeepStrat get_keepcol_strategy() { return keep_col_strat; }
        inline void set_ncols_constr(idx_t ncols_constr_) { ncols_constr = ncols_constr_; }
        inline void set_max_routes(idx_t inst_max_routes_) { inst_max_routes = inst_max_routes_; }
        inline void set_keepcol_strategy(KeepStrat keep_col_strat_) { keep_col_strat = keep_col_strat_; }
        inline void set_timelimit(double seconds) { timelimit = Timer(seconds); }
        [[nodiscard]] inline Timer &get_timelimit() { return timelimit; }

        [[nodiscard]] inline bool is_row_active(idx_t gi) const {
            assert(gi < nrows);
            return active_rows[gi];
        }

        inline void set_new_best_callback(NewBestCallback cb) { new_best_callback = cb; }

        inline void new_best_cb(GlobalSolution &sol) {
            if (new_best_callback)
                new_best_callback(*this, sol);
        }

        void inline fix_columns(const std::vector<idx_t> &idxs) {
            for (idx_t j : idxs) {
                for (idx_t i : cols[j]) {
                    active_rows[i] = false;
                }
            }
            _fix_columns(idxs);
        }

        void update_sol_costs(const std::vector<idx_t> sol_cols, real_t sol_cost) {
            for (idx_t gj : sol_cols) {
                if (cols[gj].get_solcost() > sol_cost) {
                    cols[gj].set_solcost(sol_cost);
                }
            }
            for (idx_t gj : fixed_cols) {
                if (cols[gj].get_solcost() > sol_cost) {
                    cols[gj].set_solcost(sol_cost);
                }
            }
        }

        void fix_columns(const std::vector<idx_t> &idxs, const MStar &M_star) {
            assert(fixed_cols.empty());
            assert(active_rows.size() == nrows);
            assert(M_star.size() == nrows);

            for (idx_t i = 0; i < nrows; ++i) {
                active_rows[i] = !M_star[i];
            }
            nactive_rows = M_star.get_uncovered();

            _fix_columns(idxs);
        }


        void reset_fixing() {
            assert(std::is_sorted(active_cols.begin(), active_cols.end()));
            assert(std::is_sorted(fixed_cols.begin(), fixed_cols.end()));

            active_rows.assign(nrows, true);

            // merge active and fixed columns
            active_cols.resize(cols.size());
            std::iota(active_cols.begin(), active_cols.end(), 0);

            fixed_cols.clear();
            fixed_cost = 0.0;
        }

        template <typename... _Args>
        idx_t add_column(_Args &&...args) {
            auto [j, inserted_new] = cols.add_column(std::forward<_Args>(args)...);
            assert(!inserted_new || j == cols.size() - 1);
            return j;
        }

        template <typename ColContainer>
        std::vector<idx_t> add_columns(const ColContainer &new_cols) {

            std::vector<idx_t> inserted_cols_idxs;
            inserted_cols_idxs.reserve(new_cols.size());

            for (auto &new_col : new_cols) {
                idx_t inserted_idx = add_column(new_col);
                if (inserted_idx != NOT_AN_INDEX) {
                    inserted_cols_idxs.emplace_back(inserted_idx);
                }
            }

            return inserted_cols_idxs;
        }

        template <typename ColContainer>
        std::vector<idx_t> add_columns(ColContainer &&new_cols) {

            std::vector<idx_t> inserted_cols_idxs;
            inserted_cols_idxs.reserve(new_cols.size());

            for (auto &new_col : new_cols) {
                idx_t inserted_idx = add_column(std::move(new_col));
                if (inserted_idx != NOT_AN_INDEX) {
                    inserted_cols_idxs.emplace_back(inserted_idx);
                }
            }

            return inserted_cols_idxs;
        }

        template <typename CostVec, typename SolCostVec, typename MatBegVec, typename MatValVec>
        std::vector<idx_t> add_columns(const CostVec &costs, const SolCostVec &sol_costs, const MatBegVec &matbeg,
                                       const MatValVec &matval) {
            assert(costs.size() == sol_costs.size() && costs.size() <= matbeg.size());

            idx_t ncols = costs.size();

            std::vector<idx_t> inserted_cols_idxs;
            inserted_cols_idxs.reserve(costs.size());

            for (idx_t j = 0; j < ncols - 1; ++j) {
                idx_t inserted_idx = add_column(&matval[matbeg[j]], &matval[matbeg[j + 1]], costs[j], sol_costs[j]);
                inserted_cols_idxs.emplace_back(inserted_idx);
            }

            // Matbeg might prbably contains only the starts
            idx_t inserted_idx = add_column(&matval[matbeg[ncols - 1]], &matval[matval.size() - 1], costs[ncols - 1], sol_costs[ncols - 1]);
            inserted_cols_idxs.emplace_back(inserted_idx);

            assert(cols.size() <= costs.size());
            return inserted_cols_idxs;
        }

        void fill_with_best_columns(std::vector<idx_t> &global_idxs) {
            _init_priced_cols(priced_cols);
            covering_times.reset_uncovered(nrows);
            std::sort(priced_cols.begin(), priced_cols.end(), [](const Priced_Col &c1, const Priced_Col &c2) { return c1.c_u < c2.c_u; });

            _select_C2_cols(priced_cols, covering_times, global_idxs, SUBINST_MIN_COV, SUBINST_HARD_CAP);
            _select_C3_cols(priced_cols, global_idxs, SUBINST_MIN_SOLCOST_COV, SUBINST_HARD_CAP);
        }

        real_t fill_with_best_columns(std::vector<idx_t> &global_idxs, const std::vector<real_t> &u_k) {

            real_t global_LB = _price_active_cols(u_k, priced_cols);

            covering_times.reset_uncovered(nrows);  // reset convered_rows to consider only reduced costs covering for C2
            std::sort(priced_cols.begin(), priced_cols.end(), [](const Priced_Col &c1, const Priced_Col &c2) { return c1.c_u < c2.c_u; });

            _select_C1_cols(priced_cols, covering_times, global_idxs, SUBINST_MIN_COV, SUBINST_HARD_CAP);
            _select_C2_cols(priced_cols, covering_times, global_idxs, SUBINST_MIN_COV, SUBINST_HARD_CAP);
            _select_C3_cols(priced_cols, global_idxs, SUBINST_MIN_SOLCOST_COV, SUBINST_HARD_CAP);

            return global_LB;
        }

        /**
         * @brief Prune columns maintaining only the best ones.
         *
         * @tparam Hard_cap
         * @param u_k
         * @return std::vector<idx_t> map from old indexes to new ones to translate pre-existing solutions.
         *          NOT_AN_INDEX if the column has been removed.
         */
        std::vector<idx_t> prune_instance(const std::vector<real_t> &u_k) {
            if (cols.size() > 3 * inst_max_routes) {
                std::vector<idx_t> idxs_to_keep;
                idxs_to_keep.reserve(inst_max_routes);

                _price_active_cols(u_k, priced_cols);
                covering_times.reset_uncovered(nrows);  // reset convered_rows to consider only reduced costs covering for C2

                _select_C1_cols(priced_cols, covering_times, idxs_to_keep, MAX_INDEX, inst_max_routes);
                _select_C2_cols(priced_cols, covering_times, idxs_to_keep, MAX_INDEX, inst_max_routes);
                _select_C3_cols(priced_cols, idxs_to_keep, MAX_INDEX, inst_max_routes);

                UniqueColSet new_cols;
                new_cols.reserve(idxs_to_keep.size());
                std::vector<idx_t> old_to_new_idx_map(cols.size(), NOT_AN_INDEX);

                for (idx_t gj : idxs_to_keep) {
                    old_to_new_idx_map[gj] = new_cols.size();
                    new_cols.add_column(cols[gj]);
                }
                std::swap(cols, new_cols);

                return old_to_new_idx_map;
            }

            return std::vector<idx_t>();
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
                (*this)[n].j = NOT_AN_INDEX;
                (*this)[n].c_u = (*this)[n].sol_cost = REAL_MAX;
            }
            inline bool is_selected(idx_t n) const { return (*this)[n].j == NOT_AN_INDEX; }
        };

        void _init_priced_cols(Priced_Columns &_priced_cols) {
            _priced_cols.reset(active_cols.size());

            idx_t p_idx = 0;
            for (idx_t gj : active_cols) {
                assert(gj < cols.size());

                Column &col = cols[gj];
                for (idx_t gi : col) {
                    if (active_rows[gi]) {
                        _priced_cols[p_idx++] = {gj, col.get_cost(), col.get_solcost()};
                        break;
                    }
                }
            }
            _priced_cols.resize(p_idx);
        }

        real_t _price_active_cols(const std::vector<real_t> &u_k, Priced_Columns &_priced_cols) {

            _priced_cols.reset(active_cols.size());

            // price all active columns and add their contribution to the LB
            real_t global_LB = std::accumulate(u_k.begin(), u_k.end(), static_cast<real_t>(0.0));

            idx_t p_idx = 0;
            for (idx_t gj : active_cols) {
                assert(gj < cols.size());

                Column &col = cols[gj];
                real_t c_u = col.get_cost();

                bool is_empty = true;
                for (idx_t gi : col) {
                    if (active_rows[gi]) {
                        is_empty = false;
                        c_u -= u_k[gi];  // NOTE: multipliers need to be adapted to global multipliers!!!!!
                    }
                }

                if (!is_empty) {  // check for empty columns
                    if (c_u < 0.0) {
                        global_LB += c_u;
                    }

                    _priced_cols[p_idx++] = {gj, c_u, col.get_solcost()};
                }
            }

            _priced_cols.resize(p_idx);

            return global_LB;
        }

        void _select_C1_cols(Priced_Columns &_priced_cols, MStar &_covering_times, std::vector<idx_t> &global_col_idxs, idx_t Min_cov,
                             idx_t Hard_cap) {
            if (nactive_rows == 0) {
                return;
            }

            idx_t fivem = std::min<idx_t>(Hard_cap, std::min<idx_t>(Min_cov * nactive_rows, _priced_cols.size()));
            global_col_idxs.reserve(fivem);

            std::sort(_priced_cols.begin(), _priced_cols.end(), [](const Priced_Col &c1, const Priced_Col &c2) { return c1.c_u < c2.c_u; });

            for (idx_t n = 0; n < fivem; n++) {
                assert(n < _priced_cols.size());

                if (_priced_cols.is_selected(n) || _priced_cols[n].c_u >= 0.1) {
                    continue;
                }

                idx_t gj = _priced_cols[n].j;
                assert(gj < cols.size());
                assert(std::count_if(cols[gj].begin(), cols[gj].end(), [&](idx_t i) { return active_rows[i]; }) > 0);

                global_col_idxs.emplace_back(gj);
                _covering_times.cover_rows(cols[gj]);
                _priced_cols.select(n);
            }

            SPH_DEBUG {
                [[maybe_unused]] auto old_end = global_col_idxs.end();
                assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
            }
        }

        void _select_C2_cols(Priced_Columns &_priced_cols, MStar &_covering_times, std::vector<idx_t> &global_col_idxs, idx_t Min_cov,
                             idx_t Hard_cap) {

            assert(std::is_sorted(_priced_cols.begin() + global_col_idxs.size(), _priced_cols.end(),
                                  [](const Priced_Col &c1, const Priced_Col &c2) { return c1.c_u < c2.c_u; }));

            if (nactive_rows == 0) {
                return;
            }

            idx_t min_cov = std::min<idx_t>(Min_cov, Hard_cap / nactive_rows);
            idx_t fivem = std::min<idx_t>(min_cov * nactive_rows, _priced_cols.size());
            global_col_idxs.reserve(fivem);

            // check for still-uncovered rows
            idx_t rows_to_cover = 0;
            for (idx_t gi = 0; gi < nrows; ++gi) {
                if (active_rows[gi]) {
                    _covering_times[gi] = min_cov - std::min<idx_t>(min_cov, _covering_times[gi]);
                    rows_to_cover += static_cast<idx_t>(_covering_times[gi] > 0);
                } else {
                    _covering_times[gi] = 0;
                }
            }

            for (idx_t n = global_col_idxs.size(); n < _priced_cols.size(); ++n) {
                assert(!_priced_cols.is_selected(n));

                Column &col = cols[_priced_cols[n].j];
                for (idx_t gi : col) {
                    if (_covering_times[gi] == 0) {
                        continue;
                    }

                    --_covering_times[gi];

                    if (!_priced_cols.is_selected(n)) {
                        idx_t gj = _priced_cols[n].j;
                        assert(gj < cols.size());
                        assert(std::count_if(cols[gj].begin(), cols[gj].end(), [&](idx_t i) { return active_rows[i]; }) > 0);

                        global_col_idxs.emplace_back(gj);
                        _priced_cols.select(n);
                    }

                    rows_to_cover -= static_cast<idx_t>(_covering_times[gi] == 0);
                    if (rows_to_cover == 0) {
                        assert(std::count(_covering_times.begin(), _covering_times.end(), 0) == _covering_times.size());
                        break;
                    }
                }
            }

            SPH_DEBUG {
                [[maybe_unused]] auto old_end = global_col_idxs.end();
                assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
            }
        }

        void _select_C3_cols(Priced_Columns &_priced_cols, std::vector<idx_t> &global_col_idxs, idx_t Min_cov, idx_t Hard_cap) {
            if (nactive_rows == 0) {
                return;
            }

            idx_t fivem = std::min<idx_t>(Hard_cap, std::min<idx_t>(Min_cov * nactive_rows, _priced_cols.size()));
            global_col_idxs.reserve(fivem);

            std::nth_element(_priced_cols.begin(), _priced_cols.begin() + fivem, _priced_cols.end(),
                             [](const Priced_Col &c1, const Priced_Col &c2) { return c1.sol_cost < c2.sol_cost; });

            auto min_e = std::min_element(_priced_cols.begin(), _priced_cols.begin() + fivem,
                                          [](const Priced_Col &c1, const Priced_Col &c2) { return c1.sol_cost < c2.sol_cost; });

            if (min_e->sol_cost == REAL_MAX) {
                return;
            }

            for (idx_t n = 0; n < fivem; ++n) {
                assert(n < _priced_cols.size());

                if (_priced_cols.is_selected(n) || _priced_cols[n].sol_cost == REAL_MAX) {
                    continue;
                }

                idx_t gj = _priced_cols[n].j;
                assert(gj < cols.size());
                assert(std::count_if(cols[gj].begin(), cols[gj].end(), [&](idx_t i) { return active_rows[i]; }) > 0);

                global_col_idxs.emplace_back(gj);
                _priced_cols.select(n);
            }

            SPH_DEBUG {
                [[maybe_unused]] auto old_end = global_col_idxs.end();
                assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
            }
        }

        /**
         * @brief
         * Tell "Instance::fix_columns" what to do when a column contains a
         * row covered by another fixed column.
         * Set Partitioning: keep only columns that cover uncovered rows.
         *
         */
        [[nodiscard]] bool SetPar_ActiveColTest(const UniqueCol &col, const std::vector<bool> &active_rows_) {
            for (idx_t i : col) {
                if (!active_rows_[i]) {
                    return false;
                }  // discard
            }
            return true;  // keep
        }

        /**
         * @brief
         * Tell "Instance::fix_columns" what to do when a column contains a
         * row covered by another fixed column.
         * Set Covering: keep all the columns that contain an uncovered row.
         *
         */
        [[nodiscard]] bool SetCov_ActiveColTest(const UniqueCol &col, const std::vector<bool> &active_rows_) {
            for (idx_t i : col) {
                if (active_rows_[i]) {
                    return true;
                }  // keep
            }
            return true;  // discard
        }

        inline auto keepColStrategy(Column &col, std::vector<bool> &active_rows_) {
            if (keep_col_strat == SPP)
                return SetPar_ActiveColTest(col, active_rows_);
            else
                return SetCov_ActiveColTest(col, active_rows_);
        };

        void _fix_columns(const std::vector<idx_t> &idxs) {
            idx_t iok = 0;
            for (idx_t j = 0; j < cols.size(); ++j) {
                if (keepColStrategy(cols[j], active_rows)) {
                    active_cols[iok++] = j;
                }
            }

            active_cols.resize(iok);
            fixed_cols = idxs;

            fixed_cost = 0.0;
            for (idx_t j : fixed_cols) {
                fixed_cost += cols[j].get_cost();
            }
        }


    private:
        const idx_t nrows;
        UniqueColSet cols;
        std::vector<idx_t> active_cols;
        std::vector<idx_t> fixed_cols;
        std::vector<bool> active_rows;
        idx_t nactive_rows;
        real_t fixed_cost;
        idx_t ncols_constr;
        Priced_Columns priced_cols;
        MStar covering_times;
        idx_t inst_max_routes = INST_HARD_CAP;
        KeepStrat keep_col_strat = SPP;
        NewBestCallback new_best_callback;
        Timer timelimit;
    };
}  // namespace sph

#endif