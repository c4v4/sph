#ifndef SPH_INCLUDE_SUBINSTANCE_HPP_
#define SPH_INCLUDE_SUBINSTANCE_HPP_


#include "CollectionOf.hpp"
#include "Instance.hpp"
#include "MStar.hpp"
#include "Stopwatch.hpp"
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

    [[nodiscard]] bool is_corrupted() const;
    [[nodiscard]] real_t get_global_LB(const std::vector<real_t> &u_k);
    [[nodiscard]] idx_t find_local_col_idx(idx_t gj);
    [[nodiscard]] std::vector<idx_t> get_localized_solution(const std::vector<idx_t> &glob_sol);

    std::vector<idx_t> add_columns(const Cols &new_cols);

    inline void update_sol_cost(const std::vector<idx_t> &local_sol) {
        real_t sol_cost = std::accumulate(local_sol.begin(), local_sol.end(), 0.0, [&](real_t sum, idx_t lj) { return sum + cols[lj].get_cost(); });
        sol_cost += get_fixed_cost();
        update_sol_costs(local_sol, sol_cost);
    }

    void update_sol_costs(const std::vector<idx_t> &local_sol, real_t sol_cost);

    void reset();
    NO_INLINE real_t price(const std::vector<real_t> &u_k);
    idx_t fix_columns(std::vector<idx_t> &local_idxs_to_fix, std::vector<real_t> &u_star);
    void replace_columns(const std::vector<idx_t> &glob_cols_idxs);

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

    void _init_priced_cols(Priced_Columns &_priced_cols);
    real_t _price_active_cols(const std::vector<real_t> &u_k, Priced_Columns &_priced_cols);
    NO_INLINE void _select_C0_cols(Priced_Columns &_priced_cols, std::vector<idx_t> &global_col_idxs);
    NO_INLINE void _select_C1_cols(Priced_Columns &_priced_cols, MStar &_covering_times, std::vector<idx_t> &global_col_idxs);
    NO_INLINE void _select_C2_cols(Priced_Columns &_priced_cols, MStar &_covering_times, std::vector<idx_t> &global_col_idxs);

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
