#ifndef SPH_INCLUDE_INSTANCE_HPP_
#define SPH_INCLUDE_INSTANCE_HPP_

#include <cassert>

#include "CollectionOf.hpp"
#include "MStar.hpp"
#include "Stopwatch.hpp"
#include "cft.hpp"

#define MIN_COV 4U
#define MIN_SOLCOST_COV 4U
#define HARD_CAP 15000U

class Instance {
public:
    explicit Instance(const idx_t nrows_) : nrows(nrows_), active_rows(nrows, true), nactive_rows(nrows), fixed_cost(0.0) { }

    [[nodiscard]] inline auto get_ncols() const { return cols.size(); }
    [[nodiscard]] inline idx_t get_nrows() const { return nrows; }
    [[nodiscard]] inline idx_t get_active_rows_size() const { return nactive_rows; }

    [[nodiscard]] inline auto &get_active_cols() { return active_cols; }
    [[nodiscard]] inline auto &get_fixed_cols() { return fixed_cols; }
    [[nodiscard]] inline auto &get_cols() { return cols; }
    [[nodiscard]] inline auto &get_col(idx_t idx) { return cols[idx]; }
    [[nodiscard]] inline const auto &get_col(idx_t idx) const { return cols[idx]; }
    [[nodiscard]] inline real_t get_fixed_cost() { return fixed_cost; }

    inline void set_timelimit(double seconds) { timelimit = Timer(seconds); }
    [[nodiscard]] inline Timer &get_timelimit() { return timelimit; }

    [[nodiscard]] inline bool is_row_active(idx_t gi) {
        assert(gi < nrows);
        return active_rows[gi];
    }

    void inline fix_columns(const std::vector<idx_t> &idxs) {
        for (idx_t j : idxs) {
            for (idx_t i : cols[j]) { active_rows[i] = false; }
        }

        _fix_columns(idxs);
    }

    void reset_fixing();
    void fix_columns(const std::vector<idx_t> &idxs, const MStar &M_star);
    std::vector<idx_t> add_columns(const Cols &new_cols);
    std::vector<idx_t> add_columns(const std::vector<real_t> &costs, const std::vector<real_t> &sol_costs, const std::vector<idx_t> &matbeg,
                                   const std::vector<idx_t> &matval);

private:
    void _fix_columns(const std::vector<idx_t> &idxs);


private:
    const idx_t nrows;
    Cols cols;
    std::vector<idx_t> active_cols;
    std::vector<idx_t> fixed_cols;
    std::vector<bool> active_rows;
    idx_t nactive_rows;
    real_t fixed_cost;

    Timer timelimit;
};

#endif