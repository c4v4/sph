#ifndef SPH_INCLUDE_SPHEURISTIC_HPP_
#define SPH_INCLUDE_SPHEURISTIC_HPP_

#include <cassert>

#include "Instance.hpp"
#include "Refinement.hpp"

class SPHeuristic {
public:
    explicit SPHeuristic(const idx_t nrows_) : inst(nrows_), rnd(std::mt19937()), refinement(inst, rnd) { }
    explicit SPHeuristic(const idx_t nrows_, int seed) : inst(nrows_), rnd(std::mt19937(seed)), refinement(inst, rnd) { }

    [[nodiscard]] inline auto get_ncols() const { return inst.get_ncols(); }
    [[nodiscard]] inline idx_t get_nrows() const { return inst.get_nrows(); }
    [[nodiscard]] inline auto &get_cols() { return inst.get_cols(); }
    [[nodiscard]] inline auto &get_col(idx_t idx) { return inst.get_col(idx); }
    [[nodiscard]] inline const auto &get_col(idx_t idx) const { return inst.get_col(idx); }

    void inline set_timelimit(double seconds) { inst.set_timelimit(seconds); }
    [[nodiscard]] inline Timer &get_timelimit() { return inst.get_timelimit(); }

    std::vector<idx_t> inline add_columns(const Cols &new_cols) { return inst.add_columns(new_cols); }
    std::vector<idx_t> inline add_columns(const std::vector<real_t> &costs, const std::vector<real_t> &sol_costs, const std::vector<idx_t> &matbeg,
                                          const std::vector<idx_t> &matval) {
        return inst.add_columns(costs, sol_costs, matbeg, matval);
    }

    std::vector<idx_t> inline solve([[maybe_unused]] const std::vector<idx_t> &S_init) { return refinement(S_init); }

private:
    Instance inst;
    std::mt19937 rnd;
    Refinement refinement;
};

#endif