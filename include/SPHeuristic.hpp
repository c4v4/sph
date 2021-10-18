#ifndef SPH_INCLUDE_SPHEURISTIC_HPP_
#define SPH_INCLUDE_SPHEURISTIC_HPP_

#include <cassert>

#include "Instance.hpp"
#include "Refinement.hpp"

namespace sph {

    class SPHeuristic {
    public:
        explicit SPHeuristic(const idx_t nrows_) : inst(nrows_), rnd(std::mt19937()), refinement(inst, rnd) { }
        SPHeuristic(const idx_t nrows_, int seed) : inst(nrows_), rnd(std::mt19937(seed)), refinement(inst, rnd) { }

        [[nodiscard]] inline idx_t get_ncols() const { return inst.get_ncols(); }
        [[nodiscard]] inline idx_t get_nrows() const { return inst.get_nrows(); }
        [[nodiscard]] inline UniqueColSet &get_cols() { return inst.get_cols(); }
        [[nodiscard]] inline Column &get_col(idx_t idx) { return inst.get_col(idx); }
        [[nodiscard]] inline const Column &get_col(idx_t idx) const { return inst.get_col(idx); }

        void inline set_timelimit(double seconds) { inst.set_timelimit(seconds); }

        [[nodiscard]] inline Timer &get_timelimit() { return inst.get_timelimit(); }

        template <typename UniqueColContainer>
        std::vector<idx_t> inline add_columns(const UniqueColContainer &new_cols) {
            return inst.add_columns(new_cols);
        }

        std::vector<idx_t> inline add_columns(const std::vector<real_t> &costs, const std::vector<real_t> &sol_costs, const std::vector<idx_t> &matbeg,
                                              const std::vector<idx_t> &matval) {
            return inst.add_columns(costs, sol_costs, matbeg, matval);
        }

        template <typename... _Args>
        idx_t add_column(_Args &&...args) {
            return inst.add_column(std::forward<_Args>(args)...);
        }

        template <unsigned long ROUTES_HARD_CAP = SPH_INST_HARD_CAP>
        std::vector<idx_t> inline solve([[maybe_unused]] const std::vector<idx_t> &S_init) {
            SPH_VERBOSE(0) { fmt::print(" Set Partitioning Heuristic: \n"); }
            SPH_VERBOSE(0) { fmt::print(" SP Instance size: {}x{}\n", inst.get_nrows(), inst.get_ncols()); }
            
            GlobalSolution sol = refinement.solve<ROUTES_HARD_CAP>(S_init);
            
            SPH_VERBOSE(0) { fmt::print(" Final solution value: {}\n", sol.get_cost()); }
            return sol;
        }

    private:
        Instance inst;
        std::mt19937 rnd;
        Refinement refinement;
    };

}  // namespace sph

#endif