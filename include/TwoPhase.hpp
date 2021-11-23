#ifndef SPH_INCLUDE_TWOPHASE_HPP_
#define SPH_INCLUDE_TWOPHASE_HPP_

#include <algorithm>
#include <cassert>
#include <vector>

#include "ExactSolver.hpp"
#include "Multipliers.hpp"
#include "Solution.hpp"
#include "SubGradient.hpp"
#include "SubInstance.hpp"
#include "cft.hpp"
#include "fmt/core.h"

namespace sph {

    class TwoPhase {
    public:
        TwoPhase(SubInstance& subinst_) : subinst(subinst_), subgradient(subinst_) { }

        inline GlobalSolution solve(const real_t global_UB, const GlobalSolution& S_star, const Timer& exact_time_limit) {
            return operator()(global_UB, S_star, exact_time_limit);
        }

        GlobalSolution operator()(const real_t global_UB, const GlobalSolution& S_star, const Timer& exact_time_limit) {

            LocalMultipliers u_star(subinst.get_nrows(), 0.0);

            real_t glb_UB_star = std::min<real_t>(global_UB, S_star.get_cost());
            real_t fixed_cost = subinst.get_fixed_cost();
            real_t subgrad_UB = glb_UB_star - fixed_cost;

            // 1. SUBGRADIENT PHASE
            LocalMultipliers u_k = subgradient.solve(subgrad_UB, SubGradient::u_greedy_init(subinst), subinst.get_timelimit());

            real_t lcl_LB = subgradient.get_best_LB();
            real_t glb_LB = fixed_cost + lcl_LB;

            if (fixed_cost == 0.0) {
                glo_u = GlobalMultipliers(subinst, u_k);
            }
            if (glb_LB >= glb_UB_star - HAS_INTEGRAL_COSTS || subinst.get_nrows() == 0) {
                return GlobalSolution();
            }

            // 2. HEURISTIC PHASE
            LocalSolution S(subinst.get_localized_solution(S_star));
            S = cplex_heur(S, exact_time_limit);
            return GlobalSolution(subinst, S);
        }

        inline GlobalMultipliers& get_global_u() { return glo_u; }

    private:
        LocalSolution cplex_heur(const LocalSolution& S_init, const Timer& exact_time_limit) {

            real_t fixed_cost = subinst.get_fixed_cost();
            real_t S_init_cost = S_init.compute_cost(subinst);
            SPH_VERBOSE(3) { fmt::print("    │ Initial solution value {} (global: {})\n", S_init_cost, S_init_cost + fixed_cost); }

            return exact.build_and_opt(subinst, S_init, exact_time_limit);
        }


    private:
        SubInstance& subinst;

        SubGradient subgradient;
        ExactSolver exact;

        GlobalMultipliers glo_u;
    };
}  // namespace sph

#endif