#ifndef SPH_INCLUDE_TWOPHASE_HPP_
#define SPH_INCLUDE_TWOPHASE_HPP_

#include <fmt/core.h>

#include <algorithm>
#include <cassert>
#include <vector>

#include "ExactSolver.hpp"
#include "Multipliers.hpp"
#include "Solution.hpp"
#include "SubGradient.hpp"
#include "SubInstance.hpp"
#include "cft.hpp"

#define SHORT_T_LIM 10.0
#define LONG_T_LIM(TOTAL_TIME) (TOTAL_TIME / 2.0)

class TwoPhase {
public:
    TwoPhase(SubInstance& subinst_) : subinst(subinst_), subgradient(subinst_), exact_time_limit(LONG_T_LIM(subinst.get_timelimit().seconds_until_end())) { }
    inline GlobalSolution solve(const real_t global_UB, GlobalSolution& S_star) { return operator()(global_UB, S_star); }
    GlobalSolution operator()(const real_t global_UB, GlobalSolution& S_star) {

        LocalMultipliers u_star(subinst.get_nrows(), 0.0);

        real_t glb_UB_star = std::min<real_t>(global_UB, S_star.get_cost());
        real_t fixed_cost = subinst.get_fixed_cost();
        real_t subgrad_UB = glb_UB_star - fixed_cost;

        // 1. SUBGRADIENT PHASE
        LocalMultipliers u_k = subgradient.solve(subgrad_UB, SubGradient::u_greedy_init(subinst));

        real_t lcl_LB = subgradient.get_best_LB();
        real_t glb_LB = fixed_cost + lcl_LB;

        if (fixed_cost == 0.0) { glo_u = GlobalMultipliers(subinst, u_k); }
        if (glb_LB >= glb_UB_star - HAS_INTEGRAL_COSTS) { return S_star; }

        // 2. HEURISTIC PHASE
        LocalSolution S_curr(subinst.get_localized_solution(S_star));
        cplex_heur(S_curr, S_star, glb_UB_star);

        return S_star;
    }

    inline GlobalMultipliers& get_global_u() { return glo_u; }

private:
    void cplex_heur(LocalSolution& S_init, GlobalSolution& S_star, real_t& glb_UB_star) {

        real_t fixed_cost = subinst.get_fixed_cost();
        real_t S_init_cost = S_init.compute_cost(subinst);
        IF_VERBOSE { fmt::print("│ Initial solution value {} (global: {})\n", S_init_cost, S_init_cost + fixed_cost); }

        LocalSolution S = exact.build_and_opt(subinst, S_init, exact_time_limit);
        exact_time_limit = SHORT_T_LIM;

        if (S.size() > 0) {

            real_t S_cost = S.compute_cost(subinst);

            real_t gS_cost = fixed_cost + S_cost;
            subinst.update_sol_costs(S, gS_cost);

            if (gS_cost < glb_UB_star) {

                glb_UB_star = gS_cost;
                S_star = GlobalSolution(subinst, S);
                IF_VERBOSE { fmt::print("│ ══> CPLEX improved global UB: {} (fixed {} + local-cost {})\n", S_star.get_cost(), fixed_cost, S_cost); }

            } else {
                IF_VERBOSE { fmt::print("│ ──> CPLEX Improved local UB: {} (global value {}, best is {})\n", S_cost, S_cost + fixed_cost, glb_UB_star); }
            }
        }

        IF_VERBOSE { fmt::print("└───────────────────────────────────────────────────────────────────────────────────\n\n"); }
    }


private:
    SubInstance& subinst;

    SubGradient subgradient;
    ExactSolver exact;

    double exact_time_limit;

    GlobalMultipliers glo_u;
};


#endif