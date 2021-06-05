#ifndef AC_CFT_INCLUDE_THREEPHASE_HPP
#define AC_CFT_INCLUDE_THREEPHASE_HPP

#include <fmt/core.h>

#include <algorithm>
#include <cassert>
#include <vector>

#include "ColumnFixing.hpp"
#include "Greedy.hpp"
#include "Instance.hpp"
#include "Multipliers.hpp"
#include "Solution.hpp"
#include "SubGradient.hpp"
#include "cft.hpp"

#define EXPLORING_ITERS 250

class ThreePhase {
public:
    ThreePhase(SubInstance& subinst_, MStar& covered_rows_, std::mt19937& rnd_)
        : subinst(subinst_),
          covered_rows(covered_rows_),
          subgradient(subinst),
          greedy(subinst, covered_rows_),
          col_fixing(subinst, greedy, covered_rows_),
          rnd(rnd_) { }

    std::pair<GlobalSolution, GlobalMultipliers> operator()(const real_t global_UB) {

        auto glo_u = GlobalMultipliers();
        auto u_star = LocalMultipliers(subinst.get_nrows(), 0.0);
        auto S_star = GlobalSolution(subinst, greedy(u_star));
        auto glo_UB_star = std::min(global_UB, S_star.get_cost());

        real_t u_k_LB;
        idx_t remaining_rows = subinst.get_nrows();

        auto u_k = SubGradient::u_greedy_init(subinst);
        idx_t iter = 1;

        do {

            IF_VERBOSE { fmt::print("┌─ 3-PHASE: iter {:2} ────────────────────────────────────────────────────────────────\n", iter); }

            auto S_curr = greedy(u_k);
            real_t S_curr_cost = S_curr.compute_cost(subinst);
            assert(S_curr_cost > 0);

            // 1. SUBGRADIENT PHASE
            u_k = subgradient.solve(S_curr_cost, u_k);  // use S_curr_cost because u_k refers to the current sub-problem

            // 2. HEURISTIC PHASE
            auto& u_list = subgradient.explore(S_curr_cost, u_k, EXPLORING_ITERS * iter);
            u_list.emplace_back(u_k);

            if (iter == 1) { glo_u = GlobalMultipliers(subinst, u_k); }

            auto fixed_cost = subinst.compute_fixed_cost();
            for (auto& u : u_list) {

                LocalSolution S = greedy(u);

                real_t S_cost = S.compute_cost(subinst);
                subinst.update_sol_cost(S, fixed_cost + S_cost);

                if (S_cost < S_curr_cost) {
                    S_curr = S;
                    S_curr_cost = S_cost;
                    u_star = u;
                    if (S_cost + fixed_cost < glo_UB_star) {
                        glo_UB_star = S_cost + fixed_cost;
                        S_star = GlobalSolution(subinst, S);
                        IF_VERBOSE { fmt::print("│ ══> [{:3}] Improved global UB: {}\n", &u - &u_list[0], S_star.get_cost()); }
                    } else {
                        IF_VERBOSE {
                            fmt::print("│ ──> [{:3}] Improved local UB: {} (global value {}, best is {})\n", &u - &u_list[0], S_cost,
                                       GlobalSolution(subinst, S).get_cost(), glo_UB_star);
                        }
                    }
                }
            }

            u_k_LB = subinst.get_global_LB(u_k);  // u_k.compute_lb(subinst);
            IF_VERBOSE { fmt::print("│ Active rows {}, fixed cost {}, sub-problem LB {}, current UB {}\n", remaining_rows, fixed_cost, u_k_LB, glo_UB_star); }

            if (subinst.compute_fixed_cost() + u_k_LB > glo_UB_star - 1.0) {
                IF_VERBOSE {
                    fmt::print("│ Early exit: LB >= UB - 1\n");
                    fmt::print("└───────────────────────────────────────────────────────────────────────────────────\n\n");
                }
                break;
            }

            // 3. COLUMN FIXING
            auto glo_S_curr = GlobalSolution(subinst, S_curr);
            remaining_rows = col_fixing(u_star, glo_S_curr);

            IF_VERBOSE { fmt::print("└───────────────────────────────────────────────────────────────────────────────────\n\n"); }

            u_k = SubGradient::u_perturbed_init(u_star, rnd);  // no u_star ma u_k credo che è il migliore LB prima del fixing?
            ++iter;
        } while (remaining_rows > 0 /*&& subinst.compute_fixed_cost() + u_k_LB < glo_UB_star*/);

        return std::make_pair(S_star, glo_u);
    }

private:
    SubInstance& subinst;
    MStar& covered_rows;

    SubGradient subgradient;
    Greedy greedy;
    ColumnFixing col_fixing;

    std::mt19937& rnd;
};


#endif