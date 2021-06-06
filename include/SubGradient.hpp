#ifndef AC_CFT_INCLUDE_SUBGRADIENT_HPP_
#define AC_CFT_INCLUDE_SUBGRADIENT_HPP_

#include <cassert>
#include <random>

#include "LowerBound.hpp"
#include "SubGradientUtils.hpp"

#define REAL_TOLERANCE 1E-6

class SubGradient {

public:
    /**
     * Initialize a subgradient object by setting the initial multipliers
     * @param core
     */
    explicit SubGradient(SubInstance& subinst_) : subinst(subinst_) { }

    static LocalMultipliers u_greedy_init(const SubInstance& subinst) {

        LocalMultipliers u_0(subinst.get_nrows(), std::numeric_limits<real_t>::max());

        // init multipliers
        for (auto i = 0UL; i < subinst.get_nrows(); ++i) {
            for (const auto j : subinst.get_row(i)) {
                assert(!subinst.get_rows().empty());
                const auto candidate = subinst.get_col(j).get_cost() / static_cast<real_t>(subinst.get_col(j).size());
                if (candidate < u_0[i]) { u_0[i] = candidate; }
            }
        }

        return u_0;
    }

    static LocalMultipliers u_perturbed_init(const LocalMultipliers& u_star, std::mt19937& rnd) {

        auto u_0 = LocalMultipliers(u_star.size());
        auto dist = std::uniform_real_distribution<real_t>(-0.1, 0.1);

        idx_t u0size = u_0.size();
        for (idx_t i = 0; i < u0size; i++) { u_0[i] = (1 + dist(rnd)) * u_star[i]; }

        return u_0;
    }

    LocalMultipliers solve(real_t UB, const LocalMultipliers& u_0) {
        loop<false>(UB, u_0, 10 * subinst.get_rows().size());
        return u_star;
    }

    std::vector<LocalMultipliers>& explore(real_t UB, const LocalMultipliers& u_0, idx_t trials) {
        u_list.clear();
        loop<true>(UB, u_0, trials);
        return u_list;
    }

    real_t get_best_LB() { return LB_star; }

private:
    SubInstance& subinst;

    real_t LB_star;
    LocalMultipliers u_star;
    LocalMultipliers u;

    std::vector<LocalMultipliers> u_list;  // sequence of multipliers found during exploration of u_star
    ReducedLagrMultLB<1, 1000> norm_reducer;
    DeltaLowerBound<0, 1> lb_maintainer;

    StepSizeFactor lambda;

    template <bool heuristic_phase>
    void loop(real_t UB, const LocalMultipliers& u_0, idx_t max_iter) {

        const idx_t nrows = subinst.get_rows().size();

        // PARAMETERS
        if constexpr (!heuristic_phase) { lambda = StepSizeFactor(0.1, 20); }  // step size
        auto T = PricingPeriod(10, std::min<idx_t>(1000UL, nrows / 3));        // pricing frequency
        auto time_to_exit = ExitCondition(300U);

        LB_star = std::numeric_limits<real_t>::lowest();

        u = u_0;
        u_star = u_0;

        MStar covered_rows;
        LocalSolution S;
        real_t real_LB = lb_maintainer.compute(subinst, u);

        std::vector<std::pair<idx_t, real_t>> delta_u;

        for (idx_t iter = 0; iter < max_iter; ++iter) {

            lambda.update(real_LB);

            if constexpr (heuristic_phase) {
                norm_reducer.compute_sol(subinst, S, covered_rows);
            } else {
                norm_reducer.compute_reduced_sol(subinst, S, covered_rows);
            }

            // Multipliers update:
            delta_u.clear();
            idx_t s2sum = 0;
            for (idx_t cov : covered_rows) { s2sum += (1 - static_cast<int>(cov)) * (1 - static_cast<int>(cov)); }

            if (s2sum > 0) {
                for (idx_t i = 0; i < nrows; ++i) {
                    real_t new_u = std::max<real_t>(0.0, u[i] + lambda.get() * ((UB - real_LB) / s2sum) * (1 - static_cast<int>(covered_rows[i])));
                    if (std::abs(new_u - u[i]) > REAL_TOLERANCE) {
                        delta_u.emplace_back(i, new_u - u[i]);
                        u[i] = new_u;
                    }
                }
            } else {
                IF_VERBOSE { fmt::print(" WARNING: s2sum == 0\n"); }
                if constexpr (heuristic_phase) { u_list.emplace_back(u); }
                return;
            }

            real_LB = lb_maintainer.update(subinst, delta_u);

            if (real_LB > LB_star) {
                LB_star = real_LB;
                u_star = u;
            }

            if (real_LB > UB - HAS_INTEGRAL_COSTS) {
                IF_VERBOSE { fmt::print(" WARNING: real_LB({}) > UB({}) - {}\n", real_LB, UB, HAS_INTEGRAL_COSTS); }
                if constexpr (heuristic_phase) { u_list.emplace_back(u); }
                u_star = u;
                return;
            }

            if (covered_rows.get_uncovered() == 0) {
                real_t S_cost = S.compute_cost(subinst);
                if (S_cost < UB) { UB = S_cost; }
            }

            if constexpr (heuristic_phase) {
                u_list.emplace_back(u);
            } else {

                if (time_to_exit(LB_star)) {
                    IF_VERBOSE {
                        fmt::print("[{:^4}] Lower Bound: {:.4} (best {:.4}), lambda {:.4}, S_cost {:.4}\n", iter, real_LB, LB_star, lambda.get(),
                                   S.compute_cost(subinst));
                    }

                    return;
                }

                T.inc();
                if (T.reached()) {
                    const auto global_LB = subinst.price(u);
                    T.reset(global_LB, real_LB, UB);
                    // fmt::print("Pricing: global {}, local {}\n", global_LB, real_LB);

                    LB_star = real_LB = lb_maintainer.compute(subinst, u);
                }
            }
        }
    }
};

#endif
