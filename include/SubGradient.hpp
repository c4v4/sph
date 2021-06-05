#ifndef AC_CFT_INCLUDE_SUBGRADIENT_HPP_
#define AC_CFT_INCLUDE_SUBGRADIENT_HPP_

#include <cassert>
#include <random>

#include "LowerBound.hpp"
#include "SubGradientUtils.hpp"

#define REAL_TOLERANCE 1E-10

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

private:
    SubInstance& subinst;

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

        auto LB = std::numeric_limits<real_t>::lowest();

        u = u_0;
        u_star = u_0;

        MStar covered_rows;
        LocalSolution S;
        real_t real_LB = lb_maintainer.compute(subinst, u);
        // const auto& s = subinst;
        // fmt::print("INIT: Old {}, new {}\n", lagr_mul_LB(s, u), real_LB);

        std::vector<std::pair<idx_t, real_t>> delta_u;

        for (idx_t iter = 0; iter < max_iter; ++iter) {

            lambda.update(real_LB);

            if constexpr (heuristic_phase) {
                norm_reducer.compute_sol(subinst, S, covered_rows);
            } else {
                norm_reducer.compute_reduced_sol(subinst, S, covered_rows);
                // fmt::print(stderr, "S size {} (rows: {})\n", S.size(), subinst.get_nrows());
            }

            // Multipliers update:
            delta_u.clear();
            real_t s2sum = 0.0;
            for (auto cov : covered_rows) { s2sum += static_cast<real_t>((1.0 - cov) * (1.0 - cov)); }
            for (idx_t i = 0; i < nrows; ++i) {
                real_t delta = lambda.get() * ((UB - real_LB) / s2sum) * (1.0 - covered_rows[i]);
                if (std::abs(delta) > REAL_TOLERANCE) {
                    u[i] = std::max<real_t>(0.0, u[i] + delta);
                    delta_u.emplace_back(i, delta);
                }
            }

            real_LB = lb_maintainer.update(subinst, delta_u);
            // real_LB = lb_maintainer.compute(subinst, u);
            // fmt::print("[{}]: old {}, new {}\n", iter, lagr_mul_LB(subinst, u), real_LB);

            if (real_LB >= UB) {
                fmt::print(" WARNING: real_LB > UB\n");
                u_star = u;
                return;
            }

            if (real_LB > LB) {
                LB = real_LB;
                u_star = u;
            }

            if (covered_rows.get_uncovered() == 0) {
                real_t S_cost = S.compute_cost(subinst) + subinst.compute_fixed_cost();
                if (S_cost < UB) { UB = S_cost; }
            }

            // fmt::print("[{:^4}] Lower Bound: {:.4} (best {:.4}), lambda {:.4}, S_cost {:.4}\n", iter, real_LB, LB, lambda.get(), S.compute_cost(subinst));

            if constexpr (heuristic_phase) {
                u_list.emplace_back(u);
            } else {

                if (time_to_exit(LB)) {
                    IF_VERBOSE {
                        fmt::print("[{:^4}] Lower Bound: {:.4} (best {:.4}), lambda {:.4}, S_cost {:.4}\n", iter, real_LB, LB, lambda.get(),
                                   S.compute_cost(subinst));
                    }

                    return;
                }

                T.inc();
                if (T.reached()) {
                    const auto global_LB = subinst.price(u);
                    T.reset(global_LB, real_LB, UB);
                    // fmt::print("Pricing: global {}, local {}\n", global_LB, real_LB);

                    LB = real_LB = lb_maintainer.compute(subinst, u);
                }
            }
        }
    }
};

#endif  // AC_CFT_INCLUDE_SUBGRADIENT_HPP_
