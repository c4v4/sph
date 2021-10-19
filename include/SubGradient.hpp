#ifndef SPH_INCLUDE_SUBGRADIENT_HPP_
#define SPH_INCLUDE_SUBGRADIENT_HPP_

#include <cassert>
#include <random>

#include "LowerBound.hpp"
#include "SubGradientUtils.hpp"
#include "cft.hpp"

namespace sph {


    constexpr double REAL_TOLERANCE = 1E-6;

    class SubGradient {

    public:
        /**
         * Initialize a subgradient object by setting the initial multipliers
         * @param core
         */
        explicit SubGradient(SubInstance& subinst_) : subinst(subinst_) { }

        static LocalMultipliers u_greedy_init(const SubInstance& subinst) {

            LocalMultipliers u_0(subinst.get_nrows(), REAL_MAX);

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
            auto dist = std::uniform_real_distribution<real_t>(0.9, 1.1);

            idx_t u0size = u_0.size();
            for (idx_t i = 0; i < u0size; i++) { u_0[i] = dist(rnd) * u_star[i]; }

            return u_0;
        }

        LocalMultipliers solve(real_t UB, const LocalMultipliers& u_0, Timer& time_limit) {
            size_t max_iter = 10 * subinst.get_rows().size();
            idx_t nrows = subinst.get_rows().size();

            // PARAMETERS
            lambda = StepSizeFactor(0.1, 20);                         // step size
            PricingPeriod T(10, std::min<idx_t>(1000UL, nrows / 3));  // pricing frequency
            ExitCondition exit_now(300U);

            LB_star = REAL_LOWEST;
            u = u_0;
            u_star = u_0;

            MStar covered_rows;
            LocalSolution S;
            real_t real_LB = lb_maintainer.compute(subinst, u);

            std::vector<std::pair<idx_t, real_t>> delta_u;

            for (idx_t iter = 0; iter < max_iter; ++iter) {
                if (time_limit.exceeded_tlim()) { break; }

                lambda.update(real_LB);
                norm_reducer.compute_reduced_sol(subinst, S, covered_rows);

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
                    SPH_VERBOSE(3) { fmt::print("    WARNING: s2sum == 0\n"); }
                    return u_star;
                }

                real_LB = lb_maintainer.update(subinst, delta_u);

                if (real_LB > LB_star) {
                    LB_star = real_LB;
                    u_star = u;
                }

                // fmt::print("[{:^4}] Lower Bound: {:.4} (best {:.4}), lambda {:.4}, S_cost {:.4}\n", iter, real_LB, LB_star, lambda.get(),
                // S.compute_cost(subinst));

                if (real_LB > UB - HAS_INTEGRAL_COSTS) {
                    SPH_VERBOSE(3) { fmt::print("    WARNING: real_LB({}) > UB({}) - {}\n", real_LB, UB, HAS_INTEGRAL_COSTS); }
                    u_star = u;
                    return u_star;
                }

                if (covered_rows.get_uncovered() == 0) {
                    real_t S_cost = S.compute_cost(subinst);
                    if (S_cost < UB) { UB = S_cost; }
                }

                if (exit_now(LB_star)) { return u_star; }

                T.inc();
                if (T.reached()) {
                    const auto global_LB = subinst.price(u);
                    T.reset(global_LB, real_LB, UB);
                    // fmt::print("Pricing: global {}, local {}\n", global_LB, real_LB);

                    LB_star = real_LB = lb_maintainer.compute(subinst, u);
                }
            }
            return u_star;
        }


        real_t get_best_LB() { return LB_star; }

    private:
        SubInstance& subinst;

        real_t LB_star;
        LocalMultipliers u_star;
        LocalMultipliers u;

        ReducedLagrMultLB<1, 1000> norm_reducer;
        DeltaLowerBound<0, 1> lb_maintainer;

        StepSizeFactor lambda;
    };

}  // namespace sph

#endif
