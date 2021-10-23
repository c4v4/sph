#ifndef SPH_INCLUDE_REFINEMENT_HPP_
#define SPH_INCLUDE_REFINEMENT_HPP_

#include "Instance.hpp"
#include "SubInstance.hpp"
#include "TwoPhase.hpp"
#include "cft.hpp"

namespace sph {

    constexpr double ALPHA = 1.1;
    constexpr double BETA = 0.8;
    constexpr double PI_MIN = 0.3;

    constexpr double PI_MAX = 0.9;
    constexpr unsigned TRIALS = 1000;

    constexpr double SHORT_T_LIM = 10.0;

#define LONG_T_LIM(TOTAL_TIME) (std::min(TOTAL_TIME / 2.0, 100.0))

    class Refinement {
    public:
        Refinement(Instance& inst_, std::mt19937& rnd_) : inst(inst_), subinst(inst_), two_phase(subinst), rnd(rnd_) { }

        // S_init must be a global complete solution
        template <unsigned long ROUTES_HARD_CAP, typename KeepColStrategy = SetPar_ActiveColTest>
        GlobalSolution solve([[maybe_unused]] const std::vector<idx_t>& S_init) {

            SPH_VERBOSE(1) {
                if constexpr (std::is_same_v<KeepColStrategy, SetPar_ActiveColTest>) {
                    fmt::print("  Using Set Partitioning fixing strategy (overlaps are forbidden).\n");
                } else if (std::is_same_v<KeepColStrategy, SetCov_ActiveColTest>) {
                    fmt::print("  Using Set Covering fixing strategy (overlaps are alowed).\n");
                } else {
                    fmt::print(stderr, "  Warning: Using unkown fixing strategy.\n");
                }
            }

            inst.reset_fixing();

            GlobalSolution S_star;
            if (!S_init.empty()) {
                real_t cost = 0.0;
                for (idx_t j : S_init) { cost += inst.get_col(j).get_cost(); }
                for (auto j : S_init) { S_star.push_back(j); }
                S_star.set_cost(cost);
                SPH_VERBOSE(1) {
                    fmt::print("  Found warm start of cost {}.\n", cost);
                    SPH_VERBOSE(2) {
                        for (idx_t gj : S_star) {
                            Column col = inst.get_col(gj);
                            fmt::print("   idx: {}, cost: {}, sol cost: {}\n", gj, col.get_cost(), col.get_solcost());
                        }
                        fmt::print("\n");
                    }
                }
            }

            GlobalMultipliers u_star;

            real_t pi = PI_MIN;
            real_t last_improving_pi = pi;

            unsigned trial = 0;
            Timer& global_time_limit = inst.get_timelimit();

            idx_t iter = 1;
            do {
                subinst.reset();  // 2.

                {
                    Timer iteration_timer = Timer(iter == 1 ? LONG_T_LIM(global_time_limit.seconds_until_end()) : SHORT_T_LIM);
                    GlobalSolution S = two_phase(S_star.get_cost(), S_star, iteration_timer);  // 3. & 4.

                    assert(!(std::fabs(pi - PI_MIN) > 0.001 && inst.get_fixed_cols().empty()));
                    pi *= ALPHA;  // 6.

                    if (S.get_cost() < S_star.get_cost()) {  // update best solution
                        S_star = std::move(S);               // 5.

                        last_improving_pi = pi;

                        pi = std::max(pi / (ALPHA * ALPHA), PI_MIN);  // 6.
                    }

                    if (S_star.get_cost() - 1.0 <= BETA * u_star.get_lb() || pi > PI_MAX || inst.get_active_rows_size() <= 0) {

                        if (trial++ > TRIALS) {
                            SPH_VERBOSE(2) {
                                fmt::print("   ╔═ REFINEMENT: iter {:3} ══════════════════════════════════════════════════\n", iter);
                                fmt::print("   ║ Early Exit: β(={:.1f}) * LB(={:.1f}) > UB(={:.1f}) - 1\n", BETA, u_star.get_lb(),
                                           S_star.get_cost());
                                fmt::print("   ║ Active rows {}, active cols {}, pi {:.3f}\n", inst.get_active_rows_size(),
                                           inst.get_active_cols().size(), pi);
                                fmt::print("   ║ LB {:.1f}, UB {:.1f}, UB size {}\n", u_star.get_lb(), S_star.get_cost(), S_star.size());
                                fmt::print("   ╚═════════════════════════════════════════════════════════════════════════\n\n");
                            }
                            break;
                        }

                        SPH_VERBOSE(1) { fmt::print("  Iteration: {:2}; Best: {:.1f} \n", trial, S_star.get_cost()); }

                        pi = last_improving_pi;
                        last_improving_pi = std::max(PI_MIN, last_improving_pi / ALPHA);
                    }
                }  //(S and u are potentially been moved, better to encapsulate them into a block)

                if (iter == 1) {
                    u_star = two_phase.get_global_u();
                    std::vector<idx_t> old_to_new_idx_map = inst.prune_instance<ROUTES_HARD_CAP>(u_star);

                    SPH_VERBOSE(1) { fmt::print("  Instance size: {}x{}\n", inst.get_nrows(), inst.get_ncols()); }

                    if (!old_to_new_idx_map.empty()) {
                        for (idx_t& gj : S_star) {
                            gj = old_to_new_idx_map[gj];
                            assert(gj != NOT_AN_INDEX);
                        }
                    }
                }

                inst.reset_fixing();

                if (global_time_limit.exceeded_tlim()) {
                    SPH_VERBOSE(2) {
                        fmt::print("   ╔═ REFINEMENT: iter {:3} ══════════════════════════════════════════════════\n", iter);
                        fmt::print("   ║ Timelimit exceeded\n");
                        fmt::print("   ╚═════════════════════════════════════════════════════════════════════════\n\n");
                    }
                    break;
                }

                pi = std::min<real_t>(PI_MAX, pi);

                // 7. Refinement Fix
                cols_to_fix = random_fix(S_star, pi);

                assert(cols_to_fix.size() <= S_star.size());

                std::sort(cols_to_fix.begin(), cols_to_fix.end());
                inst.fix_columns<KeepColStrategy>(cols_to_fix, covered_rows);

                SPH_VERBOSE(2) {
                    fmt::print("   ╔═ REFINEMENT: iter {:3} ══════════════════════════════════════════════════\n", iter);
                    fmt::print("   ║ Active rows {}, active cols {}, pi {:.3f}\n", inst.get_active_rows_size(),
                               inst.get_active_cols().size(), pi);
                    fmt::print("   ║ LB {:.1f}, UB {:.1f}, UB size {}\n", u_star.get_lb(), S_star.get_cost(), S_star.size());
                    fmt::print("   ╚═════════════════════════════════════════════════════════════════════════\n\n");
                }

                assert(inst.get_fixed_cost() <= S_star.get_cost());
                ++iter;
            } while (true);

            return S_star;
        }

    private:
        [[nodiscard]] std::vector<idx_t> refinement_fix(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi) {
            UniqueColSet& cols = inst.get_cols();

            covered_rows.reset_covered(cols, S_star, inst.get_nrows());

            deltas.resize(S_star.size());
            for (idx_t j = 0; j < S_star.size(); ++j) {
                auto& col = cols[S_star[j]];
                deltas[j].first = S_star[j];
                deltas[j].second = std::max<real_t>(col.compute_lagr_cost(u_star), 0.0);
                for (auto i : col) { deltas[j].second += u_star[i] * (covered_rows[i] - 1.0) / covered_rows[i]; }
            }
            std::sort(deltas.begin(), deltas.end(), [](auto& a, auto& b) { return a.second < b.second; });

            covered_rows.reset_uncovered(inst.get_nrows());
            real_t covered_fraction = 0.0;
            idx_t n = 0;
            for (; n < deltas.size() && covered_fraction < pi; ++n) {
                covered_rows.cover_rows(cols[deltas[n].first]);
                covered_fraction = static_cast<real_t>(covered_rows.get_covered()) / static_cast<real_t>(inst.get_nrows());
            }

            cols_to_fix.resize(n);
            for (idx_t j2 = 0; j2 < n; ++j2) { cols_to_fix[j2] = deltas[j2].first; }

            return cols_to_fix;
        }

        [[nodiscard]] std::vector<idx_t> random_fix(GlobalSolution S_star, real_t pi) {
            UniqueColSet& cols = inst.get_cols();

            std::shuffle(S_star.begin(), S_star.end(), rnd);

            cols_to_fix.clear();
            covered_rows.reset_uncovered(inst.get_nrows());
            real_t covered_fraction = 0.0;

            for (idx_t n = 0; n < S_star.size() && covered_fraction < pi; ++n) {
                cols_to_fix.emplace_back(S_star[n]);
                covered_rows.cover_rows(cols[S_star[n]]);
                covered_fraction = static_cast<real_t>(covered_rows.get_covered()) / static_cast<real_t>(inst.get_nrows());
            }

            return cols_to_fix;
        }

        [[nodiscard]] std::vector<idx_t> binary_tournament_fix(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi) {
            UniqueColSet& cols = inst.get_cols();

            covered_rows.reset_covered(cols, S_star, inst.get_nrows());

            deltas.resize(S_star.size());
            for (idx_t j = 0; j < S_star.size(); ++j) {
                auto& col = cols[S_star[j]];
                deltas[j].first = S_star[j];
                deltas[j].second = std::max<real_t>(col.compute_lagr_cost(u_star), 0.0);
                for (auto i : col) { deltas[j].second += u_star[i] * (covered_rows[i] - 1.0) / covered_rows[i]; }
            }

            idx_t dsize = deltas.size();

            cols_to_fix.resize(dsize);
            covered_rows.reset_uncovered(inst.get_nrows());
            real_t covered_fraction = 0.0;
            idx_t n = 0;
            for (; n < deltas.size() && dsize > 1 && covered_fraction < pi; ++n, --dsize) {
                auto& cand1 = deltas[rnd() % dsize];
                idx_t c2;
                do { c2 = rnd() % dsize; } while (cand1.first == deltas[c2].first);
                auto& winner = cand1.second < deltas[c2].second ? cand1 : deltas[c2];

                cols_to_fix[n] = winner.first;
                covered_rows.cover_rows(cols[winner.first]);
                covered_fraction = static_cast<real_t>(covered_rows.get_covered()) / static_cast<real_t>(inst.get_nrows());

                std::swap(deltas.back(), winner);
            }
            cols_to_fix.resize(n);

            return cols_to_fix;
        }

        [[nodiscard]] std::vector<idx_t> random_fix2(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi) {
            UniqueColSet& cols = inst.get_cols();
            std::uniform_real_distribution<real_t> dist(0.5, 1.5);

            covered_rows.reset_covered(cols, S_star, inst.get_nrows());

            deltas.resize(S_star.size());
            for (idx_t j = 0; j < S_star.size(); ++j) {
                auto& col = cols[S_star[j]];
                deltas[j].first = S_star[j];
                deltas[j].second = std::max<real_t>(col.compute_lagr_cost(u_star), 0.0);
                for (auto i : col) { deltas[j].second += u_star[i] * (covered_rows[i] - 1.0) / covered_rows[i]; }
                deltas[j].second *= dist(rnd);
            }
            std::sort(deltas.begin(), deltas.end(), [](auto& a, auto& b) { return a.second < b.second; });

            covered_rows.reset_uncovered(inst.get_nrows());
            real_t covered_fraction = 0.0;
            idx_t n = 0;
            for (; n < deltas.size() && covered_fraction < pi; ++n) {
                covered_rows.cover_rows(cols[deltas[n].first]);
                covered_fraction = static_cast<real_t>(covered_rows.get_covered()) / static_cast<real_t>(inst.get_nrows());
            }

            cols_to_fix.resize(n);
            for (idx_t j2 = 0; j2 < n; ++j2) { cols_to_fix[j2] = deltas[j2].first; }

            return cols_to_fix;
        }

    private:
        Instance& inst;

        SubInstance subinst;
        MStar covered_rows;
        TwoPhase two_phase;

        std::mt19937& rnd;

        // retain allocated memory (anti-RAII, should be local)
        std::vector<std::pair<idx_t, real_t>> deltas;
        std::vector<idx_t> cols_to_fix;
    };

}  // namespace sph

#endif