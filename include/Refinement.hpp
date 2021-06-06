#ifndef AC_CFT_INCLUDE_REFINEMENT_HPP_
#define AC_CFT_INCLUDE_REFINEMENT_HPP_

#include "Instance.hpp"
#include "ThreePhase.hpp"
#include "cft.hpp"

#define ALPHA (1.1)
#define BETA (1.0)
#define PI_MIN (0.3)

class Refinement {
public:
    Refinement(Instance& inst_, std::mt19937& rnd) : inst(inst_), subinst(inst_), three_phase(subinst, covered_rows, rnd), pi(PI_MIN) { }

    // S must be a global complete solution
    std::vector<idx_t> operator()([[maybe_unused]] const std::vector<idx_t>& S_init) {

        // 1.
        GlobalSolution S_star;

        if (!S_init.empty()) {
            real_t cost = 0.0;
            for (idx_t j : S_init) { cost += inst.get_col(j).get_cost(); }
            for (auto j : S_init) { S_star.push_back(j); }
            S_star.set_cost(cost);
            IF_VERBOSE { fmt::print("Found warm start with cost {}.\n", cost); }
        }

        GlobalMultipliers u_star;

        pi = PI_MIN;

        idx_t iter = 1;
        do {
            // 2.
            subinst.reset();

            {
                // 3. & 4.
                auto [S, u] = three_phase(S_star.get_cost());

                // 6.
                assert(!(std::fabs(pi - PI_MIN) > 0.001 && inst.get_fixed_cols().empty()));
                pi *= ALPHA;

                // update best solution
                if (S.get_cost() < S_star.get_cost()) {
                    S_star = std::move(S);  // 5.

                    // pi = PI_MIN;            // 6.
                    pi = std::max(pi / (ALPHA * ALPHA), PI_MIN);
                }

                // ?. update best lower bound
                if (u.get_lb() > u_star.get_lb()) { u_star = std::move(u); }
            }  //(S and u are potentially been moved, better to encapsulate them into a block)

            // 7. Refinement Fix
            inst.reset_fixing();
            auto& cols = inst.get_cols();
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
            idx_t j1 = 0;
            for (; j1 < deltas.size() && covered_fraction < pi; ++j1) {
                covered_rows.cover_rows(cols[deltas[j1].first]);
                covered_fraction = static_cast<real_t>(covered_rows.get_covered()) / static_cast<real_t>(inst.get_nrows());
            }

            cols_to_fix.resize(j1);
            for (idx_t j2 = 0; j2 < j1; ++j2) { cols_to_fix[j2] = deltas[j2].first; }

            assert(cols_to_fix.size() <= S_star.size());

            std::sort(cols_to_fix.begin(), cols_to_fix.end());
            inst.fix_columns(cols_to_fix, covered_rows);

            IF_VERBOSE {
                fmt::print("╔═ REFINEMENT: iter {:2} ═════════════════════════════════════════════════════════════\n", iter);
                fmt::print("║ Active rows {}, active cols {}, pi {}\n", inst.get_active_rows_size(), inst.get_active_cols().size(), pi);
                fmt::print("║ Best sol {}, cols {}\n", S_star.get_cost(), S_star.size());
                fmt::print("╚═══════════════════════════════════════════════════════════════════════════════════\n\n");
            }
            assert(inst.compute_fixed_cost() <= S_star.get_cost());
            ++iter;
        } while (inst.get_active_rows_size() > 0 && S_star.get_cost() - 1.0 > std::ceil(BETA * u_star.get_lb()));

        return std::move(S_star);  // messo il move solo x togliere l'errore poi lo possiam togliere
    }

private:
    Instance& inst;

    SubInstance subinst;
    MStar covered_rows;
    ThreePhase three_phase;

    real_t pi;

    // retain allocated memory (anti-RAII, should be local)
    std::vector<std::pair<idx_t, real_t>> deltas;
    std::vector<idx_t> cols_to_fix;
};

#endif