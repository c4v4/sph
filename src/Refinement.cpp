#include "Refinement.hpp"

std::vector<idx_t> Refinement::operator()([[maybe_unused]] const std::vector<idx_t>& S_init) {

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

    real_t pi = PI_MIN;
    real_t last_improving_pi = pi;

    int post_optimization_trials = POST_OPT_TRIALS;

    idx_t iter = 1;
    do {
        subinst.reset();  // 2.

        {
            GlobalSolution S = two_phase(S_star.get_cost(), S_star);  // 3. & 4.

            assert(!(std::fabs(pi - PI_MIN) > 0.001 && inst.get_fixed_cols().empty()));
            pi *= ALPHA;  // 6.

            if (S.get_cost() < S_star.get_cost()) {  // update best solution
                S_star = std::move(S);               // 5.

                last_improving_pi = pi;

                // pi = PI_MIN;            // 6.
                pi = std::max(pi / (ALPHA * ALPHA), PI_MIN);  // 6.
            }

            if (iter == 1) { u_star = two_phase.get_global_u(); }

            if (S_star.get_cost() - 1.0 <= BETA * u_star.get_lb() || pi > PI_MAX || inst.get_active_rows_size() <= 0) {

                if (post_optimization_trials <= 0) {
                    IF_VERBOSE {
                        fmt::print("╔═ REFINEMENT: iter {:2} ═════════════════════════════════════════════════════════════\n", iter);
                        fmt::print("║ Early Exit: β(={}) * LB(={}) > UB(={}) - 1\n", BETA, u_star.get_lb(), S_star.get_cost());
                        fmt::print("║ Active rows {}, active cols {}, pi {}\n", inst.get_active_rows_size(), inst.get_active_cols().size(), pi);
                        fmt::print("║ LB {}, UB {}, UB size {}\n", u_star.get_lb(), S_star.get_cost(), S_star.size());
                        fmt::print("╚═══════════════════════════════════════════════════════════════════════════════════\n\n");
                    }
                    break;
                }

                --post_optimization_trials;
                fmt::print("   POST-OPTIMIZATION REFINEMENT: iter {:2}\n", POST_OPT_TRIALS - post_optimization_trials);

                pi = last_improving_pi;
                last_improving_pi = std::max(PI_MIN, last_improving_pi / ALPHA);
            }

        }  //(S and u are potentially been moved, better to encapsulate them into a block)

        if (inst.get_timelimit().exceeded_tlim()) {
            IF_VERBOSE {
                fmt::print("╔═ REFINEMENT: iter {:2} ═════════════════════════════════════════════════════════════\n", iter);
                fmt::print("║ Timelimit exceeded\n");
                fmt::print("╚═══════════════════════════════════════════════════════════════════════════════════\n\n");
            }
            break;
        }

        pi = std::min<real_t>(PI_MAX, pi);

        // 7. Refinement Fix
        inst.reset_fixing();

        cols_to_fix = random_fix(S_star, pi);

        assert(cols_to_fix.size() <= S_star.size());

        std::sort(cols_to_fix.begin(), cols_to_fix.end());
        inst.fix_columns(cols_to_fix, covered_rows);

        IF_VERBOSE {
            fmt::print("╔═ REFINEMENT: iter {:2} ═════════════════════════════════════════════════════════════\n", iter);
            fmt::print("║ Active rows {}, active cols {}, pi {}\n", inst.get_active_rows_size(), inst.get_active_cols().size(), pi);
            fmt::print("║ LB {}, UB {}, UB size {}\n", u_star.get_lb(), S_star.get_cost(), S_star.size());
            fmt::print("╚═══════════════════════════════════════════════════════════════════════════════════\n\n");
        }

        assert(inst.get_fixed_cost() <= S_star.get_cost());
        ++iter;
    } while (true);

    return S_star;
}

std::vector<idx_t> Refinement::refinement_fix(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi) {
    Cols& cols = inst.get_cols();

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

std::vector<idx_t> Refinement::random_fix(GlobalSolution S_star, real_t pi) {
    Cols& cols = inst.get_cols();

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


std::vector<idx_t> Refinement::binary_tournament_fix(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi) {
    Cols& cols = inst.get_cols();

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

[[nodiscard]] inline std::vector<idx_t> Refinement::random_fix2(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi) {
        Cols& cols = inst.get_cols();
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