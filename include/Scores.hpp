#ifndef AC_CFT_INCLUDE_SCORE_HPP
#define AC_CFT_INCLUDE_SCORE_HPP

#include <numeric>

#include "Instance.hpp"
#include "MStar.hpp"
#include "Multipliers.hpp"
#include "Solution.hpp"
#include "cft.hpp"


class ScoreData {
public:
    ScoreData() : gamma(std::numeric_limits<real_t>::max()), mu(0), sigma(std::numeric_limits<real_t>::max()), marked(true) { }

    ScoreData(real_t gamma_, idx_t mu_, bool marked_) : gamma(gamma_), mu(mu_), sigma(_score(gamma_, mu_)), marked(marked_) { }

    ScoreData(const std::vector<real_t>& u_k, const MStar& M_star, const Column& j_col) {
        gamma = j_col.get_cost();
        mu = 0UL;
        for (const auto i : j_col) {
            if (M_star[i] == 0) {
                ++mu;
                gamma -= u_k[i];
            }
        }
        sigma = _score(gamma, mu);
        marked = false;
    }

    void update(real_t u_k) {
        gamma += u_k;
        assert(mu > 0);
        --mu;
        sigma = _score(gamma, mu);
    }

    void inserted_in_S() {
        sigma = gamma = std::numeric_limits<real_t>::max();
        mu = 0;
        marked = true;
    }

    [[nodiscard]] bool is_marked() const { return marked; }
    [[nodiscard]] real_t get_score() const { return sigma; }

private:
    static inline real_t _score(real_t gamma, idx_t mu) {
        if (mu <= 0) { return std::numeric_limits<real_t>::max(); }
        if (gamma > 0) { return gamma / static_cast<real_t>(mu); }
        return gamma * static_cast<real_t>(mu);
    }

private:
    real_t gamma;
    idx_t mu;
    real_t sigma;
    bool marked;
};


class Scores {
public:
    explicit Scores(SubInstance& subinst_) : subinst(subinst_), gam_mu_sig() { }


    // Setup for step 2
    inline void setup_scores(const LocalMultipliers& u_k, MStar& M_star) {

        auto& cols = subinst.get_cols();
        Bsize = cols.size();
        B.resize(Bsize);
        std::iota(B.begin(), B.end(), 0);
        Ssize = 0;

        gam_mu_sig.resize(cols.size());
        for (idx_t j = 0; j < gam_mu_sig.size(); ++j) { gam_mu_sig[j] = ScoreData(u_k, M_star, cols[j]); }

        tau = _compute_B_and_tau(std::min<idx_t>(subinst.get_nrows(), cols.size()));
    }

    // Step 1
    inline idx_t argmin_score() {
        idx_t j_star = B[0];
        real_t sigma_jstar = gam_mu_sig[B[0]].get_score();
        for (idx_t j = 1; j < Bsize; ++j) {
            const real_t sigma = gam_mu_sig[B[j]].get_score();
            if (sigma < sigma_jstar) {
                assert(!gam_mu_sig[B[j]].is_marked());
                sigma_jstar = sigma;
                j_star = B[j];
            }
        }

        if (sigma_jstar > tau) {
            tau = _compute_B_and_tau(std::min<idx_t>(subinst.get_nrows(), subinst.get_ncols() - Ssize));
            return argmin_score();
        }

        assert(!gam_mu_sig[j_star].is_marked());
        return j_star;
    }


    // Step 2'
    inline void update_scores(const LocalMultipliers& u_k, const idx_t j_star, const MStar& M_star) {
        assert(!gam_mu_sig[j_star].is_marked());
        gam_mu_sig[j_star].inserted_in_S();
        ++Ssize;

        auto& rows = subinst.get_rows();
        auto& j_star_col = subinst.get_col(j_star);
        for (auto i : j_star_col) {
            if (M_star[i] == 0) {
                for (auto j : rows[i]) {
                    if (!gam_mu_sig[j].is_marked()) { gam_mu_sig[j].update(u_k[i]); }
                }
            }
        }
    }


private:
    inline real_t _compute_B_and_tau(idx_t size) {
        auto last = size - 1;
        assert(last < B.size());
        Bsize = size;
        std::nth_element(B.begin(), B.begin() + last, B.end(),
                         [&](const auto a, const auto b) { return gam_mu_sig[a].get_score() < gam_mu_sig[b].get_score(); });
        return gam_mu_sig[B[last]].get_score();
    }


private:
    SubInstance& subinst;
    std::vector<ScoreData> gam_mu_sig;
    idx_t Ssize;
    idx_t Bsize;
    std::vector<idx_t> B;
    real_t tau;
};


#endif
