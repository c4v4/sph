#ifndef SPH_INCLUDE_REFINEMENT_HPP_
#define SPH_INCLUDE_REFINEMENT_HPP_

#include "Instance.hpp"
#include "SubInstance.hpp"
#include "TwoPhase.hpp"
#include "cft.hpp"

#define ALPHA 1.1
#define BETA 0.8
#define PI_MIN 0.3

#define PI_MAX 0.9
#define POST_OPT_TRIALS 100

class Refinement {
public:
    Refinement(Instance& inst_, std::mt19937& rnd_) : inst(inst_), subinst(inst_), two_phase(subinst), rnd(rnd_) { }

    // S_init must be a global complete solution
    std::vector<idx_t> operator()([[maybe_unused]] const std::vector<idx_t>& S_init);
    std::vector<idx_t> inline solve([[maybe_unused]] const std::vector<idx_t>& S_init) { return operator()(S_init); }

private:
    [[nodiscard]] std::vector<idx_t> refinement_fix(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi);
    [[nodiscard]] std::vector<idx_t> random_fix(GlobalSolution S_star, real_t pi);
    [[nodiscard]] std::vector<idx_t> binary_tournament_fix(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi);
    [[nodiscard]] std::vector<idx_t> random_fix2(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi);

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

#endif