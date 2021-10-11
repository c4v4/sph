#ifndef SPH_INCLUDE_TWOPHASE_HPP_
#define SPH_INCLUDE_TWOPHASE_HPP_

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
    GlobalSolution operator()(const real_t global_UB, GlobalSolution& S_star);
    inline GlobalMultipliers& get_global_u() { return glo_u; }

private:
    void cplex_heur(LocalSolution& S_init, GlobalSolution& S_star, real_t& glb_UB_star);


private:
    SubInstance& subinst;

    SubGradient subgradient;
    ExactSolver exact;

    double exact_time_limit;

    GlobalMultipliers glo_u;
};


#endif