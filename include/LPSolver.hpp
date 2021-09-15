#ifndef SCP_INCLUDE_LPSOLVER_HPP_
#define SCP_INCLUDE_LPSOLVER_HPP_

#include <ilcplex/ilocplex.h>

#include "Instance.hpp"
#include "Multipliers.hpp"
#include "Solution.hpp"
#include "SubGradientUtils.hpp"
#include "cft.hpp"


class LPSolver {
public:
    LPSolver(SubInstance& subinst_) : subinst(subinst_) { env = CPXopenCPLEX(NULL); };
    ~LPSolver() { CPXcloseCPLEX(&env); }

    LocalMultipliers solve(real_t UB, const LocalMultipliers& u_0) {

        const idx_t nrows = subinst.get_rows().size();

        // PARAMETERS
        idx_t max_iter = nrows * 0.1;
        auto time_to_exit = ExitCondition(300U);

        u = u_0;
        u_star = u_0;

        MStar covered_rows;
        LocalSolution S;
        real_t real_LB = LB_star = lagr_mul_LB(subinst, u);

        fmt::print("Initial LB: {}\n", real_LB);

        std::vector<std::pair<idx_t, real_t>> delta_u;

        for (idx_t iter = 0; iter < max_iter; ++iter) {
            if (subinst.get_timelimit().exceeded_tlim()) { break; }

            S.clear();
            real_LB = lp_solve(u, S);

            fmt::print("LP obj: {}, mult LB: {}\n", real_LB, lagr_mul_LB(subinst, u));

            if (real_LB > LB_star) {
                LB_star = real_LB;
                u_star = u;
            }

            fmt::print("[{:^4}] Lower Bound: {:.4} (best {:.4}), S_cost {:.4}\n", iter, real_LB, LB_star, S.compute_cost(subinst));

            if (real_LB > UB - HAS_INTEGRAL_COSTS) {
                IF_VERBOSE { fmt::print(" WARNING: real_LB({}) > UB({}) - {}\n", real_LB, UB, HAS_INTEGRAL_COSTS); }
                u_star = u;
                return u_star;
            }

            if (!S.empty()) {
                covered_rows.reset_covered(subinst.get_cols(), S, nrows);
                if (covered_rows.get_uncovered() == 0) {
                    real_t S_cost = S.compute_cost(subinst);
                    if (S_cost < UB) {
                        fmt::print("IMPROVED UB {} --> {}\n", UB, S_cost);
                        UB = S_cost;
                    }
                }
            }


            if (time_to_exit(LB_star)) {
                IF_VERBOSE { fmt::print("[{:^4}] Lower Bound: {:.4} (best {:.4}), S_cost {:.4}\n", iter, real_LB, LB_star, S.compute_cost(subinst)); }
                return u_star;
            }

            const auto global_LB = subinst.price(u);
            fmt::print("Pricing: global {}, local {}\n", global_LB, real_LB);
            LB_star = real_LB = lagr_mul_LB(subinst, u);
        }

        return u_star;
    }

    real_t lp_solve(LocalMultipliers& u_, LocalSolution& S_) {
        int status = 0;
        CPXLPptr lp = CPXcreateprob(env, &status, "lp");

        CPXchgobjsen(env, lp, CPX_MIN);
        CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_OFF);
        // CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 2);

        double zero = 0.0;
        double one = 1.0;
        char binary = 'B';
        double obj;
        for (SubInstCol& col : subinst.get_cols()) {
            obj = col.get_cost();
            status = CPXnewcols(env, lp, 1, &obj, &zero, &one, &binary, nullptr);
        }

        int zero_int = 0;
        char sense = 'E';
        std::vector<double> ones(100, 1.0);
        std::vector<CPXDIM> int_copy(100);
        for (Row& row : subinst.get_rows()) {
            ones.assign(row.size(), 1.0);
            int_copy.resize(row.size());
            for (idx_t n = 0; n < row.size(); ++n) { int_copy[n] = row[n]; }

            status = CPXaddrows(env, lp, 0, 1, row.size(), &one, &sense, &zero_int, int_copy.data(), ones.data(), nullptr, nullptr);
        }

#ifndef NDEBUG
        CPXwriteprob(env, lp, "lp_model.lp", nullptr);
#endif

        CPXchgprobtype(env, lp, CPXPROB_LP);
        status = CPXdualopt(env, lp);
        if (status) { fmt::print("Cplex dual simplex finished with error {}\n", status); }

        int stat = 0;
        obj = 0.0;
        int ncols = subinst.get_ncols();

        vals.resize(ncols + 1);
        cplex_u.resize(subinst.get_nrows());

        status = CPXsolution(env, lp, &stat, &obj, vals.data(), cplex_u.data(), nullptr, nullptr);

        u_ = LocalMultipliers(cplex_u.begin(), cplex_u.end());  // double -> real_t conversion

        if (!status && is_integer(vals)) {

            for (int i = 0; i < ncols; ++i) {
                if (vals[i] > 0.5) { S_.emplace_back(i); }
            }

            real_t sol_cost = S_.compute_cost(subinst);
            fmt::print("Solution (cost {}):\n{}\n", sol_cost, fmt::join(S_, ", "));
        }

        CPXfreeprob(env, &lp);

        return obj;
    }

    double is_integer(std::vector<double>& vals_) {
        for (idx_t i = 0; i < subinst.get_ncols(); ++i) {
            if (vals_[i] > 1E-4 && vals_[i] < 1.0 - 1E-4) { return false; }
        }
        return true;
    }

private:
    SubInstance& subinst;
    CPXENVptr env;

    std::vector<double> vals;

    real_t LB_star;
    LocalMultipliers u_star;
    LocalMultipliers u;
    std::vector<double> cplex_u;
};


#endif