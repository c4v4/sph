#ifndef SPH_INCLUDE_EXACTSOLVER_HPP_
#define SPH_INCLUDE_EXACTSOLVER_HPP_

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <ilcplex/cplex.h>

#include "Solution.hpp"
#include "SubInstance.hpp"
#include "cft.hpp"

#define RESIZE_UP(vec, sz) \
    if (vec.size() < sz) { vec.resize(sz); }

#define ASSIGN_UP(vec, sz, val) \
    if (vec.size() < sz) { vec.assign(sz, val); }

#define SET_INT(P, VAL)                                                            \
    if (int res = 0; (res = CPXsetintparam(env, P, VAL))) {                        \
        fmt::print(stderr, "Error while setting " #P " parameter at  {} \n", VAL); \
        return res;                                                                \
    }

#define SET_DBL(P, VAL)                                                           \
    if (int res = 0; (res = CPXsetdblparam(env, P, VAL))) {                       \
        fmt::print(stderr, "Error while setting " #P " parameter at {} \n", VAL); \
        return res;                                                               \
    }

class ExactSolver {
public:
    ExactSolver() : env(CPXopenCPLEX(nullptr)) { }

    ~ExactSolver() { CPXcloseCPLEX(&env); }

    LocalSolution build_and_opt(SubInstance& subinst, LocalSolution& warmstart, Timer& time_limit) {

        lp = CPXcreateprob(env, nullptr, "exact");
        int res = 0;

        if ((res = build_model(subinst))) {
            fmt::print(stderr, "Error while building the model (errno: {})\n", res);
            return LocalSolution();
        }

        if ((res = set_warmstart(warmstart))) { fmt::print(stderr, "Error while setting warmstart (errno: {})\n", res); }

        if ((res = set_CPX_params(time_limit.seconds_until_end()))) {
            fmt::print(stderr, "Error while setting parameter (errno: {})\n", res);
            return LocalSolution();
        }

        IF_DEBUG {
            if ((res = CPXwriteprob(env, lp, "model.lp", nullptr))) { fmt::print(stderr, "Error while writing problem file(errno: {})\n", res); }
        }

        if ((res = CPXmipopt(env, lp))) {
            fmt::print("Cplex finished with error {}\n", res);
            return LocalSolution();
        }

        int stat = 0;
        double obj = 0.0;
        int ncols = subinst.get_ncols();
        RESIZE_UP(dbl_vals, ncols + 1U);
        LocalSolution sol;

        if (((res = CPXsolution(env, lp, &stat, &obj, dbl_vals.data(), nullptr, nullptr, nullptr) == 0))) {
            for (int i = 0; i < ncols; ++i) {
                if (dbl_vals[i] > 0.5) { sol.emplace_back(i); }
            }

            MStar coverage(subinst.get_nrows());
            coverage.reset_covered(subinst.get_cols(), sol, subinst.get_nrows());
            if (coverage.get_uncovered() > 0) {
                fmt::print(stderr, "Error, solution does not cover all the rows!\n");
                fmt::print(stderr, " Row coverage:\n{}\n", fmt::join(coverage, ", "));
                fmt::print(stderr, " Cols value: {}\n", fmt::join(dbl_vals, ", "));
                return LocalSolution();
            }
        }

        CPXfreeprob(env, &lp);

        return sol;
    }

private:
    int build_model(SubInstance& subinst) {
        int res;

        if ((res = CPXchgobjoffset(env, lp, subinst.get_fixed_cost()))) {
            fmt::print(stderr, "Error while setting obj func constant term! (errno: {})\n", res);
            return res;
        }

        if ((res = add_variables(subinst))) {
            fmt::print(stderr, "Error while creating columns! (errno: {})\n", res);
            return res;
        }

        if ((res = add_constraints(subinst))) {
            fmt::print(stderr, "Error creating rows! (errno: {})\n", res);
            return res;
        }

        return 0;
    }

    int add_variables(SubInstance& subinst) {
        SubInstCols& cols = subinst.get_cols();
        idx_t ncols = subinst.get_ncols();

        ASSIGN_UP(ctype, ncols, 'B');
        ASSIGN_UP(ones, ncols, 1.0);
        ASSIGN_UP(lb, ncols, 0.0);
        RESIZE_UP(dbl_vals, ncols);

        std::transform(cols.begin(), cols.end(), dbl_vals.begin(), [](auto& c) { return c.get_cost(); });

        return CPXnewcols(env, lp, ncols, dbl_vals.data(), lb.data(), ones.data(), ctype.data(), nullptr);
    }

    int add_constraints(SubInstance& subinst) {
        std::vector<Row>& rows = subinst.get_rows();
        idx_t nrows = subinst.get_nrows();

        RESIZE_UP(rmatbeg, nrows);
        rmatind.clear();
        rmatind.reserve(nrows * 5U);

        int nzcount = 0;
        for (idx_t i = 0; i < nrows; ++i) {
            rmatbeg[i] = nzcount;
            nzcount += rows[i].size();
            rmatind.insert(rmatind.end(), rows[i].begin(), rows[i].end());
        }

        ASSIGN_UP(ones, static_cast<size_t>(nzcount), 1.0);
        ASSIGN_UP(sense, static_cast<size_t>(nrows), 'E');

        return CPXaddrows(env, lp, 0, nrows, nzcount, ones.data(), sense.data(), rmatbeg.data(), rmatind.data(), ones.data(), nullptr, nullptr);
    }

    int set_warmstart(LocalSolution& warmstart) {
        idx_t wsize = warmstart.size();
        if (wsize > 0) {

            int zero_int = 0;
            int effort = CPX_MIPSTART_NOCHECK;
            ASSIGN_UP(ones, wsize, 1.0);
            RESIZE_UP(rmatind, wsize);
            for (idx_t n = 0; n < wsize; ++n) { rmatind[n] = warmstart[n]; }

            return CPXaddmipstarts(env, lp, 1, wsize, &zero_int, rmatind.data(), ones.data(), &effort, nullptr);
        }

        return 0;
    }

    int set_CPX_params(double tlim) {

        if (int res = 0; (res = CPXchgobjsen(env, lp, CPX_MIN))) {
            fmt::print(stderr, "Error while setting objective function ctype! (errno: {})\n", res);
            return res;
        }

        SET_INT(CPXPARAM_ScreenOutput, CPX_OFF);
        SET_INT(CPXPARAM_MIP_Display, 2);
        SET_INT(CPXPARAM_Threads, 1);
        SET_INT(CPXPARAM_Parallel, CPX_PARALLEL_OPPORTUNISTIC);
        SET_DBL(CPXPARAM_MIP_PolishAfter_Time, tlim * 0.5);
        SET_DBL(CPXPARAM_TimeLimit, tlim);
        SET_INT(CPXPARAM_Emphasis_MIP, CPX_MIPEMPHASIS_OPTIMALITY);

        IF_VERBOSE {
            IF_DEBUG { SET_INT(CPXPARAM_ScreenOutput, CPX_ON); }
        }

        return 0;
    }

    ////////// PRIVATE FIELDS //////////

    CPXENVptr env;
    CPXLPptr lp;

    // Avoid repeated memory allocations
    std::vector<double> lb;        // variables lower bounds (ncols, 0.0)
    std::vector<double> ones;      // constaints coefficients (nzcount, 1.0)
    std::vector<char> ctype;       // binary variables (ncols, 'B')
    std::vector<int> rmatind;      // constraints indices (nzcount)
    std::vector<int> rmatbeg;      // constraints begin in rmatind or rmatval (nrows)
    std::vector<char> sense;       // constraints sense (nrows, 'E')
    std::vector<double> dbl_vals;  // store variables final values (ncols)
};


#endif