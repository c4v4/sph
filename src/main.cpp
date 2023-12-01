#include <fmt/core.h>
#include <fmt/ranges.h>

#define VERBOSE_LEVEL 1

#ifndef NDEBUG
#include "SPHeuristic.hpp"
#else
#include "../one_header_only/SPH.hpp"
#endif

#include "parsing.hpp"

static std::pair<sph::GlobalSolution, sph::real_t> run_sph(InstanceData& data, unsigned long seed, double tlim);
void run_cft();
void run_sph();
void run_sph();
void run_sph();
void run_sph();

void print_sol(sph::Instance& inst, sph::GlobalSolution& sol) {
    int route = 1;
    for (sph::idx_t j : sol) {
        fmt::print("Route #{}: {}\n", route++, fmt::join(inst.get_col(j), " "));
    }
    fmt::print("Cost {}\n", sol.get_cost());
}

int main(int argc, char** argv) {

    if (argc < 2) {
        fmt::print("Missing path to instance\n");
        return EXIT_FAILURE;
    }

    const auto path = argv[1];
    auto seed = 0UL;
    double tlim = sph::REAL_MAX;

    if (argc > 2) {
        seed = std::stoul(argv[2]);
    }
    if (argc > 3) {
        tlim = std::stod(argv[3]);
    }

    fmt::print("Command line options:\n");
    fmt::print("Instance path  : {}\n", argv[1]);
    fmt::print("Random seed    : {}\n", seed);
    fmt::print("Time limit (s) : {}\n", tlim);

    InstanceData data = parse_cvrp_instance(path);
    auto [solution, sol_cost] = run_sph(data, seed, tlim);

    fmt::print("Solution (cost {}):\n{}\n", sol_cost, fmt::join(solution, ", "));

    return 0;
}

std::pair<sph::GlobalSolution, sph::real_t> run_sph(InstanceData& data, unsigned long seed, double tlim) {

    sph::SPHeuristic sph(data.nrows, seed);
    std::vector<sph::idx_t> inserted_idxs = sph.add_columns(data.costs, data.solcosts, data.matbeg, data.matval);
    for (sph::idx_t& gj : data.warmstart) {
        gj = inserted_idxs[gj];
    }
    sph.set_timelimit(tlim);
    sph.set_ncols_constr(data.warmstart.size());
    sph.set_max_routes(50000);
    sph.set_keepcol_strategy(sph::SPP);
    sph.set_new_best_callback(print_sol);

    auto solution = sph.solve(data.warmstart);
    sph::real_t sol_cost = 0.0;
    for (auto j : solution) {
        sol_cost += sph.get_col(j).get_cost();
    }

    return {solution, sol_cost};
}