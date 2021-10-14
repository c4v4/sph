#include <fmt/core.h>

#ifndef NDEBUG
#include "SPHeuristic.hpp"
#else
#include "../one_header_only/SPH.hpp"
#endif

#include "parsing.hpp"

int main(int argc, char **argv) {

    if (argc < 2) {
        fmt::print("Missing path to instance\n");
        return EXIT_FAILURE;
    }

    const auto path = argv[1];
    auto seed = 0UL;
    double tlim = REAL_MAX;
    
    if (argc > 2) { seed = std::stoul(argv[2]); }
    if (argc > 3) { tlim = std::stoul(argv[3]); }

    const auto data = parse_cvrp_instance(path);

    SPHeuristic sph(data.nrows, seed);
    sph.add_columns(data.costs, data.solcosts, data.matbeg, data.matval);
    sph.set_timelimit(tlim);

    auto solution = sph.solve(data.warmstart);

    real_t sol_cost = 0.0;
    for (auto j : solution) { sol_cost += sph.get_col(j).get_cost(); }

    fmt::print("Solution (cost {}):\n{}\n", sol_cost, fmt::join(solution, ", "));

    return 0;
}
