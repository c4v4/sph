#include <fmt/core.h>

#include "Instance.hpp"
#include "Refinement.hpp"
#include "cft.hpp"
#include "parsing.hpp"


int main(int argc, char **argv) {

    if (argc < 2) {
        fmt::print("Missing path to instance\n");
        return EXIT_FAILURE;
    }

    const auto path = argv[1];

    auto seed = 0UL;
    if (argc > 2) { seed = std::stoul(argv[2]); }

    const auto data = parse_cvrp_instance(path);
    // const auto data = parse_rail_instance(path);
    // const auto data = parse_scp_instance(path);

    std::mt19937 rnd(seed);

    auto instance = Instance(data.nrows);
    instance.add_columns(data.costs, data.solcosts, data.matbeg, data.matval);

    Refinement cft(instance, rnd);

    auto solution = cft(data.warmstart);

    real_t sol_cost = 0.0;
    for (auto j : solution) { sol_cost += instance.get_col(j).get_cost(); }
    fmt::print("Solution (cost {}):\n{}\n", sol_cost, fmt::join(solution, ", "));

    MStar coverage;
    coverage.reset_covered(instance.get_cols(), solution, instance.get_nrows());
    fmt::print("Row coverage:\n{}\n", fmt::join(coverage, ", "));


    return 0;
}
