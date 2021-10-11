#include "Instance.hpp"

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cassert>
#include <numeric>

#include "CollectionOf.hpp"
#include "MStar.hpp"
#include "Stopwatch.hpp"
#include "cft.hpp"

void Instance::reset_fixing() {

    assert(std::is_sorted(active_cols.begin(), active_cols.end()));
    assert(std::is_sorted(fixed_cols.begin(), fixed_cols.end()));

    active_rows.assign(nrows, true);

    // merge active and fixed columns
    active_cols.resize(cols.size());
    std::iota(active_cols.begin(), active_cols.end(), 0);
    fixed_cols.clear();
    fixed_cost = 0.0;
}

void Instance::fix_columns(const std::vector<idx_t> &idxs, const MStar &M_star) {

    assert(fixed_cols.empty());
    assert(active_rows.size() == nrows);
    assert(M_star.size() == nrows);

    for (idx_t i = 0; i < nrows; ++i) { active_rows[i] = !M_star[i]; }
    nactive_rows = M_star.get_uncovered();

    _fix_columns(idxs);
}


std::vector<idx_t> Instance::add_columns(const Cols &new_cols) {
    idx_t old_ncols = cols.size();
    idx_t ncols = new_cols.size();

    std::vector<idx_t> inserted_cols_idxs;
    inserted_cols_idxs.reserve(new_cols.size());

    for (idx_t j = 0; j < ncols; ++j) {
        cols.emplace_back(new_cols[j]);
        active_cols.emplace_back(old_ncols + j);
        inserted_cols_idxs.emplace_back(old_ncols + j);
    }

    return inserted_cols_idxs;
}

std::vector<idx_t> Instance::add_columns(const std::vector<real_t> &costs, const std::vector<real_t> &sol_costs, const std::vector<idx_t> &matbeg,
                                         const std::vector<idx_t> &matval) {
    assert(costs.size() == sol_costs.size() && costs.size() == matbeg.size() - 1);

    idx_t old_ncols = cols.size();
    idx_t ncols = costs.size();

    std::vector<idx_t> inserted_cols_idxs;
    inserted_cols_idxs.reserve(costs.size());

    for (idx_t j = 0; j < ncols; ++j) {
        cols.emplace_back(matval.data() + matbeg[j], matval.data() + matbeg[j + 1], costs[j], sol_costs[j]);
        active_cols.emplace_back(old_ncols + j);
        inserted_cols_idxs.emplace_back(old_ncols + j);
    }

    assert(cols.size() == costs.size());
    return inserted_cols_idxs;
}


void Instance::_fix_columns(const std::vector<idx_t> &idxs) {

    idx_t iok = 0;
    for (idx_t j = 0; j < cols.size(); ++j) {
        auto &col = cols[j];
        for (idx_t i : col) {
            if (active_rows[i]) {
                active_cols[iok++] = j;
                break;
            }
        }
    }

    active_cols.resize(iok);
    fixed_cols = idxs;

    fixed_cost = 0.0;
    for (idx_t j : fixed_cols) { fixed_cost += cols[j].get_cost(); }
}
