#ifndef SPH_INCLUDE_SPHEURISTIC_HPP_
#define SPH_INCLUDE_SPHEURISTIC_HPP_

#include <cassert>

#include "Instance.hpp"
#include "Refinement.hpp"

namespace sph {

    class SPHeuristic {
    public:
        explicit SPHeuristic(const idx_t nrows_) : inst(nrows_), rnd(std::mt19937()), refinement(inst, rnd) { }
        SPHeuristic(const idx_t nrows_, int seed) : inst(nrows_), rnd(std::mt19937(seed)), refinement(inst, rnd) { }

        [[nodiscard]] inline idx_t get_ncols() const { return inst.get_ncols(); }
        [[nodiscard]] inline idx_t get_nrows() const { return inst.get_nrows(); }
        [[nodiscard]] inline UniqueColSet &get_cols() { return inst.get_cols(); }
        [[nodiscard]] inline Column &get_col(idx_t idx) { return inst.get_col(idx); }
        [[nodiscard]] inline const Column &get_col(idx_t idx) const { return inst.get_col(idx); }

        /**
         * @brief Set the ncols constr constraint rhs
         *
         * @param ncols_constr
         */
        inline void set_ncols_constr(idx_t ncols_constr) { inst.set_ncols_constr(ncols_constr); }

        /**
         * @brief Set the timelimit to <now> + <seconds>.
         *          Thus, to run the algorithm for, e.g., 10 seconds,
         *          call this function right before the solve method.
         *
         * @param seconds At how many second the timelimit is set.
         */
        void inline set_timelimit(double seconds) { inst.set_timelimit(seconds); }

        /**
         * @brief Get the current timelimit set
         *
         * @return Timer&
         */
        [[nodiscard]] inline Timer &get_timelimit() { return inst.get_timelimit(); }

        /**
         * @brief Call add_column for each one of the column in new_cols.
         *
         * @tparam UniqueColContainer Container type.
         * @param new_cols New columns container.
         * @return std::vector<idx_t> Vector of length costs.size(), mapping
         *          each original column index to the index inside Instance.
         */
        template <typename UniqueColContainer>
        std::vector<idx_t> inline add_columns(const UniqueColContainer &new_cols) {
            return inst.add_columns(new_cols);
        }

        /**
         * @brief Call add_column for each one of the column represented in a CPLEX-like way.
         *
         * @param costs vector of costs one for each column
         * @param sol_costs vector of solution cost of each column
         * @param matbeg vector containing at index "i" the starting
         *              position of column "i" rows inside matval
         *              vector. Thus, column "i" will start from
         *              matbeg[i] and end at matbeg[i+1]-1 (or at
         *              matval.size() for the last column).
         * @param matval vector containing all the columns in a contigous
         *              representation.
         * @return std::vector<idx_t> Vector of length costs.size(), mapping
         *          each original column index to the index inside Instance.
         */
        std::vector<idx_t> inline add_columns(const std::vector<real_t> &costs, const std::vector<real_t> &sol_costs,
                                              const std::vector<idx_t> &matbeg, const std::vector<idx_t> &matval) {
            return inst.add_columns(costs, sol_costs, matbeg, matval);
        }

        /**
         * @brief Tries to add a column into the current instance.
         *      A column is added using an heuristic criterion:
         *      for each column, two hash functions based only on
         *      the rows of the column (and not on their order) are
         *      computed. For every pair of such hashes only one
         *      column is maintained inside the column set.
         *      When a collision happens, the new column replace the
         *      old one if at least one of these conditions holds:
         *      1. Its the same column but with a better cost.
         *      2. It belongs to a better solution.
         *
         * @tparam _Args args types
         * @param args Argument to construct the column inplace.
         * @return idx_t The index of the colums at the column
         *              inserted, or that blocked the insertion.
         */
        template <typename... _Args>
        idx_t add_column(_Args &&...args) {
            return inst.add_column(std::forward<_Args>(args)...);
        }

        /**
         * @brief Call the sph main heuristic algorithm which tries to find the best
         * solution it can in the timelimit set.
         *
         * @tparam ROUTES_HARD_CAP Maximum instance size. After the first iteration is
         *          done, the lagrangian multiplier are used to filter out bad columns.
         *          Three differnt criteria are used, for each one ROUTES_HARD_CAP set
         *          the maximum amount of columns to select. Thus, in the worst case
         *          scenario, 3 * ROUTES_HARD_CAP are selected. Usually, only
         *          ~1.1 * ROUTES_HARD_CAP are selected, since the 3 criteria overlap.
         *
         * @tparam KeepColStrategy Choose betwee two fixing options:
         *          1. SetPar_ActiveColTest: when an active column is fixed, all the
         *              other active columns which overlap with at least one row, are
         *              removed from the subinstace.
         *          2. SetCov_ActiveColTest: when fixing a set of active columns,
         *              only the ones that have at least one active row are maintained.
         *              (to avoid empty columns inside the subinstance)
         *
         * @param S_init Initial solution for the algorithm.
         *
         * @return std::vector<idx_t> Best solution found.
         */
        template <unsigned long ROUTES_HARD_CAP = INST_HARD_CAP, typename KeepColStrategy = SetPar_ActiveColTest>
        std::vector<idx_t> inline solve([[maybe_unused]] const std::vector<idx_t> &S_init) {
            SPH_VERBOSE(0) { fmt::print(" Set Partitioning Heuristic: \n"); }
            SPH_DEBUG { fmt::print(" Warning: running in debug mode\n"); }
            SPH_VERBOSE(0) { fmt::print(" SP Instance size: {}x{}\n", inst.get_nrows(), inst.get_ncols()); }

            return refinement.solve<ROUTES_HARD_CAP, KeepColStrategy>(S_init);
        }

    private:
        Instance inst;
        std::mt19937 rnd;
        Refinement refinement;
    };

}  // namespace sph

#endif