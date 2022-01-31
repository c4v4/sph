// Copyright (c) 2022 Francesco Cavaliere
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef SPH_INCLUDE_SOLUTION_HPP_
#define SPH_INCLUDE_SOLUTION_HPP_
#include <vector>

#include "SubInstance.hpp"
#include "cft.hpp"

namespace sph {

    /**
     * Solution classes are simply std::vector of indexes.
     * In order to distinguish between solutions belonging different parts of the algorithm, solution of different "type" are actually
     * different classes (in this way we can always know what a solution represent exactly).
     * However, in the end they are simply a vector of indexes plus some utility functions attached to them.
     */
    class LocalSolution;
    class GlobalSolution;

    class BaseSolution : public std::vector<idx_t> {
    public:
        template <typename... Args>
        explicit BaseSolution(Args&&... _args) : std::vector<idx_t>(std::forward<Args>(_args)...) { }

        inline void mark_removal(idx_t idx) { (*this)[idx] = NOT_AN_INDEX; }

        inline void apply_removal() {
            idx_t j2 = 0;
            for (idx_t j1 = 0; j1 < this->size(); ++j1) {
                if ((*this)[j1] != NOT_AN_INDEX) {
                    (*this)[j2++] = (*this)[j1];
                }
            }
            this->resize(j2);
        }

        inline void remove(std::vector<idx_t>& col_idxs) {
            std::sort(col_idxs.begin(), col_idxs.end());  // O(n' log(n'))
            idx_t removed_counter = col_idxs.size();

            idx_t j1 = 0, j2 = 0;
            while (j2 < this->size()) {  // O(n log(n'))
                if (std::binary_search(col_idxs.begin(), col_idxs.end(), (*this)[j2])) {
                    (*this)[j2] = this->back();
                    this->pop_back();
                    if (--removed_counter == 0) {
                        break;
                    }

                } else {
                    (*this)[j1++] = (*this)[j2++];
                }
            }

            while (j2 < this->size()) {
                (*this)[j1++] = (*this)[j2++];
            }  // O(n)
            this->resize(j1);
        }

    protected:
        template <typename Inst>
        [[nodiscard]] inline real_t compute_cost(Inst& inst) const {
            real_t cost = 0.0;
            for (idx_t j : *this) {
                cost += inst.get_col(j).get_cost();
            }
            return cost;
        }
    };

    /**
     * @brief Solution used inside the 3-phase cycle. Use local indexes mappings (defined in the current subinstance) and does not
     * contain fixed columns.
     *
     */
    class LocalSolution : public BaseSolution {
    public:
        template <typename... Args>
        explicit LocalSolution(Args&&... _args) : BaseSolution(std::forward<Args>(_args)...) { }

        [[nodiscard]] inline real_t compute_cost(const SubInstance& subinst) const { return BaseSolution::compute_cost(subinst); }
    };


    /**
     * @brief Complete solution coherent with the complete instance, can be used for I/O.
     *
     */
    class GlobalSolution : public BaseSolution {
    public:
        GlobalSolution() : cost(REAL_MAX) { }
        explicit GlobalSolution(idx_t size) : BaseSolution(size), cost(REAL_MAX) { }
        GlobalSolution(idx_t size, real_t val) : BaseSolution(size, val), cost(REAL_MAX) { }
        GlobalSolution(const GlobalSolution& other) : BaseSolution(other), cost(other.cost) { }
        GlobalSolution(GlobalSolution&& other) noexcept : BaseSolution(std::move(other)), cost(other.cost) { }

        template <typename Container>
        GlobalSolution(Container&& cols, real_t cost_) : BaseSolution(std::forward<Container>(cols)) {
            cost = cost_;
        }

        GlobalSolution& operator=(const GlobalSolution& other) = default;
        GlobalSolution& operator=(GlobalSolution&& other) = default;

        GlobalSolution(const SubInstance& subinst, const LocalSolution& sol) { SetFromLocal(subinst, sol); }

        void SetFromLocal(const SubInstance& subinst, const LocalSolution& sol) {
            clear();
            if (sol.empty()) {
                cost = REAL_MAX;
                return;
            }
            const Instance& inst = subinst.get_instance();
            const UniqueColSet& gcols = inst.get_cols();
            cost = 0.0;
            reserve(sol.size() + subinst.get_fixed_cols().size() + inst.get_fixed_cols().size());
            for (idx_t j : sol) {
                cost += subinst.get_col(j).get_cost();
                emplace_back(subinst.get_global_col_idx(j));
            }
            for (idx_t j : subinst.get_fixed_cols()) {
                cost += gcols[j].get_cost();
                emplace_back(j);
            }
            for (idx_t j : inst.get_fixed_cols()) {
                cost += gcols[j].get_cost();
                emplace_back(j);
            }

            SPH_DEBUG {
                MStar cov_check;
                cov_check.reset_covered(gcols, *this, inst.get_nrows());
                assert(cov_check.get_uncovered() == 0);
            }
        }

        [[nodiscard]] inline real_t get_cost() const { return cost; }

        void set_cost(real_t cost_) { cost = cost_; }

    private:
        real_t cost;
    };

}  // namespace sph

#endif