#ifndef SPH_INCLUDE_INDEXLIST_HPP_
#define SPH_INCLUDE_INDEXLIST_HPP_

#include <algorithm>
#include <vector>

#include "cft.hpp"

namespace sph {

    using IndexList = std::vector<idx_t>;
    using Row = IndexList;

    class Column : public IndexList {

    public:
        template <typename Col>
        Column(const Col& col) : Column(col.begin(), col.end(), col.get_cost(), col.get_solcost()) { }

        template <typename Iter>
        Column(Iter beg_, Iter end_, real_t c_, real_t sol_c_) : IndexList(beg_, end_), c(c_), sol_c(sol_c_) { }

        Column(const Column& other) : IndexList(other.begin(), other.end()), c(other.c), sol_c(other.sol_c) { }
        Column(Column&& other) noexcept
            : IndexList(std::move(other)), c(std::exchange(other.c, 0)), sol_c(std::exchange(other.sol_c, 0)) { }

        Column& operator=(const Column& other) {
            IndexList::operator=(other);
            c = other.c;
            sol_c = other.sol_c;
            return *this;
        };

        Column& operator=(Column&& other) noexcept {
            IndexList::operator=(std::move(other));
            c = std::exchange(other.c, 0);
            sol_c = std::exchange(other.sol_c, 0);
            return *this;
        };

        [[nodiscard]] inline auto get_cost() const { return c; }

        [[nodiscard]] inline auto get_solcost() const { return sol_c; }
        inline void set_solcost(real_t new_c) { sol_c = new_c; }


        template <typename Multipliers>
        [[nodiscard]] inline auto compute_lagr_cost(const Multipliers& u) const {
            real_t local_c_u = c;
            for (idx_t i : *this) {
                assert(i < u.size());
                local_c_u -= u[i];
            }
            return local_c_u;
        }

        bool operator==(const Column& other) const {
            if (c != other.c || sol_c != other.sol_c || size() != other.size()) { return false; }
            for (idx_t n = 0; n < size(); ++n) {
                if ((*this)[n] != other[n]) { return false; }
            }
            return true;
        }

    private:
        real_t c = 0.0;
        real_t sol_c = 0.0;
    };

}  // namespace sph

#endif