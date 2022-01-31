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

#ifndef SPH_INCLUDE_ROUTE_HPP_
#define SPH_INCLUDE_ROUTE_HPP_

#include <algorithm>

#include "IndexList.hpp"
#include "cft.hpp"

namespace sph {

    // https://softwareengineering.stackexchange.com/questions/402542/where-do-magic-hashing-constants-like-0x9e3779b9-and-0x9e3779b1-come-from
    constexpr unsigned long INVERSE_GOLDEN_RATIO = 0x9e3779b97f4a7c15;
    constexpr unsigned long LARGE_PRIME_RATIO = 0x9e3779b97f4a7c55;

    class UniqueCol : public Column {
    public:
        template <unsigned long prime>
        static inline size_t dist_hash_op(size_t seed, idx_t idx) {
            return seed ^ (idx + prime + (seed << 6) + (seed >> 2));
        }

        template <unsigned long prime>
        static inline size_t comb_hash_op(size_t seed, idx_t idx) {
            return seed ^ (idx * prime);
        }

        template <unsigned long prime>
        struct dispHash {
            size_t operator()(const Column& col) const {
                size_t seed = col.size() * col.get_cost();
                for (idx_t idx : col) { seed = dist_hash_op<prime>(seed, idx); }
                return seed;
            }
        };

        template <unsigned long prime>
        struct combHash {
            size_t operator()(const Column& col) const {
                size_t seed = col.size();
                for (idx_t idx : col) { seed = comb_hash_op<prime>(seed, idx); }
                return seed;
            }
        };

        using CombHash1 = combHash<INVERSE_GOLDEN_RATIO>;
        using CombHash2 = combHash<LARGE_PRIME_RATIO>;

    public:
        template <typename... Args>
        UniqueCol(Args&&... _args)
            : Column(std::forward<Args>(_args)...),
              comb_hash1(CombHash1()(*this)),
              comb_hash2(CombHash2()(*this)) { }

        bool operator<(const UniqueCol& other) { return get_cost() < other.get_cost() && elem_wise_equal(other); }

        bool operator>(const UniqueCol& other) { return get_cost() > other.get_cost() && elem_wise_equal(other); }

        bool operator==(const UniqueCol& other) {
            return size() == other.size() && comb_hash1 == other.comb_hash1 && comb_hash2 == other.comb_hash2;
        }

        bool elem_wise_equal(const UniqueCol& other) { return std::equal(begin(), end(), other.begin(), other.end()); }

        [[nodiscard]] inline size_t get_comb_hash1() const { return comb_hash1; }
        [[nodiscard]] inline size_t get_comb_hash2() const { return comb_hash2; }

    private:
        size_t comb_hash1;
        size_t comb_hash2;
    };

}  // namespace sph

#endif