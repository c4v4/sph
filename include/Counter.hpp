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

#ifndef SPH_INCLUDE_COUNTER_HPP_
#define SPH_INCLUDE_COUNTER_HPP_

#include "cft.hpp"

namespace sph {

    class Counter {

    public:
        explicit Counter(idx_t max) : m(max), i(0) { }
        inline void inc() { i++; }
        [[nodiscard]] inline auto reached() const { return i >= m; }
        inline void restart() { i = 0; }
        inline void set_max(idx_t new_max) { m = new_max; }
        [[nodiscard]] inline idx_t get_max() const { return m; }
        [[nodiscard]] inline auto get() const { return i; }

    private:
        idx_t m;  // max value
        idx_t i;  // current value
    };
}  // namespace sph

#endif  // SPH_INCLUDE_COUNTER_HPP_
