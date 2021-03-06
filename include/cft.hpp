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

#ifndef SPH_INCLUDE_CFT_HPP_
#define SPH_INCLUDE_CFT_HPP_

#include <cstdint>
#include <limits>

/**
 * VERBOSE: If defined, switches on few prints.
 * VERBOSE_LEVEL [0-3]: Define how many prints: 0 few (default); 1 some; 2 a lot; 3 all
 */
#ifdef VERBOSE
#ifndef VERBOSE_LEVEL
#define VERBOSE_LEVEL 0
#endif

#define SPH_VERBOSE(A) if constexpr (VERBOSE_LEVEL >= A)
#else
#define SPH_VERBOSE(A) if constexpr (false)
#endif

#ifdef NDEBUG
#define SPH_DEBUG if constexpr (false)
#else
#define SPH_DEBUG
#endif

namespace sph {

    /* Types and constants */
    typedef uint32_t idx_t;
    typedef double real_t;

    constexpr real_t REAL_MAX = std::numeric_limits<sph::real_t>::max();
    constexpr real_t REAL_LOWEST = std::numeric_limits<sph::real_t>::lowest();
    constexpr idx_t NOT_AN_INDEX = std::numeric_limits<sph::idx_t>::max();
    constexpr idx_t MAX_INDEX = std::numeric_limits<sph::idx_t>::max();
    constexpr real_t HAS_INTEGRAL_COSTS = 1.0;  // 1.0 if yes , 0.0 if no
    constexpr real_t EPSILON = 1E-6;

    enum KeepStrat { SPP, SCP };

    class Instance;
    class GlobalSolution;
    typedef void (*NewBestCallback)(Instance &, GlobalSolution &);


}  // namespace sph

#endif  // SPH_INCLUDE_CFT_HPP_
