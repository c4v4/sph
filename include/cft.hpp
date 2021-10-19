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

#ifndef INST_HARD_CAP
    constexpr unsigned INST_HARD_CAP = 200'000U;
#endif

}  // namespace sph

#endif  // SPH_INCLUDE_CFT_HPP_
