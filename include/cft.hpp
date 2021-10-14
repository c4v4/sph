#ifndef SPH_INCLUDE_CFT_HPP_
#define SPH_INCLUDE_CFT_HPP_

#include <cstdint>
#include <limits>

/* Types and constants */
typedef uint32_t idx_t;
typedef double real_t;

#define REAL_MAX (std::numeric_limits<real_t>::max())
#define REAL_LOWEST (std::numeric_limits<real_t>::lowest())
#define REMOVED_INDEX (std::numeric_limits<idx_t>::max())
#define MAX_INDEX (std::numeric_limits<idx_t>::max())
#define HAS_INTEGRAL_COSTS 1.0  // 1.0 if yes , 0.0 if no

#ifndef INST_HARD_CAP
#define INST_HARD_CAP 200'000U
#endif

/* UTILS */
#define NO_INLINE  //__attribute__((noinline))

#ifdef VERBOSE
#define IF_VERBOSE
#else
#define IF_VERBOSE if constexpr (false)
#endif

#ifdef NDEBUG
#define IF_DEBUG if constexpr (false)
#else
#define IF_DEBUG
#endif

#endif  // SPH_INCLUDE_CFT_HPP_
