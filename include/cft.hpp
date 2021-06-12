#ifndef SCP_INCLUDE_CFT_HPP_
#define SCP_INCLUDE_CFT_HPP_

#include <cstdint>

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

typedef uint32_t idx_t;
typedef float real_t;

#define REAL_MAX (std::numeric_limits<real_t>::max())
#define REAL_LOWEST (std::numeric_limits<real_t>::lowest())

#define REMOVED_INDEX (std::numeric_limits<idx_t>::max())

#define HAS_INTEGRAL_COSTS 1.0  // 1.0 if yes , 0.0 if no

#endif  // SCP_INCLUDE_CFT_HPP_
