#ifndef SCP_INCLUDE_CFT_HPP_
#define SCP_INCLUDE_CFT_HPP_

#include <cstdint>

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

#define REMOVED_INDEX (std::numeric_limits<idx_t>::max())

#define HAS_INTEGRAL_COSTS 1.0  // 1.0 if yes , 0.0 if no

#endif  // SCP_INCLUDE_CFT_HPP_
