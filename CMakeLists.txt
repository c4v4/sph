cmake_minimum_required(VERSION 3.16)
project(SPH)

set(CMAKE_CXX_STANDARD 17)

set(WARNING_FLAGS "-Wall -Wextra -Wpedantic -Wuninitialized -Wshadow ") # -Werror
set(SANITIZERS_FLAGS "-fno-omit-frame-pointer -fsanitize=address -fsanitize-address-use-after-scope -fsanitize=undefined")
#set(OPT_FLAGS "-O3 -march=native -ffat-lto-objects -flto -g -fno-omit-frame-pointer")
set(OPT_FLAGS "-O3 -march=native -ffat-lto-objects -flto")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${SANITIZERS_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPT_FLAGS}")

set(CMAKE_VERBOSE_MAKEFILE ON)

set(SOURCE src/main.cpp)

set(LIBRARIES fmt pthread dl m)

# Find CPLEX
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(CPLEX)
if (CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIRS})
    message(${CPLEX_INCLUDE_DIRS})
    set(LIBRARIES ${LIBRARIES} ${CPLEX_LIBRARIES} pthread dl m)
    add_definitions(-DIL_STD)
endif(CPLEX_FOUND)

option(ENABLE_VERBOSE "Enable verbose output" OFF)

message("-- Build options")

if(ENABLE_VERBOSE)
    message("--- Verbose output ENABLED")
    add_definitions(-DVERBOSE)
else()
    message("--- Verbose output DISABLED")
endif()


include_directories(include)

add_executable(SPH ${SOURCE})
target_link_libraries(SPH PUBLIC ${LIBRARIES})


set(SHPHEADERS  include/cft.hpp
                include/CollectionOf.hpp
                include/IndexList.hpp
                include/MStar.hpp
                include/Stopwatch.hpp
                include/Instance.hpp
                include/TrivialHeap.hpp
                include/SubInstance.hpp
                include/Solution.hpp
                include/ExactSolver.hpp
                include/Multipliers.hpp
                include/LowerBound.hpp
                include/Counter.hpp
                include/SubGradientUtils.hpp
                include/SubGradient.hpp
                include/TwoPhase.hpp
                include/Refinement.hpp
                include/SPHeuristic.hpp)

function(cat IN_FILE OUT_FILE)
    file(READ ${IN_FILE} CONTENTS)
    file(APPEND one_header_only/SPH.hpp "\n\n/* ################################################################# */\n")
    file(APPEND one_header_only/SPH.hpp "/* #### Original Header: ${SPHHEADER} */\n")
    file(APPEND one_header_only/SPH.hpp "/* ################################################################# */\n\n")
    file(APPEND ${OUT_FILE} "${CONTENTS}")
endfunction()

file(WRITE one_header_only/SPH.hpp "/* Automatically generated one-header-only library */\n")

foreach(SPHHEADER ${SHPHEADERS})
  cat(${SPHHEADER} one_header_only/SPH.hpp)
endforeach()


function(deleteinplace IN_FILE pattern)
  file (STRINGS ${IN_FILE} LINES ENCODING "UTF-8")
  file(WRITE ${IN_FILE} "")

  foreach(LINE IN LISTS LINES)
    string(REGEX REPLACE ${pattern} "/* ${LINE} */" STRIPPED "${LINE}")
    file(APPEND ${IN_FILE} "${STRIPPED}\n")
  endforeach()
endfunction()

deleteinplace(one_header_only/SPH.hpp "\#include \".*")

