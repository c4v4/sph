cmake_minimum_required(VERSION 3.16)
project(SPH)

option(USE_CLANG "build application with clang" ON) # OFF is the default

set(CMAKE_CXX_STANDARD 17)
set(WARNING_FLAGS "-Wall -Wextra -Wpedantic -Wuninitialized -Wshadow ") # -Werror
set(SANITIZERS_FLAGS "-fno-omit-frame-pointer -fsanitize=address -fsanitize-address-use-after-scope -fsanitize=undefined")
#set(OPT_FLAGS "-O3 -march=native -ffat-lto-objects -flto -g -fno-omit-frame-pointer")
set(OPT_FLAGS "-O3 -march=native -flto -fno-exceptions")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${SANITIZERS_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPT_FLAGS}")

set(CMAKE_VERBOSE_MAKEFILE ON)
set(SOURCE src/main.cpp)
set(LIBRARIES fmt pthread dl m)
#set(LIBRARIES pthread dl m)

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

add_executable(SPH ${SOURCE} )
target_link_libraries(SPH PUBLIC ${LIBRARIES})

#ONE HEADER SECTION
add_custom_target(ONE_HEADER ALL
    COMMAND ${CMAKE_COMMAND} 
        -DTARGET_NAME=SPH
        -DTARGET_PATH=${CMAKE_CURRENT_SOURCE_DIR}
        -P ${CMAKE_SOURCE_DIR}/cmake/one_header.cmake
    COMMENT "Creating one header only library..."
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_dependencies(SPH ONE_HEADER)