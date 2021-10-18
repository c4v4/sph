function(one_header_only)
    set(SHPHEADERS  include/cft.hpp
                    include/CollectionOf.hpp
                    include/MStar.hpp
                    include/Stopwatch.hpp
                    include/IndexList.hpp
                    include/VectorSet.hpp
                    include/UniqueCol.hpp
                    include/UniqueColSet.hpp
                    include/Instance.hpp
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
        string(LENGTH "${SPHHEADER}" NAME_LENGTH)
        MATH(EXPR COUNT "37 - ${NAME_LENGTH}")
        string(REPEAT " " ${COUNT} FILLER)
        file(APPEND one_header_only/SPH.hpp "\n\n/* ######################################################################### */\n")
        file(APPEND one_header_only/SPH.hpp     "/* ######## Original Header: ${SPHHEADER} ${FILLER} ######## */\n")
        file(APPEND one_header_only/SPH.hpp     "/* ######################################################################### */\n\n")
        file(APPEND ${OUT_FILE} "${CONTENTS}")
    endfunction()

    file(WRITE one_header_only/SPH.hpp  "/* Automatically generated one-header-only library */\n\n")
    file(APPEND one_header_only/SPH.hpp "
/* 
 * Copyright (C) 2021 Francesco Cavaliere - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GPL-3 license.
 *
 * You should have received a copy of the GPL-3 license with
 * this file. If not, please write to: f.cavaliere@unibo.it, 
 * or try visit: https://github.com/c4v4/sph 
 */\n")
    file(APPEND one_header_only/SPH.hpp "
/* 
 * A subset of the of the excelent fmtlib is included in this 
 * library in header-only mode.
 * Include you local version before this file if you want to 
 * use it (can decrease compile-time).
 */\n")

    foreach(SPHHEADER ${SHPHEADERS})
    cat(${SPHHEADER} one_header_only/SPH.hpp)
    endforeach()


    function(deleteinplace IN_FILE pattern)
        message("Reading file ${IN_FILE}")
        file (STRINGS ${IN_FILE} LINES ENCODING "UTF-8")
        file(WRITE ${IN_FILE} "")

        foreach(LINE IN LISTS LINES)
            string(REGEX REPLACE ${pattern} "/* ${LINE} */" STRIPPED "${LINE}")
            string(REGEX REPLACE ".*\#include \"fmt/.*" "${LINE}" RESTORED "${STRIPPED}")
            file(APPEND ${IN_FILE} "${RESTORED}\n")
        endforeach()
    endfunction()

    deleteinplace(one_header_only/SPH.hpp "\#include \".*")
endfunction()

one_header_only()
file(COPY include/fmt DESTINATION one_header_only/)