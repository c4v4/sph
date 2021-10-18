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


    function(deleteinplace IN_FILE)
        message("Reading file ${IN_FILE}")
        file (STRINGS ${IN_FILE} LINES_RAW NEWLINE_CONSUME ENCODING "UTF-8")
        
        # CMAKE list are separated using ";". So this MF put ; everywhere and escape all the ";" 
        # originally in the text.
        # But here's the catch: if you terminate you line with \, it produces \; that is transformed into ; 
        # when the list of string is red line by line, thus removing a new-line and spwaning a ";" 
        # out of oblivion. 
        string(REGEX REPLACE "\\\\;" ";" LINES "${LINES_RAW}") 
        string(REGEX REPLACE "(\#include \"[A-Za-z0-9_\\.]*\")" "/* \\1 */" STRIPPED "${LINES}")
        file(WRITE ${IN_FILE} "${STRIPPED}")

    endfunction()

    deleteinplace(one_header_only/SPH.hpp)
endfunction()

one_header_only()
file(COPY include/fmt DESTINATION one_header_only/)