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
endfunction()

one_header_only()