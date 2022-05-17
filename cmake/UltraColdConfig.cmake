#--------------------------------------------------------
# This is the Config file for the UltraCold package
# It allows to use UltraCold in external CMake projects
#--------------------------------------------------------

# Set the compiler

set(CMAKE_C_COMPILER icc)
set(CMAKE_CXX_COMPILER icpc)

# Where to search for modules

list(APPEND CMAKE_MODULE_PATH "${ULTRACOLD_DIR}/cmake")
include(UltraColdTargets)

# Now the setup macros

macro(ULTRACOLD_SETUP_TARGET target)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qmkl")

    target_include_directories(${target} PUBLIC "${ULTRACOLD_DIR}/include")

    set(ARPACK_DIR "${ULTRACOLD_DIR}/bundled/arpack-ng")
    find_package(arpack-ng REQUIRED HINTS ${ARPACK_DIR} $ENV{ARPACK_DIR})
    target_include_directories(${target} PUBLIC ${ARPACK_DIR}/include/arpack)
    target_link_directories(${target} PUBLIC ${ARPACK_DIR}/lib)
    target_link_directories(${target} PUBLIC ${ARPACK_DIR}/lib64)

    target_link_libraries(${target} PUBLIC mkl_wrappers)
    target_link_libraries(${target} PUBLIC utilities)
    target_link_libraries(${target} PUBLIC solvers)
    target_link_libraries(${target} PUBLIC ARPACK::ARPACK)

endmacro()
