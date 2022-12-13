#--------------------------------------------------------
# This is the Config file for the UltraCold package
# It allows to use UltraCold in external CMake projects
#--------------------------------------------------------

#-------------------------------------
# Where to search for modules
#-------------------------------------

list(APPEND CMAKE_MODULE_PATH "${ULTRACOLD_DIR}/cmake")
include(UltraColdTargets)

#-------------------------------------------
# Macro for working with the base version
#-------------------------------------------

macro(ULTRACOLD_SETUP_TARGET target)

    #-------------------------------------
    # Set the compilers as the Intel ones
    #-------------------------------------

    set(CMAKE_CXX_COMPILER icpc)
    set(CMAKE_C_COMPILER icc)

    #----------------------------------------------
    # Set the correct flags for linking with mkl
    #----------------------------------------------

    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "2021.5.0.20211109")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qmkl")
    endif()

    #----------------------------------------
    # Set the fundamental include directory
    #----------------------------------------

    target_include_directories(${target} PUBLIC "${ULTRACOLD_DIR}/include")

    #-------------------------------
    # Setup and link to arpack-ng
    #-------------------------------

    set(ARPACK_DIR "${ULTRACOLD_DIR}/base/bundled/arpack-ng")
    find_package(arpack-ng REQUIRED HINTS ${ARPACK_DIR} $ENV{ARPACK_DIR})
    target_include_directories(${target} PUBLIC ${ARPACK_DIR}/include/arpack)
    target_link_directories(${target} PUBLIC ${ARPACK_DIR}/lib)
    target_link_directories(${target} PUBLIC ${ARPACK_DIR}/lib64)

    #----------------------------------------
    # Link with the remaining libraries
    #----------------------------------------

    target_link_libraries(${target} PUBLIC utilities)
    target_link_libraries(${target} PUBLIC solvers)
    target_link_libraries(${target} PUBLIC ARPACK::ARPACK)

endmacro()

#----------------------------------------
# Working with the cuda version
#----------------------------------------

macro(ULTRACOLD_SETUP_TARGET_WITH_CUDA target)

    #--------------------------------------------
    # Still need a c++ compiler. Set it still as
    # as the Intel one, plus link to MKL
    #--------------------------------------------

    set(CMAKE_CXX_COMPILER icpc)
    set(CMAKE_C_COMPILER icc)

    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "2021.5.0.20211109")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qmkl")
    endif()

    #----------------------------------------
    # Set the fundamental include directory
    #----------------------------------------

    target_include_directories(${target} PUBLIC "${ULTRACOLD_DIR}/include")

    #----------------------------------------------------------
    # Now find and link to CUDA. It is also fundamental to let
    # CUDA resolve its own symbols.
    #----------------------------------------------------------

    find_package(CUDA 11.0 REQUIRED)
    enable_language(CUDA)
    set(CMAKE_CUDA_COMPILER nvcc)
    set_target_properties(${target} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    target_include_directories(${target} PUBLIC "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
    target_link_libraries(${target} PUBLIC "${CUDA_LIBRARIES}" -lcufft)
    add_compile_definitions(ULTRACOLD_WITH_CUDA)

    #-------------------------
    # Link to other libraries
    #-------------------------

    target_link_libraries(${target} PUBLIC cudaSolvers)
    target_link_libraries(${target} PUBLIC utilities)

endmacro()
