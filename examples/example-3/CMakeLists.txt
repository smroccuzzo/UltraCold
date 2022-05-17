##########################################
# Build instructions for example-3.cpp
##########################################

#------------------------------------------------------
# Always request a reasonably recent version of cmake
#------------------------------------------------------

cmake_minimum_required(VERSION 3.15)

#-----------------------------------------------------
# Search for package UltraCold, which must have been
# installed previously at location ${ULTRACOLD_DIR}
#-----------------------------------------------------

find_package(UltraCold HINTS ${ULTRACOLD_DIR} $ENV{ULTRACOLD_DIR})

if(NOT ${UltraCold_FOUND})
    MESSAGE(FATAL_ERROR
            "\n\n"
            "******************************************************************************\n"
            "Could not locate UltraCold.\n"
            "Pass the flag -DULTRACOLD_DIR=/path/to/UltraCold to cmake"
            " or set an environment variable \"ULTRACOLD_DIR\" that contains this path.\n"
            "******************************************************************************\n"
            "\n\n"
            )
endif()

project(example-3 LANGUAGES CXX)

#--------------
# Find OpenMP
#--------------

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(example-3 example-3.cpp)
ULTRACOLD_SETUP_TARGET(example-3)