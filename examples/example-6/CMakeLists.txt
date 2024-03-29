##########################################
# Build instructions for example-6
##########################################

#------------------------------------------------------
# Always request a reasonably recent version of cmake
#------------------------------------------------------

cmake_minimum_required(VERSION 3.18)

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

project(example-6 LANGUAGES CXX CUDA)

#-----------------------------------------------------
# Choose to build either a debug or release version
#-----------------------------------------------------

if(CMAKE_BUILD_TYPE MATCHES Debug)

    message(STATUS "Building example-6 in Debug mode\n")
    add_executable(example-6-debug example-6.cpp)
    ULTRACOLD_SETUP_TARGET_WITH_CUDA(example-6-debug)

elseif(CMAKE_BUILD_TYPE MATCHES Release)

    message(STATUS "Building example-6 in Release mode\n")
    add_executable(example-6 example-6.cpp)
    ULTRACOLD_SETUP_TARGET_WITH_CUDA(example-6)

endif()
