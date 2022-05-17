############################################################################
#
# This is the top-level CMakeLists.txt file of the "UltraCold" project.
# It contains instructions to build and install all the libraries 
# provided as part of the project "UltraCold"
#
############################################################################

#--------------------------------------------------
# Set a reasonably recent minimum version of cmake
#--------------------------------------------------

cmake_minimum_required(VERSION 3.15)

#------------------------------------------
# Set the compilers as the Intel's one
#------------------------------------------

set(CMAKE_C_COMPILER icc)
set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_Fortran_COMPILER ifort)

#-------------------------
# Set the c++ standard 
#-------------------------

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#---------------------------------------------------
# Set the name, language and version of the project
#---------------------------------------------------

project(UltraCold LANGUAGES C CXX Fortran VERSION 1.0)

#---------------------------------------
# Set the default build type as Release
#---------------------------------------

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type set to ${CMAKE_BUILD_TYPE}")

#-----------------------------------------------
# If building the Debug version, append _DEBUG 
# to library names
#-----------------------------------------------

set(CMAKE_DEBUG_POSTFIX _DEBUG)

#-------------------------------------------------------
# Add cmake scripts and modules to CMake search path
#-------------------------------------------------------

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#-----------------------------------------------------------------
# Configure, build and install the bundled version of arpack-ng
#-----------------------------------------------------------------

message(STATUS "----------------------------------------------")
message(STATUS "Configuring arpack-ng...")
message(STATUS "----------------------------------------------")
set(ARPACK_DIR "${CMAKE_INSTALL_PREFIX}/bundled/arpack-ng")

execute_process(
        COMMAND ./configure_bundled.sh ${ARPACK_DIR}
        COMMAND_ECHO STDOUT
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bundled
)
message(STATUS "----------------------------------------------")
message(STATUS "arpack-ng built and installed.")
message(STATUS "----------------------------------------------")

#-----------------------------------------------------------------
# Include the directories containing the source files. 
# Each of these directories has its own CMakeLists.txt with 
# instructions to build and install the corresponding libraries
#-----------------------------------------------------------------

add_subdirectory(${CMAKE_SOURCE_DIR}/src/utilities)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/solvers)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/mkl_wrappers)

#---------------------------
# Set useful compiler flags
#---------------------------

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qmkl")

#-----------------------------------
# Include install instructions
#-----------------------------------

include(UltraColdInstall)
