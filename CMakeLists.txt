############################################################################
#
# This is the top-level CMakeLists.txt file of the "UltraCold" project.
# It contains instructions to build and install all the libraries 
# provided as part of the project "UltraCold"
#
############################################################################

#--------------------------------------------------
# Set a reasonably recent minimum version of cmake
# 3.18 is the minimum to work safely with CUDA
#--------------------------------------------------

cmake_minimum_required(VERSION 3.18)

#-----------------------------
# Set the generator as make
#-----------------------------

set (CMAKE_GENERATOR "Unix Makefiles" CACHE INTERNAL "" FORCE)

if(CMAKE_CXX_COMPILER_ID MATCHES Intel)

	#-------------------------------
	# Explicitly set the compilers
	#-------------------------------

	set(CMAKE_CXX_COMPILER icpc)
	set(CMAKE_C_COMPILER icc)
	set(CMAKE_Fortran_COMPILER ifort)

	#---------------------------
	# Set useful compiler flags
	#---------------------------

	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "2021.5.0.20211109")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl")
	else()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qmkl")
	endif()

endif()

#---------------------------------------------
# Decide whether or not to build with CUDA
#---------------------------------------------

option(ULTRACOLD_WITH_CUDA "Build UltraCold with CUDA" OFF)

if(ULTRACOLD_WITH_CUDA)
	message(STATUS "Building UltraCold with cuda libraries")
	set(CMAKE_CUDA_COMPILER nvcc)
endif()

#---------------------------------------------------
# Set the name, language and version of the project
#---------------------------------------------------

project(UltraCold LANGUAGES CUDA C CXX Fortran VERSION 1.0)

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
# Include the directories containing the source files. 
# Each of these directories has its own CMakeLists.txt with 
# instructions to build and install the corresponding libraries
#-----------------------------------------------------------------

add_subdirectory(${CMAKE_SOURCE_DIR}/source/base)

if (ULTRACOLD_WITH_CUDA)
	add_subdirectory(${CMAKE_SOURCE_DIR}/source/cuda)
endif ()

#-----------------------------------
# Include install instructions
#-----------------------------------

include(UltraColdInstall)
