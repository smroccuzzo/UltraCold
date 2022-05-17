###################################################################
# Build a library of classes wrapping useful Intel's MKL functions
###################################################################

#-------------------
# Find OpenMP
#-------------------

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#---------------------------------------
# Create the library of wrapper classes
#---------------------------------------

add_library(mkl_wrappers
    STATIC
    ${CMAKE_SOURCE_DIR}/src/mkl_wrappers/dft.cpp
)

#-------------------------------
# Set the include directories
#-------------------------------

target_include_directories(mkl_wrappers 
    INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

#----------------
# Link libraries 
#----------------

target_link_libraries(mkl_wrappers PUBLIC utilities)

#-------------------------------------------------------
# Install library and headers to their final location
#-------------------------------------------------------

install(TARGETS mkl_wrappers 
        DESTINATION lib
        EXPORT UltraColdTargets
)

install(FILES dft.hpp DESTINATION include)