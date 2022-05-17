##################################################################################
# Build a library containing all solver classes based on Intel's MKL 
##################################################################################

#--------------
# Find OpenMP
#--------------

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#-------------------
# Find arpack-ng
#-------------------

find_package(arpack-ng REQUIRED HINTS ${ARPACK_DIR} $ENV{ARPACK_DIR})

#-------------------------------
# Create the solvers library
#-------------------------------

add_library(solvers
    STATIC
    ${CMAKE_SOURCE_DIR}/src/solvers/gp_solver.cpp
    ${CMAKE_SOURCE_DIR}/src/solvers/bogolyubov_solver.cpp
    ${CMAKE_SOURCE_DIR}/src/solvers/dipolar_bogolyubov_solver.cpp
    ${CMAKE_SOURCE_DIR}/src/solvers/dipolar_gp_solver.cpp
)

#-------------------------------
# Set the include directories
#-------------------------------

target_include_directories(solvers
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

#----------------
# Link libraries
#----------------

target_link_libraries(solvers mkl_wrappers)

#------------------------------------
# Link additional internal libraries
#------------------------------------

target_link_libraries(solvers utilities)

#---------------------
# Link to arpack-ng
#---------------------

target_include_directories(solvers PUBLIC ${ARPACK_DIR}/include/arpack)
target_link_directories(solvers PUBLIC ${ARPACK_DIR}/lib)
target_link_libraries(solvers ARPACK::ARPACK)

#-------------------------------------------------------
# Install library and headers to their final location
#-------------------------------------------------------

install(TARGETS solvers
        DESTINATION lib
        EXPORT UltraColdTargets
)

install(FILES gp_solvers.hpp DESTINATION include)
install(FILES bogolyubov_solvers.hpp DESTINATION include)
