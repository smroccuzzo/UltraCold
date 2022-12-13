###################################################################
# This files contains instructions for cmake to install UltraCold
# to its final location, specified in CMAKE_INSTALL_PREFIX
###################################################################

#---------------------------------------------------------
# Install the fundamental header to the final destination
#---------------------------------------------------------

install(FILES ${CMAKE_SOURCE_DIR}/source/UltraCold.hpp DESTINATION include)

#------------------------------------------------
# Install the config file to its final location
#------------------------------------------------

install(EXPORT UltraColdTargets
        FILE UltraColdTargets.cmake
        DESTINATION cmake
        )

include(CMakePackageConfigHelpers)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/UltraColdTargets.cmake"
        INSTALL_DESTINATION "cmake"
        NO_SET_AND_CHECK_MACRO
        NO_CHECK_REQUIRED_COMPONENTS_MACRO
        )
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/UltraColdTargets.cmake
        DESTINATION cmake
        )
install(FILES
        ${CMAKE_SOURCE_DIR}/cmake/UltraColdConfig.cmake
        DESTINATION cmake
        )
export(EXPORT UltraColdTargets
        FILE "${CMAKE_CURRENT_BINARY_DIR}/UltraColdTargets.cmake")