########################################################################
# Build a library containing functions and classes of general utility
########################################################################

add_library(utilities
        STATIC
        ${CMAKE_SOURCE_DIR}/source/base/utilities/HardwareInspector.cpp
        ${CMAKE_SOURCE_DIR}/source/base/utilities/InputParser.cpp
        ${CMAKE_SOURCE_DIR}/source/base/utilities/DataWriter.cpp
        ${CMAKE_SOURCE_DIR}/source/base/utilities/Vector.cpp
        ${CMAKE_SOURCE_DIR}/source/base/utilities/DFtCalculator.cpp
        ${CMAKE_SOURCE_DIR}/source/base/utilities/mesh_fourier_space.cpp
        ${CMAKE_SOURCE_DIR}/source/base/utilities/Colors.cpp
        )

target_include_directories(utilities
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
        )

#-------------------------------------------------------
# Install library and headers to their final location
#-------------------------------------------------------

install(TARGETS utilities
        DESTINATION lib
        EXPORT UltraColdTargets
        )

install(FILES HardwareInspector.hpp DESTINATION include)
install(FILES InputParser.hpp DESTINATION include)
install(FILES Vector.hpp DESTINATION include)
install(FILES DataWriter.hpp DESTINATION include)
install(FILES mesh_fourier_space.hpp DESTINATION include)
install(FILES DFtCalculator.hpp DESTINATION include)
install(FILES Colors.hpp DESTINATION include)
