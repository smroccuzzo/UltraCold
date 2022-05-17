########################################################################
# Build a library containing functions and classes of general utility
########################################################################

add_library(utilities
        STATIC
        ${CMAKE_SOURCE_DIR}/src/utilities/hardware_inspector.cpp
        ${CMAKE_SOURCE_DIR}/src/utilities/input_parser.cpp
        ${CMAKE_SOURCE_DIR}/src/utilities/data_out.cpp
        ${CMAKE_SOURCE_DIR}/src/utilities/vector.cpp
        ${CMAKE_SOURCE_DIR}/src/utilities/mesh_fourier_space.cpp
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

install(FILES hardware_inspector.hpp DESTINATION include)
install(FILES input_parser.hpp DESTINATION include)
install(FILES vector.hpp DESTINATION include)
install(FILES data_out.hpp DESTINATION include)
install(FILES mesh_fourier_space.hpp DESTINATION include)
