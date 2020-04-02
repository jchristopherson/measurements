# Define the MEASUREMENTS source files
set(measurements_sources
    measurements_core.f90
    measurement_stats.f90
    measurements_c_api.f90
)

# Build the library
add_library(measurements ${measurements_sources})
target_link_libraries(measurements 
    ${linalg_LIBRARIES} 
    ${ferror_LIBRARIES}
    ${nonlin_LIBRARIES}
)

# ------------------------------------------------------------------------------
# INSTALLATION INSTRUCTIONS
# ------------------------------------------------------------------------------
# Define target information
set_property(TARGET measurements PROPERTY VERSION ${MEASUREMENTS_VERSION})
set_property(TARGET measurements PROPERTY SOVERSION ${MEASUREMENTS_MAJOR_VERSION})
set_property(TARGET measurements PROPERTY INTERFACE_measurements_MAJOR_VERSION ${MEASUREMENTS_MAJOR_VERSION})
set_property(TARGET measurements APPEND PROPERTY COMPATIBLE_INTERFACE_STRING measurements_MAJOR_VERSION)

# Locate the "include" directories
set(measurements_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

# Define the installation instructions
install(TARGETS measurements EXPORT measurementsTargets
    RUNTIME DESTINATION measurements/bin
    LIBRARY DESTINATION measurements/lib
    ARCHIVE DESTINATION measurements/lib
    CONFIGURATIONS Release
    INCLUDES DESTINATION measurements/include
)
install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include
    DESTINATION ${CMAKE_INSTALL_PREFIX}/measurements
    COMPONENT Devel
)

# Include the documentation
install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/html DESTINATION measurements/doc)
install(FILES ${PROJECT_SOURCE_DIR}/doc/refman.pdf DESTINATION measurements/doc)

# Define the version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/measurementsConfigVersion.cmake"
    VERSION ${MEASUREMENTS_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(EXPORT measurementsTargets
    FILE "${CMAKE_BINARY_DIR}/measurementsTargets.cmake"
)

# Define the configuration file
configure_file(
    "${PROJECT_SOURCE_DIR}/measurementsConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/measurementsConfig.cmake"
    COPYONLY
)

set(ConfigPackageLocation measurements/lib/cmake/measurements)
install(
    EXPORT measurementsTargets
    FILE measurementsTargets.cmake
    DESTINATION ${ConfigPackageLocation}
)
install(
    FILES
        "${CMAKE_BINARY_DIR}/measurementsConfig.cmake"
        "${CMAKE_BINARY_DIR}/measurementsConfigVersion.cmake"
    DESTINATION
        ${ConfigPackageLocation}
    COMPONENT
        Devel
)