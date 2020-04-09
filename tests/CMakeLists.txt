# MEASUREMENT_STATS TESTS
set(measurement_stats_files
    measurement_stats_tests.f90
    stats_tests.f90
    special_fcn_tests.f90
)

# Build the executable
add_executable(MeasurementStatsTests ${measurement_stats_files})
target_link_libraries(MeasurementStatsTests
    measurements
    ${linalg_LIBRARIES}
)

# Run the tests
add_custom_command(
    OUTPUT measurement_stats_tests_output
    DEPENDS ${measurement_stats_files}
    COMMAND MeasurementStatsTests
)
add_custom_target(
    run_measurement_stats_tests ALL 
    DEPENDS measurement_stats_tests_output
)

# ------------------------------------------------------------------------------
# MEASUREMENTS_INTERP TESTS
set(measurement_interp_files
    measurement_interp_tests.f90
    interp_tests.f90
)

# Build the executable
add_executable(MeasurementInterpTests ${measurement_interp_files})
target_link_libraries(MeasurementInterpTests
    measurements
)

# Run the tests
add_custom_command(
    OUTPUT measurement_interp_tests_output
    DEPENDS ${measurement_interp_files}
    COMMAND MeasurementInterpTests
)
add_custom_target(
    run_measurement_interp_tests ALL 
    DEPENDS measurement_interp_tests_output
)

# ------------------------------------------------------------------------------
# C API MEASUREMENT STATS TESTS
# ------------------------------------------------------------------------------
set(c_measurement_stats_files
    c_stats_tests.c
    c_array.c
)

# Locate headers
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/tests)

# Build the executable
add_executable(C_MeasurementStatsTests ${c_measurement_stats_files})
target_link_libraries(C_MeasurementStatsTests
    measurements
    ${linalg_LIBRARIES}
)

# Run the tests
add_custom_command(
    OUTPUT c_measurement_stats_tests_output
    DEPENDS ${c_measurement_stats_files}
    COMMAND C_MeasurementStatsTests
)
add_custom_target(
    run_c_measurement_stats_tests ALL 
    DEPENDS c_measurement_stats_tests_output
)
