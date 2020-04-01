# MEASUREMENT_STATS TESTS
set(measurement_stats_files
    measurement_stats_tests.f90
    stats_tests.f90
)

# Build the executable
add_executable(MeasurementStatsTests ${measurement_stats_files})
target_link_libraries(MeasurementStatsTests
    measurements
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
