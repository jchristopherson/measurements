# Peak Detection Example
add_executable(peak_detect_example peak_detect_example.f90)
target_link_libraries(peak_detect_example measurements)

# Interpolation Example
add_executable(interpolation_example interpolation_example.f90)
target_link_libraries(interpolation_example measurements)

# C API Interpolation Example
include_directories(${PROJECT_SOURCE_DIR}/include)
add_executable(c_interpolation_example c_interpolation_example.c)
target_link_libraries(c_interpolation_example measurements)

# ANOVA Example
add_executable(anova_example anova_example.f90)
target_link_libraries(anova_example measurements)
