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

# ANOVA Example 2
add_executable(anova_example_2 anova_example_2.f90)
target_link_libraries(anova_example_2 measurements)

# C API ANOVA Example
add_executable(c_anova_example c_anova_example.c)
target_link_libraries(c_anova_example measurements)

# Distribution Example
find_package(fplot)
find_package(fcore)
add_executable(distribution_example distribution_example.f90)
target_link_libraries(distribution_example measurements fplot fcore)

# Special Function Example
add_executable(specials_example specials_example.f90)
target_link_libraries(specials_example measurements fplot fcore)
