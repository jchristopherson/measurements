! measurement_stats_tests.f90

program main
    use iso_fortran_env
    use stats_tests
    implicit none

    ! Local Variables
    logical :: local, overall

    ! Initialization
    overall = .true.

    ! Introduction
    print '(A)', ""
    print '(A)', "MEASUREMENT STATS TESTS UNDERWAY..."

    ! Tests
    local = mean_test()
    if (.not.local) overall = .false.

    local = median_test_even()
    if (.not.local) overall = .false.

    local = median_test_odd()
    if (.not.local) overall = .false.

    local = variance_test()
    if (.not.local) overall = .false.

    local = std_dev_test()
    if (.not.local) overall = .false.

    local = range_test()
    if (.not.local) overall = .false.

    local = z_score_test()
    if (.not.local) overall = .false.

    local = confidence_interval_test()
    if (.not.local) overall = .false.

    local = normal_distribution_test()
    if (.not.local) overall = .false.

    local = t_distribution_test()
    if (.not.local) overall = .false.

    ! beta

    local = beta_distribution_test()
    if (.not.local) overall = .false.

    ! incomplete beta

    local = f_distribution_test()
    if (.not.local) overall = .false.

    ! End
    if (overall) then
        print '(A)', "MEASUREMENT STATS TESTS COMPLETED SUCCESSFULLY."
    else
        print '(A)', "MEASUREMENT STATS TESTS FAILED."
    end if
    print '(A)', ""
end program
