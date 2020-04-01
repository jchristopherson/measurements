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

    ! End
    if (overall) then
        print '(A)', "MEASUREMENT STATS TESTS COMPLETED SUCCESSFULLY."
    else
        print '(A)', "MEASUREMENT STATS TESTS FAILED."
    end if
    print '(A)', ""
end program
