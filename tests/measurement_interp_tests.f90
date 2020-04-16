! measurement_interp_tests.f90

program main
    use interp_tests
    use regression_test
    implicit none

    ! Local Variables
    logical :: local, overall

    ! Initialization
    overall = .true.

    ! Introduction
    print '(A)', ""
    print '(A)', "MEASUREMENT INTERPOLATION TESTS UNDERWAY..."

    ! Tests
    local = linear_interp_test()
    if (.not.local) overall = .false.

    local = polynomial_interp_test()
    if (.not.local) overall = .false.

    local = spline_interp_test()
    if (.not.local) overall = .false.

    local = llsq_mimo_test()
    if (.not.local) overall = .false.

    local = llsq_miso_test()
    if (.not.local) overall = .false.

    local = peak_detect_test()
    if (.not.local) overall = .false.

    ! --------------------------------------------------------------------------
    ! End
    if (overall) then
        print '(A)', "MEASUREMENT INTERPOLATION TESTS COMPLETED SUCCESSFULLY."
    else
        print '(A)', "MEASUREMENT INTERPOLATION TESTS FAILED."
    end if
    print '(A)', ""
end program
