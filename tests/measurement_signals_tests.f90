! measurement_signals_tests.f90

program main
    use signals_tests
    implicit none

    ! Local Variables
    logical :: local, overall

    ! Initialization
    overall = .true.

    ! Introduction
    print '(A)', ""
    print '(A)', "MEASUREMENT SIGNALS TESTS UNDERWAY..."

    ! Tests
    local = fourier_transform_test()
    if (.not.local) overall = .false.

    local = periodogram_test()
    if (.not.local) overall = .false.

    local = fourier_frequency_test()
    if (.not.local) overall = .false.

    local = low_pass_filter_test()
    if (.not.local) overall = .false.

    ! --------------------------------------------------------------------------
    ! End
    if (overall) then
        print '(A)', "MEASUREMENT SIGNALS TESTS COMPLETED SUCCESSFULLY."
    else
        print '(A)', "MEASUREMENT SIGNALS TESTS FAILED."
    end if
    print '(A)', ""
end program
