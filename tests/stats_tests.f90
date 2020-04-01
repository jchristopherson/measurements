! stats_tests.f90

module stats_tests
    use iso_fortran_env
    use measurements_core
    implicit none
contains
! ------------------------------------------------------------------------------
function mean_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 100000
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), ans, computed, delta

    ! Initialization
    rst = .true.
    call random_number(x)

    ! Compute the actual solution the old-school way
    ans = sum(x) / npts

    ! Utilize the actual method
    computed = mean(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "MEAN_TEST FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
