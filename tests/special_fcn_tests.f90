! special_fcn_tests.f90

module special_fcn_tests
    use iso_fortran_env
    use measurements_core
    implicit none

contains
! ------------------------------------------------------------------------------
! beta test
function beta_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8
    integer(int32), parameter :: npts = 10

    ! Local Variables
    integer(int32) :: i
    real(real64) :: a(npts), b(npts), x(npts), ans(npts), delta

    ! Initialization
    rst = .true.
    a = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0, 10.0d0]
    b = [10.0d0, 9.0d0, 8.0d0, 7.0d0, 6.0d0, 5.0d0, 4.0d0, 3.0d0, 2.0d0, 1.0d0]
    ans = [0.1d0, 0.0111111111111111000d0, 0.0027777777777777800d0, &
        0.0011904761904761900d0, 0.0007936507936507940d0, &
        0.0007936507936507940d0, 0.0011904761904761900d0, &
        0.0027777777777777800d0, 0.0111111111111111000d0, 0.1d0]

    ! Process
    x = beta(a, b)

    ! Tests
    do i = 1, npts
        delta = ans(i) - x(i)
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "BETA_TEST FAILED."
            print *, "Expected: ", ans(i)
            print *, "Computed: ", x(i)
            print *, "Difference (Expected - Computed): ", delta
        end if
    end do
end function

! ------------------------------------------------------------------------------
! regularized beta test
function reg_beta_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8
    integer(int32), parameter :: npts = 10
    real(real64), parameter :: a = 5.0d0
    real(real64), parameter :: b = 3.0d0

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(npts), ans(npts), f(npts), delta

    ! Initialization
    x = [0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, 0.9d0, 1.0d0]
    ans = [1.765d-4, 4.672d-3, 2.87955d-2, 9.6256d-2, 2.265625d-1, 4.19904d-1, &
        6.470694999999997d-1, 8.519679999999999d-1, 9.743085d-01, 1.0d0]

    ! Process
    f = regularized_beta(x, a, b)

    ! Tests
    do i = 1, npts
        delta = ans(i) - f(i)
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "REG_BETA_TEST FAILED."
            print *, "Expected: ", ans(i)
            print *, "Computed: ", f(i)
            print *, "Difference (Expected - Computed): ", delta
            print *, "X: ", x(i)
        end if
    end do
end function

! No need to test the incomplete beta function as the regularized beta function
! utilizes the incomplete beta function.

! ------------------------------------------------------------------------------
end module
