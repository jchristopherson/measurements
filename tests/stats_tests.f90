! stats_tests.f90

module stats_tests
    use iso_fortran_env
    use measurements_core
    use linalg_core
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
function median_test_even() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), x1, x2, ans, computed, delta

    ! Initialization
    rst = .true.

    ! Populate an array with random values, and then sort into ascending
    ! order
    call random_number(x)
    call sort(x, .true.)

    ! Compute the answer
    x1 = x(npts / 2)
    x2 = x(npts / 2 + 1)
    ans = 0.5d0 * (x1 + x2)

    ! Process
    computed = median(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "MEDIAN_TEST_EVEN FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function median_test_odd() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 10001
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), ans, computed, delta

    ! Initialization
    rst = .true.

    ! Populate an array with random values, and then sort into ascending
    ! order
    call random_number(x)
    call sort(x, .true.)

    ! Compute the answer
    ans = x(floor(npts / 2.0d0) + 1)

    ! Process
    computed = median(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "MEDIAN_TEST_ODD FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function variance_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), ans, computed, delta, avg

    ! Initialization
    rst = .true.

    ! Create a random data set
    call random_number(x)

    ! Compute the solution
    avg = mean(x)
    ans = sum((x - avg)**2) / (npts - 1.0d0)

    ! Process
    computed = variance(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "VARIANCE_TEST FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function std_dev_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), ans, computed, delta, avg

    ! Initialization
    rst = .true.

    ! Create a random data set
    call random_number(x)

    ! Compute the solution
    avg = mean(x)
    ans = sqrt(sum((x - avg)**2) / (npts - 1.0d0))

    ! Process
    computed = standard_deviation(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "STANDARD_DEVIATION_TEST FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function range_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), ans, computed, delta

    ! Initialization
    rst = .true.

    ! Create a random data set
    call random_number(x)

    ! Compute the solution
    ans = maxval(x) - minval(x)

    ! Process
    computed = data_range(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "RANGE_TEST FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function z_score_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameter
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: tol = 1.0d-8
    real(real64), parameter :: z80 = 1.281551565545d0
    real(real64), parameter :: z90 = 1.644853626951d0
    real(real64), parameter :: z95 = 1.959963984540d0
    real(real64), parameter :: z99 = 2.575829303549d0
    real(real64), parameter :: z999 = 3.290526731492d0

    ! Local Variables
    real(real64) :: ans, computed, delta

    ! Initialization
    rst = .true.

    ! 80%
    ans = z80
    computed = z_score(0.8d0)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "Z_SCORE_TEST 80% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! 90%
    ans = z90
    computed = z_score(0.9d0)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "Z_SCORE_TEST 90% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! 95%
    ans = z95
    computed = z_score(0.95d0)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "Z_SCORE_TEST 95% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! 99%
    ans = z99
    computed = z_score(0.99d0)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "Z_SCORE_TEST 99% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! 99.9%
    ans = z999
    computed = z_score(0.999d0)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "Z_SCORE_TEST 99.9% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
! confidence interval test

! ------------------------------------------------------------------------------
! normal distribution test

! ------------------------------------------------------------------------------
! t distribution test

! ------------------------------------------------------------------------------
end module
