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
function confidence_interval_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: alpha = 0.95d0
    integer(int32), parameter :: npts = 10
    real(real64), parameter :: ans = 0.196943590756784000d0
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: zval, delta, computed, x(npts)

    ! Initialization
    rst = .true.
    x = [0.1266904993777120d0, 0.3827711768054340d0, 0.1370953805850570d0, &
        0.5852213531153070d0, 0.2267533281658030d0, 0.0999861358308985d0, &
        0.5851003510284570d0, 0.8136628645855180d0, 0.7400357894369070d0, &
        0.9787774755208680d0]

    ! Process
    zval = z_score(alpha)
    computed = confidence_interval(x, zval)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "CONFIDENCE_INTERVAL_TEST FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
! normal distribution test
function normal_distribution_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: avg = 5.264049719320700d-1
    real(real64), parameter :: sigma = 2.110340893383650d-1
    integer(int32), parameter :: npts = 21
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: delta, x(npts), ans(npts), computed(npts)
    integer(int32) :: i

    ! Initialization
    rst = .true.
    x = [-5.0d0, -4.5d0, -4.0d0, -3.5d0, -3.0d0, -2.5d0, -2.0d0, -1.5d0, &
        -1.0d0, -0.5d0, 0.0d0, 0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, &
        3.5d0, 4.0d0, 4.5d0, 5.0d0]
    ans = [2.306272166314290d-149, 1.229690115816890d-123, &
        2.392009083889280d-100, 1.697508606279160d-79, 4.394840953499330d-61, &
        4.151034651360500d-45, 1.430380472227950d-31, 1.798161951816050d-20, &
        8.246849202359420d-12, 1.379841874725850d-5, 8.422724695133530d-02, &
        1.875676375640870d0, 1.523860560727160d-1, 4.516630549640910d-5, &
        4.883890565618130d-11, 1.926634346379860d-19, 2.772775386838060d-30, &
        1.455834812310320d-43, 2.788634062219750d-59, 1.948735843727060d-77, &
        4.968169870317840d-98]

    ! Process
    computed = normal_distribution(avg, sigma, x)

    ! Test
    do i = 1, npts
        delta = ans(i) - computed(i)
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "NORMAL_DISTRIBUTION_TEST FAILED."
            print *, "Expected: ", ans(i)
            print *, "Computed: ", computed(i)
            print *, "Difference (Expected - Computed): ", delta
            print *, "Index: ", i
        end if
    end do
end function

! ------------------------------------------------------------------------------
! t distribution test

! ------------------------------------------------------------------------------
! beta function test

! ------------------------------------------------------------------------------
! beta distribution test

! ------------------------------------------------------------------------------
! incomplete beta function test

! ------------------------------------------------------------------------------
! f distribution test

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
