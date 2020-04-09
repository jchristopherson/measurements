! interp_tests.f90

module interp_tests
    use iso_fortran_env
    use measurements_core
    implicit none
contains
! ------------------------------------------------------------------------------
function linear_interp_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8
    integer(int32), parameter :: npts = 21
    integer(int32), parameter :: ni = 15

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(npts), y(npts), xi(ni), yi(ni), ans(ni), delta
    type(linear_interp) :: interp

    ! Initialization
    rst = .true.
    x = [0.0d0, 0.05d0, 0.1d0, 0.15d0, 0.2d0, 0.25d0, 0.3d0, 0.35d0, 0.4d0, &
        0.45d0, 0.5d0, 0.55d0, 0.6d0, 0.65d0, 0.7d0, 0.75d0, 0.8d0, 0.85d0, &
        0.9d0, 0.95d0, 1.0d0]
    y = [0.0d0, 0.30901699437494700d0, 0.58778525229247300d0, &
        0.80901699437494700d0, 0.95105651629515400d0, 1.0d0, &
        0.95105651629515400d0, 0.80901699437494700d0, 0.58778525229247300d0, &
        0.30901699437494800d0, 0.0d0, -0.30901699437494700d0, &
        -0.58778525229247300d0, -0.80901699437494700d0, &
        -0.95105651629515400d0, -1.0d0, -0.95105651629515300d0, &
        -0.80901699437494700d0, -0.58778525229247200d0, &
        -0.30901699437494600d0, 0.0d0]
    xi = [0.0d0, 0.075d0, 0.15d0, 0.225d0, 0.3d0, 0.375d0, 0.45d0, 0.525d0, &
        0.6d0, 0.675d0, 0.75d0, 0.825d0, 0.9d0, 0.975d0, 1.0d0]
    ans = [0.0d0, 0.4484011230d0, 0.8090169940d0, 0.9755282580d0, &
        0.9510565160d0, 0.6984011230d0, 0.3090169940d0, -0.1545084970d0, &
        -0.5877852520d0, -0.8800367550d0, -1.0d0, -0.88003675500d0, &
        -0.58778525200d0, -0.154508496999999d0, 0.0d0]

    ! Process
    call interp%initialize(x, y)
    yi = interp%interpolate(xi)

    ! Test
    do i = 1, ni
        delta = ans(i) - yi(i)
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "LINEAR_INTERP_TEST FAILED."
            print *, "Expected: ", ans(i)
            print *, "Computed: ", yi(i)
            print *, "Difference (Expected - Computed): ", delta
            print *, "Index: ", i
        end if
    end do
end function

! ------------------------------------------------------------------------------
function polynomial_interp_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 5.0d-3
    integer(int32), parameter :: npts = 21
    integer(int32), parameter :: ni = 15

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(npts), y(npts), xi(ni), yi(ni), ans(ni), delta
    type(polynomial_interp) :: interp

    ! Initialization
    rst = .true.
    x = [0.0d0, 0.05d0, 0.1d0, 0.15d0, 0.2d0, 0.25d0, 0.3d0, 0.35d0, 0.4d0, &
        0.45d0, 0.5d0, 0.55d0, 0.6d0, 0.65d0, 0.7d0, 0.75d0, 0.8d0, 0.85d0, &
        0.9d0, 0.95d0, 1.0d0]
    y = [0.0d0, 0.30901699437494700d0, 0.58778525229247300d0, &
        0.80901699437494700d0, 0.95105651629515400d0, 1.0d0, &
        0.95105651629515400d0, 0.80901699437494700d0, 0.58778525229247300d0, &
        0.30901699437494800d0, 0.0d0, -0.30901699437494700d0, &
        -0.58778525229247300d0, -0.80901699437494700d0, &
        -0.95105651629515400d0, -1.0d0, -0.95105651629515300d0, &
        -0.80901699437494700d0, -0.58778525229247200d0, &
        -0.30901699437494600d0, 0.0d0]
    xi = [0.0d0, 0.075d0, 0.15d0, 0.225d0, 0.3d0, 0.375d0, 0.45d0, 0.525d0, &
        0.6d0, 0.675d0, 0.75d0, 0.825d0, 0.9d0, 0.975d0, 1.0d0]
    ans = [0.0d0, 0.45420421588342d0, 0.80901699400000d0, 0.98462842510122d0, &
        0.95105651600000d0, 0.70761184418907d0, 0.30901699400000d0, &
        -0.15649633470076d0, -0.58778525200000d0, -0.89256206037553d0, &
        -1.0d0, -0.892562060375530d0, -0.587785252000000d0, &
        -0.158386880700761d0, 0.0d0]

    ! Process
    call interp%initialize(x, y, 3)
    yi = interp%interpolate(xi)

    ! Test
    do i = 1, ni
        delta = ans(i) - yi(i)
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "POLYNOMIAL_INTERP_TEST FAILED."
            print *, "Expected: ", ans(i)
            print *, "Computed: ", yi(i)
            print *, "Difference (Expected - Computed): ", delta
            print *, "Index: ", i
        end if
    end do
end function

! ------------------------------------------------------------------------------
function spline_interp_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 2.0d-3
    integer(int32), parameter :: npts = 21
    integer(int32), parameter :: ni = 15

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(npts), y(npts), xi(ni), yi(ni), ans(ni), delta
    type(spline_interp) :: interp

    ! Initialization
    rst = .true.
    x = [0.0d0, 0.05d0, 0.1d0, 0.15d0, 0.2d0, 0.25d0, 0.3d0, 0.35d0, 0.4d0, &
        0.45d0, 0.5d0, 0.55d0, 0.6d0, 0.65d0, 0.7d0, 0.75d0, 0.8d0, 0.85d0, &
        0.9d0, 0.95d0, 1.0d0]
    y = [0.0d0, 0.30901699437494700d0, 0.58778525229247300d0, &
        0.80901699437494700d0, 0.95105651629515400d0, 1.0d0, &
        0.95105651629515400d0, 0.80901699437494700d0, 0.58778525229247300d0, &
        0.30901699437494800d0, 0.0d0, -0.30901699437494700d0, &
        -0.58778525229247300d0, -0.80901699437494700d0, &
        -0.95105651629515400d0, -1.0d0, -0.95105651629515300d0, &
        -0.80901699437494700d0, -0.58778525229247200d0, &
        -0.30901699437494600d0, 0.0d0]
    xi = [0.0d0, 0.075d0, 0.15d0, 0.225d0, 0.3d0, 0.375d0, 0.45d0, 0.525d0, &
        0.6d0, 0.675d0, 0.75d0, 0.825d0, 0.9d0, 0.975d0, 1.0d0]
    ans = [0.0d0, 0.45395574322220d0, 0.80901699400000d0, 0.98766310413990d0, &
        0.95105651600000d0, 0.70708838900300d0, 0.30901699400000d0, &
        -0.15643039704573d0, -0.58778525200000d0, -0.89098339112519d0, -1.0d0, &
        -0.890981711444173d0, -0.587785252000000d0, -0.156516060777801d0, 0.0d0]

    ! Process
    call interp%initialize(x, y)
    yi = interp%interpolate(xi)

    ! Test
    do i = 1, ni
        delta = ans(i) - yi(i)
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "SPLINE_INTERP_TEST FAILED."
            print *, "Expected: ", ans(i)
            print *, "Computed: ", yi(i)
            print *, "Difference (Expected - Computed): ", delta
            print *, "Index: ", i
        end if
    end do
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module