! measurement_stats.f90

submodule (measurements_core) measurement_stats
    use linalg_core
    use nonlin_core
    use nonlin_solve
contains
! ------------------------------------------------------------------------------
pure module function is_monotonic(x) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    logical :: rst

    ! Local Variables
    integer(int32) :: i
    logical :: ascend

    ! Quick Return
    rst = .true.
    if (size(x) <= 1) return
    if (x(2) == x(1)) then
        rst = .false.
        return
    end if

    ! Process
    ascend = x(2) > x(1)
    if (ascend) then
        do i = 2, size(x)
            if (x(i) <= x(i-1)) then
                rst = .false.
                exit
            end if
        end do
    else
        do i = 2, size(x)
            if (x(i) >= x(i-1)) then
                rst = .false.
                exit
            end if
        end do
    end if
end function

! ------------------------------------------------------------------------------
pure module function mean(x) result(z)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64) :: z

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(x)
    if (n == 0) then
        z = zero
    else
        z = x(1)
        do i = 2, n
            z = z + (x(i) - z) / i
        end do
    end if
end function

! ------------------------------------------------------------------------------
module function median(x) result(z)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64) :: z

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: half = 0.5d0

    ! Local Variables
    integer(int32) :: n, nmid, nmidp1, iflag
    real(real64), allocatable, dimension(:) :: y

    ! Initialization
    n = size(x)
    nmid = n / 2
    nmidp1 = nmid + 1
    iflag = n - 2 * nmid

    ! Sort the data into ascending order
    y = x
    call sort(y, .true.)

    ! Compute the median
    if (iflag == 0) then
        z = half * (y(nmid) + y(nmidp1))
    else
        z = y(nmidp1)
    end if
end function

! ------------------------------------------------------------------------------
pure module function variance(x) result(v)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64) :: v

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: oldMean, newMean

    ! Process
    n = size(x)
    if (n <= 1) then
        v = zero
    else
        oldMean = x(1)
        v = zero
        do i = 2, n
            newMean = oldMean + (x(i) - oldMean) / i
            v = v + (x(i) - oldMean) * (x(i) - newMean)
            oldMean = newMean
        end do
        v = v / (n - one)
    end if
end function

! ------------------------------------------------------------------------------
pure module function standard_deviation(x) result(s)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64) :: s

    ! Process
    s = sqrt(variance(x))
end function

! ------------------------------------------------------------------------------
pure module function data_range(x) result(r)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64) :: r

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: maxx, minx

    ! Initialization
    r = 0.0d0
    n = size(x)

    ! Quick Return
    if (n <= 1) return

    ! Process
    maxx = x(1)
    minx = x(1)
    do i = 2, n
        if (x(i) > maxx) maxx = x(i)
        if (x(i) < minx) minx = x(i)
    end do
    r = maxx - minx
end function

! ------------------------------------------------------------------------------
module function z_score(c, err) result(z)
    ! Arguments
    real(real64), intent(in) :: c
    class(errors), intent(inout), optional, target :: err
    real(real64) :: z

    ! Local Variables
    type(fcn1var_helper) :: obj
    procedure(fcn1var), pointer :: fcn
    type(brent_solver) :: solver
    type(value_pair) :: lim
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (c <= 0.0d0 .or. c >= 1.0d0) then
        z = 0.0d0
        call errmgr%report_error("z_score", "The confidence level must " // &
            "be a value between 0 and 1 (nonlinclusive).", &
            M_INVALID_INPUT_ERROR)
        return
    end if

    ! Set solver tolerances
    call solver%set_fcn_tolerance(1.0d-12)

    ! Compute the solution
    fcn => zfun
    call obj%set_fcn(fcn)
    z = 1.0d0
    lim%x1 = 0.0d0
    lim%x2 = 1.0d1
    call solver%solve(obj, z, lim, err = errmgr)
    
contains
    ! Compute the solution to: alpha = erf(z / sqrt(2)) for z.
    function zfun(x) result(f)
        real(real64), intent(in) :: x
        real(real64) :: f
        f = c - erf(x / sqrt(2.0d0))
    end function
end function

! ------------------------------------------------------------------------------
pure module function confidence_interval(x, zval) result(ci)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in) :: zval
    real(real64) :: ci

    ! Local Variables
    real(real64) :: stdev
    integer(int32) :: n

    ! Compute the standard deviation
    stdev = standard_deviation(x)

    ! Compute the interval
    n = size(x)
    ci = zval * stdev / sqrt(real(n))
end function

! ------------------------------------------------------------------------------
pure elemental module function normal_distribution(mu, sigma, x, comp) result(f)
    ! Arguments
    real(real64), intent(in) :: mu, sigma, x
    logical, intent(in), optional :: comp
    real(real64) :: f

    ! Constants
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    logical :: check

    ! Initialization
    if (present(comp)) then
        check = comp
    else
        check = .false.
    end if

    ! Process
    if (check) then
        f = 0.5d0 * (1.0d0 + erf((x - mu) / (sigma * sqrt(2.0d0))))
    else
        f = (1.0d0 / (sigma * sqrt(2.0d0 * pi))) * &
            exp(-0.5 * ((x - mu) / sigma)**2)
    end if
end function

! ------------------------------------------------------------------------------
pure elemental module function t_distribution(dof, t, comp) result(f)
    ! Arguments
    integer(int32), intent(in) :: dof
    real(real64), intent(in) :: t
    logical, intent(in), optional :: comp
    real(real64) :: f

    ! Constants
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    real(real64) :: arg1
    logical :: check

    ! Initialization
    if (present(comp)) then
        check = comp
    else
        check = .false.
    end if

    ! Process
    if (check) then
        ! TO DO: Need hypergeometric function first
    else
        arg1 = 0.5d0 * (dof + 1)
        f = (gamma(arg1) / (gamma(0.5d0 * dof) * sqrt(dof * pi))) * &
            (1.0d0 + t**2 / dof)**(-arg1)
    end if
end function

! ------------------------------------------------------------------------------
pure elemental module function beta(a, b) result(z)
    ! Arguments
    real(real64), intent(in) :: a, b
    real(real64) :: z

    ! Process
    ! REF: https://en.wikipedia.org/wiki/Beta_function
    z = gamma(a) * gamma(b) / gamma(a + b)
end function

! ------------------------------------------------------------------------------
pure elemental module function beta_distribution(a, b, x, comp) result(z)
    ! Arguments
    real(real64), intent(in) :: a, b, x
    logical, intent(in), optional :: comp
    real(real64) :: z

    ! Local Variables
    logical :: check

    ! Initialization
    if (present(comp)) then
        check = comp
    else
        check = .false.
    end if

    ! Process
    if (check) then
        z = beta_distribution(a, b, x) / beta(a, b)
    else
        z = x**(a - 1.0d0) * (1.0d0 - x)**(b - 1.0d0) / beta(a, b)
    end if
end function

! ------------------------------------------------------------------------------
pure elemental module function incomplete_beta(x, a, b) result(z)
    ! Arguments
    real(real64), intent(in) :: x, a, b
    real(real64) :: z

    ! Process
    ! REF: https://en.wikipedia.org/wiki/Beta_function
    z = beta_distribution(a, b, x) * beta(a, b)
end function

! ------------------------------------------------------------------------------
! incomplete gamma function

! ------------------------------------------------------------------------------
! hypergeometric function
! Need PSI before HYGFX
! https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
! https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
pure elemental module function f_distribution(d1, d2, x, comp) result(z)
    ! Arguments
    real(real64), intent(in) :: d1, d2, x
    logical, intent(in), optional :: comp
    real(real64) :: z

    ! Local Variables
    real(real64) :: arg, betaReg
    logical :: check

    ! Initialization
    if (present(comp)) then
        check = comp
    else
        check = .false.
    end if

    ! Process
    if (check) then
        arg = d1 * x / (d1 * x + d2)
        betaReg = beta_distribution(arg, 0.5d0 * d1, 0.5d0 * d2) / &
            beta(0.5d0 * d1, 0.5d0 * d2)
    else
        arg = ((d1 * x)**d1) * (d2**d2) / ((d1 * x + d2)**(d1 + d2))
        z = sqrt(arg) / (x * beta(0.5d0 * d1, 0.5d0 * d2))
    end if
end function

! ------------------------------------------------------------------------------
! f test

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule
