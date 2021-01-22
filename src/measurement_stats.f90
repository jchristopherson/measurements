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
    ! Compute the solution to: C = erf(z / sqrt(2)) for z.
    function zfun(x) result(f)
        real(real64), intent(in) :: x
        real(real64) :: f
        f = c - erf(x / sqrt(2.0d0))
    end function
end function

! ------------------------------------------------------------------------------
module function t_score(c, n, err) result(t)
    ! Arguments
    real(real64), intent(in) :: c
    integer(int32), intent(in) :: n
    class(errors), intent(inout), optional, target :: err
    real(real64) :: t

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
        t = 0.0d0
        call errmgr%report_error("t_score", "The confidence level must " // &
            "be a value between 0 and 1 (nonlinclusive).", &
            M_INVALID_INPUT_ERROR)
        return
    end if

    ! Set solver tolerances
    call solver%set_fcn_tolerance(1.0d-12)

    ! Compute the solution
    fcn => tfun
    call obj%set_fcn(fcn)
    t = 1.0d0
    lim%x1 = 0.0d0
    lim%x2 = 1.0d1
    call solver%solve(obj, t, lim, err = errmgr)
    
contains
    ! Compute the solution to: C = 1 - 1/2 Ix(v/2, 1/2) for t, where
    ! x = v / (v + t**2).
    function tfun(x) result(f)
        ! Arguments
        real(real64), intent(in) :: x
        real(real64) :: f
        
        ! Local Variables
        real(real64) :: v, xx
        
        ! Process
        v = n - 1.0d0
        xx = v / (v + x**2)
        f = c - (1.0d0 - 0.5d0 * regularized_beta(xx, 0.5d0 * v, 0.5d0))
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
pure elemental module function normal_distribution(mu, sigma, x) result(f)
    ! Arguments
    real(real64), intent(in) :: mu, sigma, x
    real(real64) :: f

    ! Constants
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Process
    f = (1.0d0 / (sigma * sqrt(2.0d0 * pi))) * &
        exp(-0.5 * ((x - mu) / sigma)**2)
end function

! ------------------------------------------------------------------------------
pure elemental module function t_distribution(dof, t) result(f)
    ! Arguments
    real(real64), intent(in) :: dof
    real(real64), intent(in) :: t
    real(real64) :: f

    ! Constants
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    real(real64) :: arg

    ! Process
    arg = 0.5d0 * (dof + 1)
    f = (gamma(arg) / (gamma(0.5d0 * dof) * sqrt(dof * pi))) * &
        (1.0d0 + t**2 / dof)**(-arg)
end function

! ------------------------------------------------------------------------------
pure elemental module function beta_distribution(a, b, x) result(z)
    ! Arguments
    real(real64), intent(in) :: a, b, x
    real(real64) :: z

    ! Process
    z = x**(a - 1.0d0) * (1.0d0 - x)**(b - 1.0d0) / beta(a, b)
end function

! ------------------------------------------------------------------------------
pure elemental module function f_distribution(d1, d2, x) result(z)
    ! Arguments
    real(real64), intent(in) :: d1, d2, x
    real(real64) :: z

    ! Local Variables
    real(real64) :: arg

    ! Process
    arg = ((d1 * x)**d1) * (d2**d2) / ((d1 * x + d2)**(d1 + d2))
    z = sqrt(arg) / (x * beta(0.5d0 * d1, 0.5d0 * d2))
end function

! ------------------------------------------------------------------------------
pure module function t_test(x1, x2, method) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x1, x2
    integer(int32), intent(in), optional :: method
    type(statistic) :: rst

    ! Local Variables
    integer(int32) :: flag

    ! Initialization
    flag = EQUAL_VARIANCE_ASSUMPTION
    if (present(method)) flag = method

    ! In the event that flag = PAIRED_DATA_SET_ASSUMPTION, and 
    ! size(x1) /= size(x2), default to EQUAL_VARIANCE_ASSUMPTION.
    if (flag == PAIRED_DATA_SET_ASSUMPTION .and. size(x1) /= size(x2)) then
        flag = EQUAL_VARIANCE_ASSUMPTION
    end if

    ! Process
    select case (flag)
        case (UNEQUAL_VARIANCE_ASSUMPTION)
            rst = t_test_no_same_var(x1, x2)
        case (PAIRED_DATA_SET_ASSUMPTION)
            rst = t_test_paired(x1, x2)
        case default
            rst = t_test_same_var(x1, x2)
    end select
end function

! ------------------------------------------------------------------------------
pure module function f_test(x1, x2) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x1, x2
    type(statistic) :: rst

    ! Local Variables
    integer(int32) :: d1, d2
    real(real64) :: var1, var2

    ! Process
    d1 = size(x1) - 1
    d2 = size(x2) - 1
    var1 = variance(x1)
    var2 = variance(x2)
    rst%value = var1 / var2
    rst%probability = ftest_probability(rst%value, d1, d2)
end function

! ------------------------------------------------------------------------------
pure module function remove_nans(x) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64), allocatable, dimension(:) :: rst

    ! Local Variables
    integer(int32) :: i, j, n
    real(real64), allocatable, dimension(:) :: buffer

    ! Process
    n = size(x)
    allocate(buffer(n))
    j = 0
    do i = 1, n
        if (.not.isnan(x(i))) then
            j = j + 1
            buffer(j) = x(i)
        end if
    end do
    rst = buffer(1:j)
end function

! ------------------------------------------------------------------------------
pure module function remove_zeros(x, tol) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in), optional :: tol
    real(real64), allocatable, dimension(:) :: rst

    ! Local Variables
    integer(int32) :: i, j, n
    real(real64), allocatable, dimension(:) :: buffer
    real(real64) :: t

    ! Process
    n = size(x)
    allocate(buffer(n))
    t = 2.0d0 * epsilon(t)
    if (present(tol)) t = tol
    j = 0
    do i = 1, n
        if (abs(x(i)) > t) then
            j = j + 1
            buffer(j) = x(i)
        end if
    end do
    rst = buffer(1:j)
end function

! ------------------------------------------------------------------------------
module function r_squared(y, ym, err) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: y, ym
    class(errors), intent(inout), optional, target :: err
    real(real64) :: rst

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: esum, vt
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Setting up error handling
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Ensure the input arrays are properly sized
    n = size(y)
    if (size(ym) /= n) then
        call errmgr%report_error("r_squared", "The input arrays must " // &
            "be the same size.", M_ARRAY_SIZE_ERROR)
        return
    end if

    ! Compute the sum of the errors squared
    esum = 0.0d0
    do i = 1, n
        esum = esum + (y(i) - ym(i))**2
    end do

    ! Compute the total variance
    vt = variance(y) * (n - 1.0d0)

    ! Compute the r-squared as 1 - esum / vt
    rst = 1.0d0 - esum / vt
end function

! ------------------------------------------------------------------------------
pure module function ftest_probability(f, dof1, dof2) result(rst)
    ! Arguments
    real(real64), intent(in) :: f
    integer(int32), intent(in) :: dof1, dof2

    ! Local Variables
    real(real64) :: arg, d1, d2

    ! Process
    d1 = real(min(dof1, dof2), real64)
    d2 = real(max(dof1, dof2), real64)
    arg = d1 * f / (d1 * f + d2)
    rst = regularized_beta(arg, 0.5d0 * d1, 0.5d0 * d2)
end function

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
pure function t_test_same_var(x1, x2) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x1, x2
    type(statistic) :: rst

    ! Local Variables
    integer(int32) :: n1, n2
    real(real64) :: var1, var2, avg1, avg2, svar, df

    ! Process
    n1 = size(x1)
    n2 = size(x2)
    df = n1 + n2 - 2.0d0
    avg1 = mean(x1)
    avg2 = mean(x2)
    var1 = variance(x1)
    var2 = variance(x2)
    svar = ((n1 - 1.0d0) * var1 + (n2 - 1.0d0) * var2) / df
    rst%value = (avg1 - avg2) / sqrt(svar * (1.0d0 / n1 + 1.0d0 / n2))
    rst%probability = regularized_beta(&
        df / (df + rst%value**2), &
        0.5d0 * df, &
        0.5d0)
end function

! ------------------------------------------------------------------------------
pure function t_test_no_same_var(x1, x2) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x1, x2
    type(statistic) :: rst

    ! Local Variables
    integer(int32) :: n1, n2
    real(real64) :: var1, var2, avg1, avg2, df

    ! Process
    n1 = size(x1)
    n2 = size(x2)
    avg1 = mean(x1)
    avg2 = mean(x2)
    var1 = variance(x1)
    var2 = variance(x2)
    rst%value = (avg1 - avg2) / sqrt(var1 / n1 + var2 / n2)
    df = (var1 / n1 + var2 / n2)**2 / &
        ((var1 / n1)**2 / (n1 - 1.0d0) + (var2 / n2)**2 / (n2 - 1.0d0))
    rst%probability = regularized_beta(&
        df / (df + rst%value**2), &
        0.5d0 * df, &
        0.5d0)
end function

! ------------------------------------------------------------------------------
pure function t_test_paired(x1, x2) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x1, x2
    type(statistic) :: rst

    ! Local Variables
    integer(int32) :: j, n
    real(real64) :: var1, var2, avg1, avg2, df, sd, cov

    ! Process
    n = size(x1)
    avg1 = mean(x1)
    avg2 = mean(x2)
    var1 = variance(x1)
    var2 = variance(x2)
    df = n - 1.0d0
    cov = 0.0d0
    do j = 1, n
        cov = cov + (x1(j) - avg1) * (x2(j) - avg2)
    end do
    cov = cov / df
    sd = sqrt((var1 + var2 - 2.0d0 * cov) / n)
    rst%value = (avg1 - avg2) / sd
    rst%probability = regularized_beta(&
        df / (df + rst%value**2), &
        0.5d0 * df, &
        0.5d0)
end function

! ------------------------------------------------------------------------------
end submodule
