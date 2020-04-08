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
! ANOVA type GR&R
pure module function anova(x, alpha) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:,:,:) :: x
    real(real64), intent(in), optional :: alpha
    type(process_variance) :: rst

    ! Local Variables
    integer(int32) :: npart, nrepeats, nops, partDOF, opDOF, partOpDOF, &
        rptDOF, totalDOF
    real(real64) :: a, overallAvg, ssqPart, ssqOp, ssqWithin, ssqTotal, &
        ssqPartOp, msqPart, msqOp, msqPartOp, msqRpt, msqTotal, opFtest, &
        partOpProb

    ! Initialization
    npart = size(x, 1)
    nrepeats = size(x, 2)
    nops = size(x, 3)
    a = 5.0d-2
    if (present(alpha)) a = alpha

    ! Quick Return
    if (npart < 2 .or. nrepeats < 2 .or. nops < 1) then
        rst%measurement_variance = 0.0d0
        rst%part_variance = 0.0d0
        rst%total_variance = 0.0d0
        rst%equipment_variance = 0.0d0
        rst%operator_variance = 0.0d0
        rst%operator_by_part_variance = 0.0d0
        return
    end if

    ! Compute the overall mean
    overallAvg = mean(pack(x, .true.))

    ! Compute the SSQ terms
    ssqPart = ssq_part(x, overallAvg)
    ssqOp = ssq_operator(x, overallAvg)
    ssqWithin = ssq_repeat(x)
    ssqTotal = ssq_total(x, overallAvg)
    ssqPartOp = ssqTotal - ssqPart - ssqOp - ssqWithin

    ! Compute the DOF terms
    partDOF = npart - 1
    opDOF = nops - 1
    partOpDOF = partDOF * opDOF
    rptDOF = npart * nops * (nrepeats - 1)
    totalDOF = npart * nops * nrepeats - 1

    ! Compute the mean of the squared differences
    msqPart = ssqPart / partDOF
    msqOp = ssqOp / opDOF
    msqPartOp = ssqPartOp / partOpDOF
    msqTotal = ssqTotal / totalDOF
    msqRpt = ssqWithin / rptDOF

    ! Compute the F tests
    opFtest = msqPartOp / msqRpt

    ! Compute the probability terms
    partOpProb = 1.0d0 - f_distribution(real(partOpDOF, real64), &
        real(rptDOF, real64), opFtest)

    ! Compute the variance terms
    if (partOpProb > a) then
        rst%part_variance = max(0.0d0, (msqPart - msqPartOp) / &
            (nops * nrepeats))
        rst%operator_variance = max(0.0d0, (msqOp - msqPartOp) / &
            (npart * nrepeats))
    else
        rst%part_variance = max(0.0d0, (msqPart - msqRpt) / (nops * nrepeats))
        rst%operator_variance = max(0.0d0, (msqOp - msqRpt) / &
            (npart * nrepeats))
    end if
    rst%operator_by_part_variance = max(0.0d0, (msqPartOp - msqRpt) / nrepeats)
    rst%equipment_variance = msqRpt
    rst%measurement_variance = rst%equipment_variance + &
        rst%operator_variance + rst%operator_by_part_variance
    rst%total_variance = rst%measurement_variance + rst%part_variance
end function

! ------------------------------------------------------------------------------
! Control Chart Type GR&R
! - https://www.spcforexcel.com/knowledge/measurement-systems-analysis/three-methods-analyze-gage-rr-studies
! - https://www.qualitydigest.com/inside/twitter-ed/problems-gauge-rr-studies.html
pure module function control_chart_variance(x) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:,:,:) :: x
    type(process_variance) :: rst

    ! Local Variables
    integer(int32) :: i, j, k, npart, nrepeats, nops
    real(real64), allocatable, dimension(:) :: opsMeans, rptMeans, partMeans, &
        repRng
    real(real64) :: opsRange, partRange, avgRange

    ! Initialization
    npart = size(x, 1)
    nrepeats = size(x, 2)
    nops = size(x, 3)
    rst%operator_by_part_variance = 0.0d0   ! Always 0 for this analysis

    ! Quick Return
    if (npart < 2 .or. nrepeats < 2 .or. nops < 1) then
        rst%measurement_variance = 0.0d0
        rst%part_variance = 0.0d0
        rst%total_variance = 0.0d0
        rst%equipment_variance = 0.0d0
        rst%operator_variance = 0.0d0
        return
    end if

    ! Local Memory Allocation
    allocate(opsMeans(nops))
    allocate(rptMeans(nrepeats))
    allocate(partMeans(npart))
    allocate(repRng(npart * nops))

    ! Compute the mean (average) terms
    do i = 1, nops
        opsMeans(i) = mean(reshape(x(:,:,i), [npart * nrepeats]))
    end do
    do i = 1, nrepeats
        rptMeans(i) = mean(reshape(x(:,i,:), [npart * nops]))
    end do
    do i = 1, npart
        partMeans(i) = mean(reshape(x(i,:,:), [nrepeats * nops]))
    end do

    ! Compute the range terms
    opsRange = data_range(opsMeans)
    partRange = data_range(partMeans)
    j = 0
    do k = 1, nops
        do i = 1, npart
            j = j + 1
            repRng(j) = data_range(x(i,:,k))
        end do
    end do
    avgRange = data_range(repRng)

    ! Compute the variance terms
    rst%equipment_variance = (avgRange / d2(nrepeats))**2
    if (nops > 1) then
        rst%operator_variance = (opsRange / d2(nops))**2 - &
            (real(nops) / real(npart * nops * nrepeats)) * &
            rst%equipment_variance
    else
        rst%operator_variance = 0.0d0
    end if 
    rst%measurement_variance = rst%equipment_variance + rst%operator_variance
    rst%part_variance = (partRange / d2(npart)**2) - &
        (real(npart) / real(npart * nops * nrepeats)) * &
        rst%equipment_variance
contains
    ! REF: 
    ! - http://www.bessegato.com.br/UFJF/resources/table_of_control_chart_constants_old.pdf
    ! - https://r-bar.net/control-chart-constants-tables-explanations/#XmR_Sec
    pure function d2(n) result(d)
        ! Arguments
        integer(int32), intent(in) :: n
        real(real64) :: d

        ! Define the table and associated indices
        integer(int32), parameter, dimension(29) :: index = [ &
            2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, &
            17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
        real(real64), parameter, dimension(29) :: values = [ &
            1.128d0, 1.693d0, 2.059d0, 2.326d0, 2.534d0, 2.704d0, &
            2.847d0, 2.970d0, 3.078d0, 3.173d0, 3.258d0, 3.336d0, &
            3.407d0, 3.472d0, 3.532d0, 3.588d0, 3.640d0, 3.689d0, &
            3.735d0, 3.778d0, 3.819d0, 3.858d0, 3.895d0, 3.931d0, &
            3.964d0, 3.997d0, 4.027d0, 4.057d0, 4.086d0]

        ! Local Variables
        integer(int32) :: i

        ! Ensure the index is within bounds
        if (n < 2) then
            d = 0.0d0
            return
        else if (n > size(index)) then
            d = values(size(index))
            return
        end if

        ! Find the matching array index, and return the correct value
        do i = 1, size(index)
            if (n == index(i)) exit
        end do
        d = values(i)
    end function
end function

! ------------------------------------------------------------------------------
pure module function compute_grr(k, x, usl, lsl) result(rst)
    ! Arguments
    real(real64), intent(in) :: k, usl, lsl
    type(process_variance), intent(in) :: x
    type(grr_results) :: rst

    ! Process
    rst%tolerance_range = abs(usl - lsl)
    rst%pt_ratio = k * sqrt(x%measurement_variance) / rst%tolerance_range
    rst%ptv_ratio = sqrt(x%measurement_variance) / sqrt(x%total_variance)
end function

! ------------------------------------------------------------------------------
pure elemental module function discrimination_ratio(tv, mv) result(x)
    real(real64), intent(in) :: tv, mv
    real(real64) :: x
    x = sqrt(2.0d0 * (tv / mv) - 1.0d0)
end function

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
pure module function ssq_part(x, xmean) result(ssq)
    ! Arguments
    real(real64), intent(in), dimension(:,:,:) :: x
    real(real64), intent(in) :: xmean
    real(real64) :: ssq

    ! Local Variables
    integer(int32) :: i, j, k, l, nop, npart, nrep
    real(real64) :: xi
    real(real64), allocatable, dimension(:) :: y

    ! Initialization
    npart = size(x, 1)
    nrep = size(x, 2)
    nop = size(x, 3)
    allocate(y(nrep * nop))

    ! Cycle over each part, and compute the sum of the differences
    ssq = 0.0d0
    do i = 1, npart
        ! Collect the data whose mean we are to evaluate, 
        ! and then compute the mean
        l = 0
        do j = 1, nrep
            do k = 1, nop
                l = l + 1
                y(l) = x(i,j,k)
            end do
        end do
        xi = mean(y)

        ! Compute and include the difference
        ssq = ssq + (xi - xmean)**2
    end do

    ! Output
    ssq = nop * nrep * ssq
end function

! ------------------------------------------------------------------------------
pure module function ssq_operator(x, xmean) result(ssq)
    ! Arguments
    real(real64), intent(in), dimension(:,:,:) :: x
    real(real64), intent(in) :: xmean
    real(real64) :: ssq

    ! Local Variables
    integer(int32) :: i, j, k, l, nop, npart, nrep
    real(real64) :: xi
    real(real64), allocatable, dimension(:) :: y

    ! Initialization
    npart = size(x, 1)
    nrep = size(x, 2)
    nop = size(x, 3)
    allocate(y(npart * nrep))

    ! Cycle over each operator
    ssq = 0.0d0
    do k = 1, nop
        ! Collect the data whose mean we are to evaluate, 
        ! and then compute the mean
        l = 0
        do j = 1, nrep
            do i = 1, npart
                l = l + 1
                y(l) = x(i,j,k)
            end do
        end do
        xi = mean(y)

        ! Compute and include the difference
        ssq = ssq + (xi - xmean)**2
    end do

    ! Output
    ssq = npart * nrep * ssq
end function

! ------------------------------------------------------------------------------
pure module function ssq_repeat(x) result(ssq)
    ! Arguments
    real(real64), intent(in), dimension(:,:,:) :: x
    real(real64) :: ssq

    ! Local Variables
    integer(int32) :: i, j, k, nop, npart, nrep
    real(real64) :: xik
    
    ! Initialization
    npart = size(x, 1)
    nrep = size(x, 2)
    nop = size(x, 3)

    ! Cycle over each observation
    ssq = 0.0d0
    do k = 1, nop
        do i = 1, npart
            xik = mean(x(i,:,k))
            do j = 1, nrep
                ssq = ssq + (x(i,j,k) - xik)**2
            end do
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure module function ssq_total(x, xmean) result(ssq)
    ! Arguments
    real(real64), intent(in), dimension(:,:,:) :: x
    real(real64), intent(in) :: xmean
    real(real64) :: ssq

    ! Local Variables
    integer(int32) :: i, j, k, nop, npart, nrep
    
    ! Initialization
    npart = size(x, 1)
    nrep = size(x, 2)
    nop = size(x, 3)

    ! Cycle over each part, and compute the results
    ssq = 0.0d0
    do k = 1, nop
        do j = 1, nrep
            do i = 1, npart
                ssq = ssq + (x(i,j,k) - xmean)**2
            end do
        end do
    end do
end function

! ------------------------------------------------------------------------------
end submodule
