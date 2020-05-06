! measurements_array.f90

submodule (measurements_core) measurements_array
    use real_transform_routines
contains
! ------------------------------------------------------------------------------
module function cumulative_sum(x, err) result(y)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: y

    ! Local Variables
    integer(int32) :: i, n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Process
    n = size(x)
    allocate(y(n), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("cumulative_sum", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
        return
    end if
    y(1) = x(1)
    do i = 2, n
        y(i) = y(i-1) + x(i)
    end do
end function

! ------------------------------------------------------------------------------
module function difference(x, err) result(dx)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: dx

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: i, n, flag
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)

    ! Quick Return
    if (n <= 1) then
        allocate(dx(1))
        dx = 0.0d0
        return
    end if

    ! Process
    allocate(dx(n-1), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("difference", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
        return
    end if
    do i = 2, n
        dx(i) = x(i) - x(i-1)
    end do
end function

! ------------------------------------------------------------------------------
module function unwrap(x, cutoff, err) result(p)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in), optional :: cutoff
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: p

    ! Arguments
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: m, flag
    real(real64), allocatable, dimension(:) :: dp, dp_corr
    real(real64) :: c
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = size(x)
    if (present(cutoff)) then
        c = cutoff
    else
        c = pi
    end if

    ! Local Memory Allocation
    allocate(p(m), stat = flag)
    if (flag == 0) allocate(dp(m-1), stat = flag)
    if (flag == 0) allocate(dp_corr(m-1), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("unwrap", "Insufficient memory available.", &
            M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Process
    dp = difference(x, errmgr)
    if (errmgr%has_error_occurred()) return

    dp_corr = dp / (2.0d0 * pi)
    dp_corr = round(dp_corr, -1)

    do i = 1, m - 1
        if (abs(dp(i)) < c) dp_corr(i) = 0.0d0
    end do

    p(1) = x(1)
    p(2:n) = x(2:n) - (2.0d0 * pi) * cumulative_sum(dp_corr)
end function

! ------------------------------------------------------------------------------
module function finite_difference(x, y, err) result(dydx)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x, y
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: dydx

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: i, n, flag
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)

    ! Input Check
    if (size(y) /= n) then
        call errmgr%report_error("finite_difference", &
            "The input arrays must be the same size.", &
            M_ARRAY_SIZE_ERROR)
        return
    end if

    ! Local Memory Allocation
    allocate(dydx(n), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("finite_difference", &
            "Insufficient memory available.", &
            M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Quick Return
    if (n == 1) then
        dydx = 0.0d0
        return
    end if

    ! Forward Difference
    dydx(1) = (y(2) - y(1)) / (x(2) - x(1))

    ! Central Difference
    do i = 2, n - 1
        dydx(i) = (y(i + 1) - y(i - 1)) / (x(i + 1) - x(i - 1))
    end do

    ! Backward Difference
    dydx(n) = (y(n) - y(n - 1)) / (x(n) - x(n - 1))
end function

! ------------------------------------------------------------------------------
module function integrate(x, y, c, err) result(f)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x, y
    real(real64), intent(in), optional :: c
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: f

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: i, n, flag
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)

    ! Input Check
    if (size(y) /= n) then
        call errmgr%report_error("integrate", &
            "The input arrays must be the same size.", &
            M_ARRAY_SIZE_ERROR)
        return
    end if

    ! Local Memory Allocation
    allocate(f(n), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("integrate", &
            "Insufficient memory available.", &
            M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Process
    if (present(c)) then
        f(1) = c
    else
        f(1) = 0.0d0
    end if
    do i = 2, n
        f(i) = f(i-1) + (x(i) - x(i-1)) * y(i)
    end do
end function

! ------------------------------------------------------------------------------
module function trapz_integrate(x, y, err) result(f)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x, y
    class(errors), intent(inout), optional, target :: err
    real(real64) :: f

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: i, n
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)

    ! Input Check
    if (size(y) /= n) then
        call errmgr%report_error("trapz_integrate", &
            "The input arrays must be the same size.", &
            M_ARRAY_SIZE_ERROR)
        return
    end if

    ! Process
    f = 0.0d0
    do i = 1, n - 1
        f = f + 0.5d0 * (y(i + 1) + y(i)) * (x(i + 1) - x(i))
    end do
end function

! ------------------------------------------------------------------------------
module function remove_dc_offset(x, err) result(y)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: y

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: n, lwsave, lwork, flag
    real(real64), allocatable, dimension(:) :: wsave, work
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    lwork = n
    lwsave = n + int(log(real(n, real64)) / log(2.0d0), int32) + 4

    ! Local Memory Allocation
    allocate(y(n), stat = flag)
    if (flag == 0) allocate(work(lwork), stat = flag)
    if (flag == 0) allocate(wsave(lwsave), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("remove_dc_offset", &
            "Insufficient memory available.", &
            M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Initialize the transform
    call rfft1i(n, wsave, lwsave, flag)

    ! Compute the transform and zero out the DC term
    y = x
    call rfft1f(n, 1, y, n, wsave, lwsave, work, lwork, flag)
    y(1) = 0.0d0

    ! Compute the inverse transform to arrive back at the signal
    call rfft1b(n, 1, y, n, wsave, lwsave, work, lwork, flag)
end function

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
!> @brief Rounds a number to the required precision, but rounds 0.5 down to
!! the next lowest value.
!!
!! @param[in] x The value on which to operate.
!! @param[in] p The precision with which to round.
!! @return The result of the operation.
pure elemental function round(x, p) result(y)
    ! Arguments
    real(real64), intent(in) :: x
    integer(int32), intent(in) :: p
    real(real64) :: y

    ! Local Variables
    real(real64) :: scale, val

    ! Compute the scaling factor
    scale = (1.0d1)**(-p)

    ! Apply the scaling factor, and round accordingly
    val = x * scale + 0.49d0
    y = floor(val) / scale
end function

! ------------------------------------------------------------------------------
end submodule
