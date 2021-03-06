! measurements_interp.f90

submodule (measurements_core) measurements_interp
contains
! ******************************************************************************
! INTERP_MANAGER MEMBERS
! ------------------------------------------------------------------------------
module subroutine im_init(this, x, y, order, err)
    ! Arguments
    class(interp_manager), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x, y
    integer(int32), intent(in), optional :: order
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: i, n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = 256) :: errmsg

    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(order)) then
        this%m_order = order
    else
        ! Default to a first order (linear) interpolation
        this%m_order = 1
    end if
    this%m_savedIndex = 1
    this%m_indexCheck = 1
    n = size(x)
    if (size(y) /= n) then
        ! ERROR
        write(errmsg, '(AI0AI0A)') &
            "Expected the dependent variable array to be of length ", &
            size(x), ", but found an array of length ", size(y), "."
        call errmgr%report_error("im_init", trim(errmsg), &
            M_ARRAY_SIZE_ERROR)
        return
    end if
    if (.not.is_monotonic(x)) then
        ! ERROR: Non-monotonic Array
        call errmgr%report_error("im_init", &
            "The supplied independent data array was not monotonic.", &
            M_NONMONOTONIC_ARRAY_ERROR)
        return
    end if

    if (allocated(this%m_x)) deallocate(this%m_x)
    if (allocated(this%m_y)) deallocate(this%m_y)

    allocate(this%m_x(n), stat = flag)
    if (flag == 0) allocate(this%m_y(n), stat = flag)
    if (flag /= 0) then
        ! ERROR
        call errmgr%report_error("im_init", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Copy the data
    do i = 1, n
        this%m_x(i) = x(i)
        this%m_y(i) = y(i)
    end do
end subroutine

! ------------------------------------------------------------------------------
module function im_locate(this, pt, err) result(j)
    ! Arguments
    class(interp_manager), intent(inout) :: this
    real(real64), intent(in) :: pt
    class(errors), intent(inout), optional, target :: err
    integer :: j

    ! Local Variables
    integer(int32) :: n, m, jhi, jmid, jlo
    logical :: ascnd
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    j = 0
    n = size(this%m_x)
    m = this%m_order + 1
    ascnd = this%m_x(n) >= this%m_x(1)
    jlo = 1
    jhi = n
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Ensure data has been defined
    if (.not.allocated(this%m_x) .or. .not.allocated(this%m_y)) then
        call errmgr%report_error("im_locate", "No data has been defined.", &
            M_NO_DATA_DEFINED_ERROR)
        return
    end if

    ! Process
    do while (jhi - jlo > 1)
        jmid = (jhi + jlo) / 2
        if (pt >= this%m_x(jmid) .eqv. ascnd) then
            jlo = jmid
        else
            jhi = jmid
        end if
    end do

    ! Check to see if we should use a more efficient search approach next
    ! time
    this%m_correlated = abs(jlo - this%m_savedIndex) <= this%m_indexCheck
    this%m_savedIndex = jlo

    ! Output
    ! j = max(1, min(n + 1 - m, jlo - (m - 1) / 2))
    if (pt == this%m_x(1)) then
        j = 1
    else if (pt == this%m_x(n)) then
        j = n - 1
    else
        j = jlo
    end if
end function

! ------------------------------------------------------------------------------
module function im_hunt(this, pt, err) result(j)
    ! Arguments
    class(interp_manager), intent(inout) :: this
    real(real64), intent(in) :: pt
    class(errors), intent(inout), optional, target :: err
    integer(int32) :: j

    ! Local Variables
    integer(int32) :: jlo, jmid, jhi, inc, n, m
    logical :: ascnd
    class(errors), pointer :: errmgr
    type(errors), target :: deferr


    ! Initialization
    j = 0
    n = size(this%m_x)
    m = this%m_order + 1
    jlo = this%m_savedIndex
    inc = 1
    ascnd = this%m_x(n) > this%m_x(1)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Ensure data has been defined
    if (.not.allocated(this%m_x) .or. .not.allocated(this%m_y)) then
        call errmgr%report_error("im_hunt", "No data has been defined.", &
            M_NO_DATA_DEFINED_ERROR)
        return
    end if

    ! Process
    if (jlo < 1 .or. jlo > n) then
        jlo = 1
        jhi = n
    else
        if (pt >= this%m_x(jlo) .eqv. ascnd) then
            do
                jhi = jlo + inc
                if (jhi >= n) then
                    jhi = n
                    exit
                else if (pt < this%m_x(jhi) .eqv. ascnd) then
                    exit
                else
                    jlo = jhi
                    inc = inc + inc
                end if
            end do
        else
            jhi = jlo
            do
                jlo = jlo - inc
                if (jlo <= 1) then
                    jlo = 1
                    exit
                else if (pt >= this%m_x(jlo) .eqv. ascnd) then
                    exit
                else
                    jhi = jlo
                    inc = inc + inc
                end if
            end do
        end if
    end if

    ! The hunt is done, so begin the final bisection phase
    do while (jhi - jlo > 1)
        jmid = (jhi + jlo) / 2
        if (pt >= this%m_x(jmid) .eqv. ascnd) then
            jlo = jmid
        else
            jhi = jmid
        end if
    end do

    ! Check to see if we should hunt or locate the next time around
    this%m_correlated = abs(jlo - this%m_savedIndex) <= this%m_indexCheck
    this%m_savedIndex = jlo

    ! Output
    ! j = max(1, min(n + 1 - m, jlo - (m - 1) / 2))
    if (pt == this%m_x(1)) then
        j = 1
    else if (pt == this%m_x(n)) then
        j = n - 1
    else
        j = jlo
    end if
end function

! ------------------------------------------------------------------------------
module function im_perform(this, pt, err) result(yy)
    ! Arguments
    class(interp_manager), intent(inout) :: this
    real(real64), intent(in) :: pt
    class(errors), intent(inout), optional, target :: err
    real(real64) :: yy

    ! Local Variables
    integer(int32) :: jlo

    ! Process
    if (this%m_correlated) then
        jlo = this%hunt(pt, err)
    else
        jlo = this%locate(pt, err)
    end if
    yy = this%raw_interp(jlo, pt)
end function

! ------------------------------------------------------------------------------
module function im_perform_array(this, pts, err) result(yy)
    ! Arguments
    class(interp_manager), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: pts
    class(errors), intent(inout), optional, target :: err
    real(real64), dimension(size(pts)) :: yy

    ! Local Variables
    integer(int32) :: i, jlo

    ! Process
    do i = 1, size(pts)
        if (this%m_correlated) then
            jlo = this%hunt(pts(i), err)
        else
            jlo = this%locate(pts(i), err)
        end if
        yy(i) = this%raw_interp(jlo, pts(i))
    end do
end function

! ------------------------------------------------------------------------------
pure module function im_get_num_pts(this) result(n)
    class(interp_manager), intent(in) :: this
    integer(int32) :: n
    n = 0
    if (allocated(this%m_x) .and. allocated(this%m_y)) n = size(this%m_x)
end function

! ------------------------------------------------------------------------------
pure module function im_get_x(this, ind) result(x)
    class(interp_manager), intent(in) :: this
    integer(int32), intent(in) :: ind
    real(real64) :: x
    x = 0.0d0
    if (allocated(this%m_x)) x = this%m_x(ind)
end function

! ------------------------------------------------------------------------------
pure module function im_get_y(this, ind) result(y)
    class(interp_manager), intent(in) :: this
    integer(int32), intent(in) :: ind
    real(real64) :: y
    y = 0.0d0
    if (allocated(this%m_y)) y = this%m_y(ind)
end function


! ******************************************************************************
! LINEAR_INTERP MEMBERS
! ------------------------------------------------------------------------------
!> @brief Performs the actual linear interpolation.
!!
!! @param[in,out] this The linear_interp_mgr instance.
!! @param[in] jlo The array index below which @p pt is found in x.
!! @param[in] pt The independent variable value to interpolate.
!!
!! @return The interpolated value.
module function li_raw_interp(this, jlo, pt) result(yy)
    ! Arguments
    class(linear_interp), intent(inout) :: this
    integer(int32), intent(in) :: jlo
    real(real64), intent(in) :: pt
    real(real64) :: yy

    ! Process
    if (this%m_x(jlo) == this%m_x(jlo+1)) then
        yy = this%m_y(jlo)
    else
        yy = this%m_y(jlo) + ((pt - this%m_x(jlo)) / &
            (this%m_x(jlo+1) - this%m_x(jlo))) * &
            (this%m_y(jlo+1) - this%m_y(jlo))
    end if
end function

! ******************************************************************************
! POLYNOMIAL_INTERP MEMBERS
! ------------------------------------------------------------------------------
module subroutine pi_init(this, x, y, order, err)
    ! Arguments
    class(polynomial_interp), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x, y
    integer(int32), intent(in), optional :: order
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: m, flag, odr
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(order)) then
        odr = order
    else
        odr = 1
    end if

    ! Input Checking
    if (odr < 1) then
        call errmgr%report_error("pi_init", &
            "A polynomial order greater than or equal to 1 must " // &
            "be specified.", M_INVALID_INPUT_ERROR)
        return
    end if

    ! Memory Allocation
    call im_init(this, x, y, odr, err)
    if (allocated(this%m_c)) deallocate(this%m_c)
    if (allocated(this%m_d)) deallocate(this%m_d)
    m = odr + 1
    allocate(this%m_c(m), stat = flag)
    if (flag == 0) allocate(this%m_d(m), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("pi_init", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
        return
    end if
end subroutine

! ------------------------------------------------------------------------------
module function pi_raw_interp(this, jlo, pt) result(yy)
    ! Arguments
    class(polynomial_interp), intent(inout) :: this
    integer(int32), intent(in) :: jlo
    real(real64), intent(in) :: pt
    real(real64) :: yy

    ! Local Variables
    integer(int32) :: n, nn
    real(real64) :: dy

    ! Process
    nn = jlo + this%m_order
    n = min(size(this%m_x), nn)
    call polint(this%m_x(jlo:n), this%m_y(jlo:n), pt, yy, dy)
end function

! ******************************************************************************
! SPLINE_INTERP MEMBERS
! ------------------------------------------------------------------------------
module subroutine penta_solve(a1, a2, a3, a4, a5, b, x)
    ! Arguments
    real(real64), intent(in), dimension(:) :: a1, a5
    real(real64), intent(inout), dimension(:) :: a2, a3, a4, b
    real(real64), intent(out), dimension(:) :: x

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: xmult

    ! Initialization
    n = size(a1)

    ! Process
    do i = 2, n - 1
        xmult = a2(i) / a3(i - 1)
        a3(i) = a3(i) - xmult * a4(i - 1)
        a4(i) = a4(i) - xmult * a5(i - 1)
        b(i) = b(i) - xmult * b(i - 1)
        xmult = a1(i + 1) - xmult * a4(i - 1)
        a2(i + 1) = a2(i + 1) - xmult * a4(i - 1)
        a3(i + 1) = a3(i + 1) - xmult * a5(i - 1)
        b(i + 1) = b(i + 1) - xmult * b(i - 1)
    end do

    xmult = a2(n) / a3(n - 1)
    a3(n) = a3(n) - xmult * a4(n - 1)
    x(n) = (b(n) - xmult * b(n - 1)) / a3(n)
    x(n - 1) = (b(n - 1) - a4(n - 1) * x(n)) / a3(n - 1)
    do i = n - 2, 1, -1
        x(i) = (b(i) - a4(i) * x(i + 1) - a5(i) * x(i + 2)) / a3(i)
    end do
end subroutine

! ------------------------------------------------------------------------------
module function si_raw_interp(this, jlo, pt) result(yy)
    ! Arguments
    class(spline_interp), intent(inout) :: this
    integer(int32), intent(in) :: jlo
    real(real64), intent(in) :: pt
    real(real64) :: yy

    ! Parameters
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: three = 3.0d0
    real(real64), parameter :: six = 6.0d0

    ! Local Variables
    integer(int32) :: right
    real(real64) :: dt, h

    ! Initialization
    right = jlo + 1
    dt = pt - this%m_x(jlo)
    h = this%m_x(right) - this%m_x(jlo)

    ! Process
    yy = this%m_y(jlo) + dt * ((this%m_y(right) - this%m_y(jlo)) / h - &
        (this%m_ypp(right) / six + this%m_ypp(jlo) / three) * h + &
        dt * (half * this%m_ypp(jlo) + &
        dt * ((this%m_ypp(right) - this%m_ypp(jlo)) / (six * h))))
end function

! ------------------------------------------------------------------------------
module subroutine si_second_deriv(this, ibcbeg, ybcbeg, ibcend, ybcend, err)
    ! Arguments
    class(spline_interp), intent(inout) :: this
    integer(int32), intent(in) :: ibcbeg, ibcend
    real(real64), intent(in) :: ybcbeg, ybcend
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: three = 3.0d0
    real(real64), parameter :: six = 6.0d0

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: i, n, flag
    real(real64), allocatable, dimension(:) :: a1, a2, a3, a4, a5, b

    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(this%m_x)

    ! Allocate Memory
    if (allocated(this%m_ypp)) deallocate(this%m_ypp)
    allocate(this%m_ypp(n), stat = flag)
    if (flag == 0) allocate(a1(n), stat = flag)
    if (flag == 0) allocate(a2(n), stat = flag)
    if (flag == 0) allocate(a3(n), stat = flag)
    if (flag == 0) allocate(a4(n), stat = flag)
    if (flag == 0) allocate(a5(n), stat = flag)
    if (flag == 0) allocate(b(n), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("si_second_deriv", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Zero out the matrix
    do i = 1, n
        a1(i) = zero
        a2(i) = zero
        a3(i) = zero
        a4(i) = zero
        a5(i) = zero
    end do

    ! Set the first equation
    select case (ibcbeg)
    case (SPLINE_QUADRATIC_OVER_INTERVAL)
        b(1) = zero
        a3(1) = one
        a4(1) = -one
    case (SPLINE_KNOWN_FIRST_DERIVATIVE)
        b(1) = (this%m_y(2) - this%m_y(1)) / &
            (this%m_x(2) - this%m_x(1)) - ybcbeg
        a3(1) = (this%m_x(2) - this%m_x(1)) / three
        a4(1) = (this%m_x(2) - this%m_x(1)) / six
    case (SPLINE_KNOWN_SECOND_DERIVATIVE)
        b(1) = ybcbeg
        a3(1) = one
        a4(1) = zero
    case (SPLINE_CONTINUOUS_THIRD_DERIVATIVE)
        b(1) = zero
        a3(1) = this%m_x(2) - this%m_x(3)
        a4(1) = this%m_x(3) - this%m_x(1)
        a5(1) = this%m_x(1) - this%m_x(2)
    case default
        b(1) = zero
        a3(1) = one
        a4(1) = one
    end select

    ! Set the intermediate equations
    do i = 2, n - 1
        b(i) = (this%m_y(i+1) - this%m_y(i)) / &
            (this%m_x(i+1) - this%m_x(i)) - &
            (this%m_y(i) - this%m_y(i-1)) / (this%m_x(i) - this%m_x(i-1))
        a2(i) = (this%m_x(i+1) - this%m_x(i)) / six
        a3(i) = (this%m_x(i+1) - this%m_x(i-1)) / three
        a4(i) = (this%m_x(i) - this%m_x(i-1)) / six
    end do

    ! Set the last equation
    select case (ibcend)
    case (SPLINE_QUADRATIC_OVER_INTERVAL)
        b(n) = zero
        a2(n) = -one
        a3(n) = one
    case (SPLINE_KNOWN_FIRST_DERIVATIVE)
        b(n) = ybcend - (this%m_y(n) - this%m_y(n-1)) / &
            (this%m_x(n) - this%m_x(n-1))
        a2(n) = (this%m_x(n) - this%m_x(n-1)) / six
        a3(n) = (this%m_x(n) - this%m_x(n-1)) / three
    case (SPLINE_KNOWN_SECOND_DERIVATIVE)
        b(n) = ybcend
        a2(n) = zero
        a3(n) = one
    case (SPLINE_CONTINUOUS_THIRD_DERIVATIVE)
        b(n) = zero
        a1(n) = this%m_x(n-1) - this%m_x(n)
        a2(n) = this%m_x(n) - this%m_x(n-2)
        a3(n) = this%m_x(n-2) - this%m_x(n-1)
    case default
        b(n) = zero
        a2(n) = -one
        a3(n) = one
    end select

    ! Define the second derivative
    if (n == 2 .and. ibcbeg == SPLINE_QUADRATIC_OVER_INTERVAL .and. &
            ibcend == SPLINE_QUADRATIC_OVER_INTERVAL) then
        ! Deal with the special case of N = 2, and IBCBEG = IBCEND = 0
        this%m_ypp(1) = zero
        this%m_ypp(2) = zero
    else
        ! Solve the linear system
        call penta_solve(a1, a2, a3, a4, a5, b, this%m_ypp)
    end if
end subroutine

! ------------------------------------------------------------------------------
module subroutine si_init_1(this, x, y, order, err)
    ! Arguments
    class(spline_interp), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x, y
    integer(int32), intent(in), optional :: order
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: dummy

    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(order)) dummy = order ! Avoids complaining by the compiler

    ! Initialize the base object
    call im_init(this, x, y, 3, err)

    ! Evaluate the second derivatives
    call this%compute_diff2(SPLINE_QUADRATIC_OVER_INTERVAL, zero, &
        SPLINE_QUADRATIC_OVER_INTERVAL, zero, errmgr)
end subroutine

! ------------------------------------------------------------------------------
module subroutine si_init_2(this, x, y, ibcbeg, ybcbeg, ibcend, ybcend, err)
    ! Arguments
    class(spline_interp), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x, y
    integer(int32), intent(in), optional :: ibcbeg, ibcend
    real(real64), intent(in), optional :: ybcbeg, ybcend
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: ibeg, iend
    real(real64) :: ybeg, yend
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    ibeg = SPLINE_QUADRATIC_OVER_INTERVAL
    iend = SPLINE_QUADRATIC_OVER_INTERVAL
    ybeg = zero
    yend = zero
    if (present(ibcbeg)) ibeg = ibcbeg
    if (present(ybcbeg)) ybeg = ybcbeg
    if (present(ibcend)) iend = ibcend
    if (present(ybcend)) yend = ybcend

    ! Input Check
    if (ibeg /= SPLINE_CONTINUOUS_THIRD_DERIVATIVE .and. &
        ibeg /= SPLINE_KNOWN_SECOND_DERIVATIVE .and. &
        ibeg /= SPLINE_KNOWN_FIRST_DERIVATIVE .and. &
        ibeg /= SPLINE_QUADRATIC_OVER_INTERVAL) &
            ibeg = SPLINE_QUADRATIC_OVER_INTERVAL
    if (iend /= SPLINE_CONTINUOUS_THIRD_DERIVATIVE .and. &
        iend /= SPLINE_KNOWN_SECOND_DERIVATIVE .and. &
        iend /= SPLINE_KNOWN_FIRST_DERIVATIVE .and. &
        iend /= SPLINE_QUADRATIC_OVER_INTERVAL) &
            iend = SPLINE_QUADRATIC_OVER_INTERVAL

    ! Initialize the base object
    call im_init(this, x, y, 3, err)

    ! Evaluate the second derivatives
    call this%compute_diff2(ibeg, ybeg, iend, yend, errmgr)
end subroutine

! ------------------------------------------------------------------------------
module function si_diff1(this, pt, err) result(yy)
    ! Arguments
    class(spline_interp), intent(inout) :: this
    real(real64), intent(in) :: pt
    class(errors), intent(inout), optional, target :: err
    real(real64) :: yy

    ! Parameters
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: three = 3.0d0
    real(real64), parameter :: six = 6.0d0

    ! Local Variables
    integer(int32) :: jlo,right
    real(real64) :: dt, h

    ! Process
    if (this%m_correlated) then
        jlo = this%hunt(pt, err)
    else
        jlo = this%locate(pt, err)
    end if
    right = jlo + 1
    dt = pt - this%m_x(jlo)
    h = this%m_x(right) - this%m_x(jlo)
    yy = (this%m_y(right) - this%m_y(jlo)) / h - &
        (this%m_ypp(right) / six + this%m_ypp(jlo) / three) * h + &
        dt * (this%m_ypp(jlo) + &
        dt * (half * (this%m_ypp(right) - this%m_ypp(jlo)) / h))
end function

! ------------------------------------------------------------------------------
module function si_diff1_array(this, pts, err) result(yy)
    ! Arguments
    class(spline_interp), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: pts
    class(errors), intent(inout), optional, target :: err
    real(real64), dimension(size(pts)) :: yy

    ! Parameters
    real(real64), parameter :: half = 0.5d0
    real(real64), parameter :: three = 3.0d0
    real(real64), parameter :: six = 6.0d0

    ! Local Variables
    integer(int32) :: i, jlo,right
    real(real64) :: dt, h

    ! Process
    do i = 1, size(pts)
        if (this%m_correlated) then
            jlo = this%hunt(pts(i), err)
        else
            jlo = this%locate(pts(i), err)
        end if
        right = jlo + 1
        dt = pts(i) - this%m_x(jlo)
        h = this%m_x(right) - this%m_x(jlo)
        yy(i) = (this%m_y(right) - this%m_y(jlo)) / h - &
            (this%m_ypp(right) / six + this%m_ypp(jlo) / three) * h + &
            dt * (this%m_ypp(jlo) + &
            dt * (half * (this%m_ypp(right) - this%m_ypp(jlo)) / h))
    end do
end function

! ------------------------------------------------------------------------------
module function si_diff2(this, pt, err) result(yy)
    ! Arguments
    class(spline_interp), intent(inout) :: this
    real(real64), intent(in) :: pt
    class(errors), intent(inout), optional, target :: err
    real(real64) :: yy

    ! Local Variables
    integer(int32) :: jlo,right
    real(real64) :: dt, h

    ! Process
    if (this%m_correlated) then
        jlo = this%hunt(pt, err)
    else
        jlo = this%locate(pt, err)
    end if
    right = jlo + 1
    dt = pt - this%m_x(jlo)
    h = this%m_x(right) - this%m_x(jlo)
    yy = this%m_ypp(jlo) + dt * (this%m_ypp(right) - this%m_ypp(jlo)) / h
end function

! ------------------------------------------------------------------------------
module function si_diff2_array(this, pts, err) result(yy)
    ! Arguments
    class(spline_interp), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: pts
    class(errors), intent(inout), optional, target :: err
    real(real64), dimension(size(pts)) :: yy

    ! Local Variables
    integer(int32) :: i, jlo,right
    real(real64) :: dt, h

    ! Process
    do i = 1, size(pts)
        if (this%m_correlated) then
            jlo = this%hunt(pts(i), err)
        else
            jlo = this%locate(pts(i), err)
        end if
        right = jlo + 1
        dt = pts(i) - this%m_x(jlo)
        h = this%m_x(right) - this%m_x(jlo)
        yy(i) = this%m_ypp(jlo) + &
            dt * (this%m_ypp(right) - this%m_ypp(jlo)) / h
    end do
end function

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
subroutine polint(xa, ya, x, y, dy)
    ! Arguments
    real(real64), intent(in), dimension(:) :: xa, ya
    real(real64), intent(in) :: x
    real(real64), intent(out) :: y, dy

    ! Local Variables
    integer(int32) :: i, m, ns, n
    real(real64) :: den, dif, dift, ho, hp, w
    real(real64), allocatable, dimension(:) :: c, d

    ! Initialization
    n = size(xa)
    ns = 1
    dif = abs(x - xa(1))
    allocate(c(n))
    allocate(d(n))

    ! Find the nearest index
    do i = 1, n
        dift = abs(x - xa(i))
        if (dift < dif) then
            ns = i
            dif = dift
        end if
        c(i) = ya(i)
        d(i) = ya(i)
    end do
    
    ! The initial approximation
    y = ya(ns)

    ! The update procedure
    ns = ns - 1
    do m = 1, n - 1
        do i = 1, n - m
            ho = xa(i) - x
            hp = xa(i+m) - x
            w = c(i+1) - d(i)
            den = (ho - hp)
            den = w / den
            d(i) = hp * den
            c(i) = ho * den
        end do
        if (2 * ns < n - m) then
            dy = c(ns + 1)
        else
            dy = d(ns)
            ns = ns - 1
        end if
        y = y + dy
    end do
end subroutine

! ------------------------------------------------------------------------------
end submodule
