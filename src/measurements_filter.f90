! measurements_filter.f90

submodule (measurements_core) measurements_filter
    use real_transform_routines, only : rfft1i, rfft1b, rfft1f
contains
! ------------------------------------------------------------------------------
module function low_pass_filter(x, fs, cutoff, err) result(y)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in) :: fs, cutoff
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: y

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: i, m, n, flag
    real(real64) :: s
    complex(real64), allocatable, dimension(:) :: f
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    if (mod(n, 2) == 0) then
        m = n / 2 + 1
    else
        m = (n + 1) / 2
    end if

    ! Input Check - ensure cutoff is > 0 & cutoff < fs / 2
    if (cutoff <= 0.0d0 .or. cutoff >= 0.5d0 * fs) then
        ! ERROR:
        call errmgr%report_error("low_pass_filter", &
            "The cutoff frequency must be positive, and be less than " // &
            "the Nyquist frequency.", M_INVALID_INPUT_ERROR)
        return
    end if

    ! Allocate memory
    allocate(f(m), stat = flag)
    if (flag == 0) allocate(y(n), stat = flag)
    if (flag /= 0) then
        ! ERROR:
        call errmgr%report_error("low_pass_filter", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
    end if
    y = x

    ! Construct the filter function
    do i = 1, m
        s = fourier_frequency(fs, i - 1, n)
        if (s <= cutoff) then
            f(i) = cmplx(1.0d0, 0.0d0, real64)
        else
            f(i) = cmplx(0.0d0, 0.0d0, real64)
        end if
    end do

    ! Apply the filter
    call fourier_filter(y, f, errmgr)
end function

! ------------------------------------------------------------------------------
module function high_pass_filter(x, fs, cutoff, err) result(y)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in) :: fs, cutoff
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: y

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: i, m, n, flag
    real(real64) :: s
    complex(real64), allocatable, dimension(:) :: f
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    if (mod(n, 2) == 0) then
        m = n / 2 + 1
    else
        m = (n + 1) / 2
    end if

    ! Input Check - ensure cutoff is > 0 & cutoff < fs / 2
    if (cutoff <= 0.0d0 .or. cutoff >= 0.5d0 * fs) then
        ! ERROR:
        call errmgr%report_error("high_pass_filter", &
            "The cutoff frequency must be positive, and be less than " // &
            "the Nyquist frequency.", M_INVALID_INPUT_ERROR)
        return
    end if

    ! Allocate memory
    allocate(f(m), stat = flag)
    if (flag == 0) allocate(y(n), stat = flag)
    if (flag /= 0) then
        ! ERROR:
        call errmgr%report_error("high_pass_filter", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
    end if
    y = x

    ! Construct the filter function
    do i = 1, m
        s = fourier_frequency(fs, i - 1, n)
        if (s >= cutoff) then
            f(i) = cmplx(1.0d0, 0.0d0, real64)
        else
            f(i) = cmplx(0.0d0, 0.0d0, real64)
        end if
    end do

    ! Apply the filter
    call fourier_filter(y, f, errmgr)
end function

! ------------------------------------------------------------------------------
module function band_pass_filter(x, fs, cutoff1, cutoff2, err) result(y)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in) :: fs, cutoff1, cutoff2
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: y

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: i, m, n, flag
    real(real64) :: s
    complex(real64), allocatable, dimension(:) :: f
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    if (mod(n, 2) == 0) then
        m = n / 2 + 1
    else
        m = (n + 1) / 2
    end if

    ! Input Check - ensure cutoff is > 0 & cutoff < fs / 2
    if ((cutoff1 <= 0.0d0 .or. cutoff1 >= 0.5d0 * fs) .or. &
        (cutoff2 <= 0.0d0 .or. cutoff2 >= 0.5d0 * fs)) &
    then
        ! ERROR:
        call errmgr%report_error("band_pass_filter", &
            "The cutoff frequency must be positive, and be less than " // &
            "the Nyquist frequency.", M_INVALID_INPUT_ERROR)
        return
    end if

    ! Allocate memory
    allocate(f(m), stat = flag)
    if (flag == 0) allocate(y(n), stat = flag)
    if (flag /= 0) then
        ! ERROR:
        call errmgr%report_error("band_pass_filter", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
    end if
    y = x

    ! Construct the filter function
    do i = 1, m
        s = fourier_frequency(fs, i - 1, n)
        if (s >= min(cutoff1, cutoff2) .or. s <= max(cutoff1, cutoff2)) then
            f(i) = cmplx(1.0d0, 0.0d0, real64)
        else
            f(i) = cmplx(0.0d0, 0.0d0, real64)
        end if
    end do

    ! Apply the filter
    call fourier_filter(y, f, errmgr)
end function

! ------------------------------------------------------------------------------
module function band_stop_filter(x, fs, cutoff1, cutoff2, err) result(y)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in) :: fs, cutoff1, cutoff2
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: y

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: i, m, n, flag
    real(real64) :: s
    complex(real64), allocatable, dimension(:) :: f
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    if (mod(n, 2) == 0) then
        m = n / 2 + 1
    else
        m = (n + 1) / 2
    end if

    ! Input Check - ensure cutoff is > 0 & cutoff < fs / 2
    if ((cutoff1 <= 0.0d0 .or. cutoff1 >= 0.5d0 * fs) .or. &
        (cutoff2 <= 0.0d0 .or. cutoff2 >= 0.5d0 * fs)) &
    then
        ! ERROR:
        call errmgr%report_error("band_stop_filter", &
            "The cutoff frequency must be positive, and be less than " // &
            "the Nyquist frequency.", M_INVALID_INPUT_ERROR)
        return
    end if

    ! Allocate memory
    allocate(f(m), stat = flag)
    if (flag == 0) allocate(y(n), stat = flag)
    if (flag /= 0) then
        ! ERROR:
        call errmgr%report_error("band_stop_filter", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
    end if
    y = x

    ! Construct the filter function
    do i = 1, m
        s = fourier_frequency(fs, i - 1, n)
        if (s <= min(cutoff1, cutoff2) .or. s >= max(cutoff1, cutoff2)) then
            f(i) = cmplx(1.0d0, 0.0d0, real64)
        else
            f(i) = cmplx(0.0d0, 0.0d0, real64)
        end if
    end do

    ! Apply the filter
    call fourier_filter(y, f, errmgr)
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
! Applies a function to the Fourier transformed version of X such that 
! H = X * Y, and then computes the inverse Fourier transform to arrive at
! the filtered signal.
!
! - x: On input, the N-element signal to filter.  On output, the N-element
!   filtered signal.
! - y: The M-element filter to apply.  If N is even, M = N / 2 + 1; else,
!   if N is odd, then M = (N + 1) / 2.
! - err: An error handling mechanism.  The following errors are possible:
!   - M_OUT_OF_MEMORY_ERROR
subroutine fourier_filter(x, y, err)
    ! Arguments
    real(real64), intent(inout), dimension(:) :: x
    complex(real64), intent(in), dimension(:) :: y
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: i, n, m, lwsave, lwork, flag, nend
    real(real64) :: nd, a, b, c, d
    real(real64), allocatable, dimension(:) :: wsave, work

    ! Initialization
    n = size(x)
    m = size(y)
    nd = real(n, real64)
    lwsave = n + int(log(nd) / log(2.0d0), int32) + 4
    lwork = n
    allocate(wsave(lwsave), stat = flag)
    if (flag == 0) allocate(work(lwork), stat = flag)
    if (flag /= 0) then
        call err%report_error("filter_by_mask", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Initialize and compute the transform
    call rfft1i(n, wsave, lwsave, flag)
    call rfft1f(n, 1, x, n, wsave, lwsave, work, lwork, flag)

    ! Multiply with y.  Remember, complex multiplication is as follows:
    ! (a + ib) * (c + id) = (ac - bd) + i(ad + bc)
    if (mod(n, 2) == 0) then
        ! N is even, and the transform vector length is N / 2 + 1.  
        ! Remember, the DC and Nyquist terms are always real-valued for
        ! an even-length data set.
        nend = m - 1
        x(n) = x(n) * real(y(m), real64)
    else
        ! N is odd, and the transform vector length is (N + 1) / 2
        nend = m
    end if
    x(1) = x(1) * real(y(1), real64)    ! The DC term
    do i = 2, nend
        a = x(2 * i - 2)
        b = x(2 * i - 1)
        c = real(y(i), real64)
        d = aimag(y(i))
        x(2 * i - 2) = a * c - b * d
        x(2 * i - 1) = a * d + b * c
    end do

    ! Compute the inverse transform to obtain the filtered signal
    call rfft1b(n, 1, x, n, wsave, lwsave, work, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
end submodule
