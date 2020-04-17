! measurements_spectrum.f90

submodule (measurements_core) measurements_spectrum
    use real_transform_routines, only : rfft1i, rfft1b, rfft1f

    ! Constants
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    !> @brief This type provides a mechanism allowing straight-forward signal
    !! overlapping and averaging for periodogram estimates.
    type spectrum_register
        integer(int32) :: window_size   ! Must be a power of 2
        integer(int32) :: segment_count = 0
        real(real64), allocatable, dimension(:) :: segments ! window_size
        real(real64), allocatable, dimension(:) :: buffer   ! window_size/2 + 1
        ! Workspace Arrays - FFT Operations
        real(real64), allocatable, dimension(:) :: wsave
        real(real64), allocatable, dimension(:) :: work
    contains
        procedure, public :: initialize => initialize_spectrum_register
        procedure, private :: add_segment => add_spectrum_segment
        procedure, public :: add => overlap_segment
    end type
contains
! ------------------------------------------------------------------------------
pure module function is_power_of_two(n) result(rst)
    ! Arguments
    integer(int32), intent(in) :: n
    logical :: rst

    ! Process
    rst = (n /= 0) .and. (iand(n, n - 1) == 0)
    ! This is equivalent to the C form: n != 0 && ((n & (n - 1)) == 0)
end function

! ------------------------------------------------------------------------------
pure module function next_power_of_two(x) result(n)
    ! Arguments
    integer(int32), intent(in) :: x
    integer(int32) :: n

    ! Process
    n = ceiling(log(real(x, real64)) / log(2.0d0), int32)
end function

! ------------------------------------------------------------------------------
pure module function pad_with_zeros(x) result(xp)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64), allocatable, dimension(:) :: xp

    ! Local Variables
    integer(int32) :: nx, n

    ! Initialization
    nx = size(x)

    ! Process
    if (is_power_of_two(nx)) then
        ! Simply return the array
        xp = x
    else
        ! Find the next higher power of two & pad with zeros
        n = 2**next_power_of_two(nx)
        allocate(xp(n))
        xp(1:nx) = x
        xp(nx+1:n) = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
module function fourier_transform(x, err) result(f)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    complex(real64), allocatable, dimension(:) :: f

    ! Local Variables
    integer(int32) :: i, n, lwsave, lwork, flag, nxfrm, nend
    real(real64), allocatable, dimension(:) :: s, wsave, work
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    lwork = n
    lwsave = n + int(log(real(n, real64)) / log(2.0d0), int32) + 4
    if (mod(n, 2) == 0) then
        nxfrm = n / 2 + 1
        nend = nxfrm - 1
    else
        nxfrm = (n + 1) / 2
        nend = nxfrm
    end if
    allocate(work(lwork), stat = flag)
    if (flag == 0) allocate(wsave(lwsave))
    if (flag == 0) allocate(s(lwsave))
    if (flag == 0) allocate(f(nxfrm), stat = flag)
    if (flag /= 0) then
        ! ERROR
        call errmgr%report_error("fourier_transform", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Initialize the transform
    call rfft1i(n, wsave, lwsave, flag)

    ! Compute the FFT
    s = x   ! Creating a copy so we don't overwrite X
    call rfft1f(n, 1, s, n, wsave, lwsave, work, lwork, flag)

    ! Extract the magnitude and phase as a complex value
    f(1) = cmplx(s(1), 0.0d0, real64)
    do i = 2, nend
        f(i) = cmplx(s(2*i-2), s(2*i-1), real64)
    end do

    ! On even-value length transforms, the Nyquist term is always real-valued
    if (nend /= nxfrm) f(nxfrm) = s(n)  

    ! End
    return
end function

! ------------------------------------------------------------------------------
module function periodogram(x, winfun, nfft, err) result(p)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    procedure(window_function), pointer, intent(in) :: winfun
    integer(int32), intent(in) :: nfft
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: p

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    type(spectrum_register) :: reg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Ensure nfft is a power of two
    if (.not.is_power_of_two(nfft)) then
        call errmgr%report_error("periodogram", "The window size " // &
            "parameter (NFFT) was expected to be an integer power of two.", &
            M_INVALID_INPUT_ERROR)
        return
    end if
    
    ! Construct the spectrum register routine
    call reg%initialize(nfft, errmgr)
    if (errmgr%has_error_occurred()) return

    call reg%add(x, winfun, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the periodogram
    p = reg%segments / reg%segment_count
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ******************************************************************************
! WINDOW ROUTINES
! ------------------------------------------------------------------------------
pure module function rectangular_window(j, n) result(x)
    integer(int32), intent(in) :: j, n
    real(real64) :: x
    x = 1.0d0
end function

! ------------------------------------------------------------------------------
pure module function hann_window(j, n) result(x)
    integer(int32), intent(in) :: j, n
    real(real64) :: x
    x = 0.5d0 * (1.0d0 - cos(2.0d0 * pi * j / n))
end function

! ------------------------------------------------------------------------------
pure module function hamming_window(j, n) result(x)
    integer(int32), intent(in) :: j, n
    real(real64) :: x
    x = 0.54d0 - 0.46d0 * cos(2.0d0 * pi * j / n)
end function

! ------------------------------------------------------------------------------
pure module function welch_window(j, n) result(x)
    integer(int32), intent(in) :: j, n
    real(real64) :: x
    x = 1.0d0 - ((j - 0.5d0 * n) / (0.5d0 * n))**2
end function

! ------------------------------------------------------------------------------
pure module function blackman_harris_window(j, n) result(x)
    integer(int32), intent(in) :: j, n
    real(real64) :: x
    x = 0.35875d0 - &
        0.48829d0 * cos(2.0d0 * pi * j / n) + &
        0.14128d0 * cos(4.0d0 * pi * j / n) - &
        0.01168d0 * cos(6.0d0 * pi * j / n)
end function

! ******************************************************************************
! SPECTRUM_REGISTER ROUTINES
! ------------------------------------------------------------------------------
!> @brief Initializes the spectrum_register object.
!!
!! @param[in,out] this The spectrum_register object.
!! @param[in] winsize The window size.  This must be an integer power of two.
!! @param[in,out] err The error handling object.  Possible errors are as 
!!  follows.
!!  - M_OUT_OF_MEMORY_ERROR: Insufficient memory available.
subroutine initialize_spectrum_register(this, winsize, err)
    ! Arguments
    class(spectrum_register), intent(inout) :: this
    integer(int32), intent(in) :: winsize
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: n, lwork, lwsave, flag

    ! Determine the length of the transformed data set
    n = winsize / 2 + 1

    ! Allocate space for arrays
    if (allocated(this%segments)) deallocate(this%segments)
    if (allocated(this%buffer)) deallocate(this%buffer)
    if (allocated(this%wsave)) deallocate(this%wsave)
    if (allocated(this%work)) deallocate(this%work)
    lwork = winsize
    lwsave = winsize + int(log(real(winsize, real64)) / log(2.0d0), int32) + 4
    allocate(this%segments(n), stat = flag)
    if (flag == 0) allocate(this%buffer(winsize), stat = flag)
    if (flag == 0) allocate(this%wsave(lwsave), stat = flag)
    if (flag == 0) allocate(this%work(lwork), stat = flag)
    if (flag /= 0) then
        call err%report_error("initialize_spectrum_register", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Initialize the transform
    call rfft1i(winsize, this%wsave, lwsave, flag)

    ! Populate the segment with all zeros
    this%segments = 0.0d0

    ! Additional Initialization
    this%window_size = winsize
    this%segment_count = 0
end subroutine

! ------------------------------------------------------------------------------
!> @brief Adds a window-sized data segment to the spectrum_register object.
!! The array length is not verified.
!!
!! @param[in,out] this The spectrum_register object.
!! @param[in] x A window-length input array.
!! @param[in] winfun The window function to apply.
subroutine add_spectrum_segment(this, x, winfun)
    ! Arguments
    class(spectrum_register), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x
    procedure(window_function), intent(in), pointer :: winfun

    ! Local Variables
    integer(int32) :: i, j, lwsave, lwork, n, m, flag
    real(real64) :: w, fac, sumw

    ! Initialization
    n = this%window_size
    m = n / 2 + 1
    lwsave = size(this%wsave)
    lwork = size(this%work)

    ! Fill the buffer with windowed data
    sumw = 0.0d0
    do i = 1, n
        j = i - 1
        w = winfun(j, n)
        this%buffer(i) = w * x(i)
        sumw = sumw + w
    end do
    fac = 0.5d0 * n / sumw

    ! Compute the FFT
    call rfft1f(n, 1, this%buffer, n, this%wsave, lwsave, this%work, lwork, &
        flag)

    ! Include in the register
    this%segments(1) = this%segments(1) + 0.5d0 * fac * (this%buffer(1))**2
    do i = 2, n / 2
        this%segments(i) = this%segments(i) + &
            fac * ((this%buffer(2*i-2))**2 + (this%buffer(2*i-1))**2)
    end do
    this%segments(m) = this%segments(m) + 0.5d0 * fac * (this%buffer(n))**2

    ! Increment the counter
    this%segment_count = this%segment_count + 1
end subroutine

! ------------------------------------------------------------------------------
!> @brief Overlaps window_size segments from the original data set to construct
!! the spectrum register.
!!
!! @param[in,out] this The spectrum_register object.
!! @param[in] x The signal.  If less in length than the overall window size the
!!  data set is padded with zeros to achieve the appropriate length.
!! @param[in] winfun The window function to apply.
!! @param[in,out] err The error handling object.  Possible errors are as 
!!  follows.
!!  - M_OUT_OF_MEMORY_ERROR: Insufficient memory available.
subroutine overlap_segment(this, x, winfun, err)
    ! Arguments
    class(spectrum_register), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x
    procedure(window_function), intent(in), pointer :: winfun
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: i, k, noff, nt, nk, m, flag
    real(real64) :: del
    real(real64), allocatable, dimension(:) :: y

    ! Initialization
    nt = size(x)
    m = this%window_size / 2 + 1
    nk = (nt - 1) / m
    if (nk > 1) then
        del = (nt - this%window_size) / (nk - 1.0d0)
    else
        del = 0.0d0
    end if
    allocate(y(this%window_size), stat = flag)
    if (flag /= 0) then
        call err%report_error("overlap_segment", &
            "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Process
    if (nt < this%window_size) then
        ! Pad with zeros
        y(1:nt) = x
        y(nt+1:this%window_size) = 0.0d0
        call this%add_segment(y, winfun)
    else
        ! Overlap
        do k = 0, nk - 1
            noff = int(k * del + 0.5d0, int32)
            do i = 1, this%window_size
                y(i) = x(noff + i)
            end do
            call this%add_segment(y, winfun)
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
end submodule
