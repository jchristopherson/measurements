! measurements_spectrum.f90

submodule (measurements_core) measurements_spectrum
    use real_transform_routines, only : rfft1i, rfft1b, rfft1f
contains
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
    if (nend /= nxfrm) f(nxfrm) = s(n)

    ! End
    return
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule
