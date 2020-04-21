! measurements_filter.f90

submodule (measurements_core) measurements_filter
    use real_transform_routines, only : rfft1i, rfft1b, rfft1f
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
    ! Applies a mask to the Fourier transformed version of X such that 
    ! H = X * F, and then computes the inverse Fourier transform to arrive at
    ! the filtered signal.
    !
    ! - x: On input, the N-element signal to filter.  On output, the N-element
    !   filtered signal.
    ! - mask: The M-element filter to apply.  If N is even, M = N / 2 + 1; else,
    !   if N is odd, then M = (N + 1) / 2.
    ! - err: An error handling mechanism.  The following errors are possible:
    !   - M_OUT_OF_MEMORY_ERROR
    subroutine filter_by_mask(x, mask, err)
        ! Arguments
        real(real64), intent(inout), dimension(:) :: x
        complex(real64), intent(in), dimension(:) :: mask
        class(errors), intent(inout) :: err

        ! Local Variables
        integer(int32) :: i, n, m, lwsave, lwork, flag
        real(real64) :: nd, a, b, c, d
        real(real64), allocatable, dimension(:) :: wsave, work

        ! Initialization
        n = size(x)
        m = size(mask)
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

        ! Multiply with the mask.  Remember, complex multiplication is
        ! as follows:
        ! (a + ib) * (c + id) = (ac - bd) + i(ad + bc)
        x(1) = x(1) * real(mask(1), real64)
        if (mod(n, 2) == 0) then
            ! N is even, and the transform vector length is N / 2 + 1.  
            ! Remember, the DC and Nyquist terms are always real-valued for
            ! an even-length data set.
            do i = 1, n / 2
                a = x(2 * i - 2)
                b = x(2 * i - 1)
                c = real(mask(i), real64)
                d = aimag(mask(i))
                x(2 * i - 2) = a * c - b * d
                x(2 * i - 1) = a * d + b * c
            end do
            x(m) = x(m) * real(mask(m), real64)
        else
            ! N is odd, and the transform vector length is (N + 1) / 2
            do i = 1, n / 2
                a = x(2 * i - 2)
                b = x(2 * i - 1)
                c = real(mask(i), real64)
                d = aimag(mask(i))
                x(2 * i - 2) = a * c - b * d
                x(2 * i - 1) = a * d + b * c
            end do
        end if

        ! Compute the inverse transform to obtain the filtered signal
        call rfft1b(n, 1, x, n, wsave, lwsave, work, lwork, flag)
    end subroutine

! ------------------------------------------------------------------------------
end submodule
