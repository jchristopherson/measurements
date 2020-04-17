! signals_tests.f90

module signals_tests
    use iso_fortran_env
    use measurements_core
    implicit none
contains
! ------------------------------------------------------------------------------
    function fourier_transform_test() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter :: npts = 1024
        real(real64), parameter :: fs = 1024.0d0
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        real(real64), parameter :: f1 = 5.0d0
        real(real64), parameter :: f2 = 50.0d0
        real(real64), parameter :: x1 = 3.0d0
        real(real64), parameter :: x2 = 0.5d0
        integer(int32), parameter :: ind1 = 6
        integer(int32), parameter :: ind2 = 51
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        integer(int32) :: i
        real(real64) :: t, dt, x(npts), delta
        complex(real64), allocatable, dimension(:) :: f

        ! Initialization
        rst = .true.
        dt = 1.0d0 / fs
        t = 0.0d0
        x(1) = 0.0d0
        do i = 1, npts
            t = t + dt
            x(i) = x1 * sin(2.0d0 * pi * f1 * t) + x2 * sin(2.0d0 * pi * f2 * t)
        end do

        ! Compute the Fourier transform of X
        f = fourier_transform(x)

        ! Check
        if (size(f) /= npts / 2 + 1) then
            rst = .false.
            print '(A)', "FOURIER_TRANSFORM FAILED."
            print '(AI0AI0A)', "Expected an array of length ", npts / 2 + 1, &
                ", but found an array of length ", size(f), "."
        end if

        delta = abs(f(ind1)) - x1
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "FOURIER_TRANSFORM FAILED."
            print *, "Expected: ", x1
            print *, "Computed: ", abs(f(ind1))
            print *, "Difference (Expected - Computed): ", delta
        end if

        delta = abs(f(ind2)) - x2
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "FOURIER_TRANSFORM FAILED."
            print *, "Expected: ", x2
            print *, "Computed: ", abs(f(ind2))
            print *, "Difference (Expected - Computed): ", delta
        end if
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
