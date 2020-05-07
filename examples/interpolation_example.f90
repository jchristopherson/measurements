! interpolation_example.f90

program main
    use iso_fortran_env
    use measurements_core
    implicit none

    ! Parameters
    integer(int32), parameter :: knotpts = 9
    integer(int32), parameter :: npts = 1000

    ! Local Variables
    integer(int32) :: i
    real(real64) :: dx, x(knotpts), y(knotpts), xi(npts), yi(npts), yp(npts)
    type(spline_interp) :: spline
    type(linear_interp) :: linear

    ! Define a data set:
    x = [-4.0d0, -3.0d0, -2.0d0, -1.0d0, 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]
    y = [0.0d0, 0.15d0, 1.12d0, 2.36d0, 2.36d0, 1.46d0, 0.49d0, 0.06d0, 0.0d0]

    ! Define the interpolation points
    dx = (maxval(x) - minval(x)) / (npts - 1.0d0)
    do i = 1, npts
        if (i == 1) then
            xi(i) = minval(x)
        else
            xi(i) = xi(i-1) + dx
        end if
    end do

    ! Compute the spline interpolation
    call spline%initialize(x, y)
    yi = spline%interpolate(xi)

    ! Compute the linear interpolation
    call linear%initialize(x, y)
    yp = linear%interpolate(xi)

    ! Print out the results
    do i = 1, npts
        print '(F8.5AF8.5AF8.5)', xi(i), achar(9), yi(i), achar(9), yp(i)
    end do
end program
