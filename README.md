# measurements
A library that provides routines supporting calculations related to measurement activities.

## Dependencies
The measurements library depends upon the following libraries.
- [FERROR](https://github.com/jchristopherson/ferror)
- [LINALG](https://github.com/jchristopherson/linalg)
- [NONLIN](https://github.com/jchristopherson/nonlin)

This library also makes use of the [modern_fftpack](https://github.com/jlokimlin/modern_fftpack) library, but embeds it as a submodule.

## Peak Detection Example
The following example highlights basic usage of the peak-detection functionallity.

```fortran
program main
    use iso_fortran_env
    use measurements_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: dt = 1.0d-3
    real(real64), parameter :: threshold = 1.0d-2

    ! Variables
    real(real64) :: t(npts), x(npts)
    type(peak_info) :: pks
    integer(int32) :: i

    ! Create a waveform
    do i = 1, npts
        if (i == 1) then
            t(i) = 0.0d0
        else
            t(i) = t(i-1) + dt
        end if
        x(i) = exp(-2.0d0 * t(i)) * sin(15.0d0 * t(i))
    end do

    ! Locate the peaks
    pks = peak_detect(x, threshold)

    ! Print the peak and valley information
    print '(A)', "Peaks (indices, t, x)"
    do i = 1, size(pks%max_values)
        print '(AI0AF7.5AF7.5)', achar(9), &
            pks%max_value_indices(i), achar(9), &
            t(pks%max_value_indices(i)), achar(9), &
            pks%max_values(i)
    end do

    print '(A)', "Valleys (indices, t, x)"
    do i = 1, size(pks%min_values)
        print '(AI0AF7.5AF8.5)', achar(9), &
            pks%min_value_indices(i), achar(9), &
            t(pks%min_value_indices(i)), achar(9), &
            pks%min_values(i)
    end do
end program
```
The above program produces the following output.
```text
Peaks (indices, t, x)
        97      0.09600 0.81826
        516     0.51500 0.35404
        935     0.93400 0.15319
Valleys (indices, t, x)
        306     0.30500 -0.53823
        725     0.72400 -0.23288
```
The following plot illustrates the peak detection results.
![](images/peak_detect_example.png?raw=true)

## Interpolation Example
The following example illustrates the most basic use of the linear and spline interpolation routines.
```fortran
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
```
The library also exposes a C API.  The following code is the C equivalent to the Fortran code shown above.
```c
#include <stdio.h>
#include "measurements.h"

#define NINTERP 1000
#define NPTS 9

int main() {
    // Local Variables
    int i, flag;
    double dx;
    double x[NPTS] = { -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0 };
    double y[NPTS] = { 0.0, 0.15, 1.12, 2.36, 2.36, 1.46, 0.49, 0.06, 0.0 };
    double xi[NINTERP], yi[NINTERP], yp[NINTERP];

    // Define the interpolation points
    dx = (x[NPTS-1] - x[0]) / (NINTERP - 1.0);
    for (i = 0; i < NINTERP; ++i) xi[i] = i == 0 ? x[0] : xi[i - 1] + dx;

    // Compute the spline interpolation
    flag = c_spline(NPTS, x, y, NINTERP, xi, yi, 
        SPLINE_QUADRATIC_OVER_INTERVAL, 0.0, 
        SPLINE_QUADRATIC_OVER_INTERVAL, 0.0);
    
    // Compute the linear interpolation
    flag = c_interpolate(1, NPTS, x, y, NINTERP, xi, yp);

    // Print out the results
    for (i = 0; i < NINTERP; ++i) {
        printf("%f\t%f\t%f\n", xi[i], yi[i], yp[i]);
    }

    // End
    return 0;
}
```
The following plot illustrates the interpolation results.
![](images/interpolation_example.png?raw=true)
