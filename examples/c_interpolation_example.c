// c_interpolation_example.c

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