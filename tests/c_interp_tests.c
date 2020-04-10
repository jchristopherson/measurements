// c_interp_tests.c

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "measurements.h"
#include "c_array.h"
#include "c_test_macros.h"

bool llsq_mimo_test();
bool llsq_miso_test();

int main() {
    // Local Variables
    bool local, overall;

    // Initialization
    overall = true;

    // Introduction
    printf("\nC API MEASUREMENT INTERPOLATION TESTS UNDERWAY...\n");

    // Testing
    local = llsq_mimo_test();
    if (!local) overall = false;

    local = llsq_miso_test();
    if (!local) overall = false;

    // End
    if (overall) {
        printf("C API MEASUREMENT INTERPOLATION TESTS COMPLETED SUCCESSFULLY.\n\n");
    }
    else {
        printf("C API MEASUREMENT INTERPOLATION TESTS FAILED.\n\n");
    }
}

bool llsq_mimo_test() {
    // Local Variables
    bool rst;
    const double tol = 1.0e-8;
    const int dof = 2;
    const int npts = 34;
    int i, j, flag;
    double delta, a[4], x[68], y[68];
    double xt[] = {0.0, 0.38905, 0.77816, 0.97269, 1.16714, 
        1.556, 1.94484, 0.9726, -0.00001, 0.0, -0.388886, 
        -0.77775, -0.97215, -1.16654, -1.55533, -1.9441, -0.97171, 
        0.00004, 0.0, -0.00044, -0.0013, -0.0024, -0.00382, 
        -0.00528, -0.00257, 0.00015, 0.0, 0.00144, 0.00306, 
        0.00446, 0.00567, 0.00688, 0.00451, -0.00002, 0.0, 
        0.00122, 0.00259, 0.0029, 0.00314, 0.00338, 0.00356, 
        0.00477, -0.00001, 0.0, 0.00021, 0.00051, 0.00069, 
        0.00088, 0.0013, 0.00175, 0.00058, 0.00003, 0.0, 
        0.27156, 0.54329, 0.81507, 1.08682, 1.35881, 0.81553, 
        0.0001, 0.0, -0.27145, -0.54312, -0.81493, -1.0868, 
        -1.35879, -0.81548, 0.0};
    double yt[] = {0.0, 3000.0, 6000.0, 7500.0, 9000.0, 12000.0, 
        15000.0, 7500.0, 0.0, 0.0, -3000.0, -6000.0, -7500.0, 
        -9000.0, -12000.0, -15000.0, -7500.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 67.7908728, 135.5817456, 203.3726184, 
        271.1634912, 338.954364, 203.3726184, 0.0, 0.0, 
        -67.7908728, -135.5817456, -203.3726184, -271.1634912, 
        -338.954364, -203.3726184, 0.0};
    double ans[] = {7713.710356038485, -0.2154564584637752, 
        33.5206925537885, 249.4768502933202};
    transpose(npts, dof, xt, x);
    transpose(npts, dof, yt, y);
    
    // Compute the solution
    flag = c_linear_least_squares_mimo(dof, dof, npts, x, dof, y, dof, a, dof);

    // Test
    rst = true;
    for (j = 0; j < dof; ++j) {
        for (i = 0; i < dof; ++i) {
            delta = ans[INDEX(i,j,dof)] - a[INDEX(i,j,dof)];
            if (fabs(delta) > tol) {
                rst = false;
                printf("LLSQ_MIMO TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\nRow: %i\nColumn: %i\n",
                    ans[INDEX(i,j,dof)], a[INDEX(i,j,dof)], delta, i, j);
            }
        }
    }

    // End
    return rst;
}

bool llsq_miso_test() {
    // Local Variables
    bool rst;
    const int dof = 2;
    const int npts = 18;
    const double tol = 1.0e-8;
    int i, flag;
    double a[2], x[36], delta;
    double xt[] = {0.0, 0.38905, 0.77816, 0.97269, 1.16714, 1.556, 
        1.94484, 0.9726, -1.0e-05, 0.0, -0.388886, -0.77775, 
        -0.97215, -1.16654, -1.55533, -1.9441, -0.97171, 4.0e-05, 
        0.0, 0.0, 0.00122, 0.00259, 0.0029, 0.00314, 0.00338, 
        0.00356, 0.00477, -1.0e-05, 0.0, 0.00021, 0.00051, 
        0.00069, 0.00088, 0.0013, 0.00175, 0.00058, 3.0e-05};
    double y[] = {0.0, 3000., 6000., 7500., 9000., 12000., 15000., 
        7500., 0., 0., -3000., -6000., -7500., -9000., 
        -12000., -15000., -7500., 0.0};
    double ans[] = {7714.3632680523087, -860.96082867646840};
    transpose(npts, dof, xt, x);

    // Compute A
    flag = c_linear_least_squares_miso(dof, npts, x, dof, y, a);

    // Test
    for (i = 0; i < dof; ++i) {
        delta = ans[i] - a[i];
        if (fabs(delta) > tol) {
            rst = false;
            printf("LLSQ_MISO TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\nRow: %i\n",
                ans[i], a[i], delta, i);
        }
    }

    // End
    return rst;
}
