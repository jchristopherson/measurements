// c_interp_tests.c

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "measurements.h"
#include "c_array.h"
#include "c_test_macros.h"

bool poly_interp_test();
bool spline_interp_test();
bool llsq_mimo_test();
bool llsq_miso_test();
bool peak_detect_test();

int main() {
    // Local Variables
    bool local, overall;

    // Initialization
    overall = true;

    // Introduction
    printf("\nC API MEASUREMENT INTERPOLATION TESTS UNDERWAY...\n");

    // Testing
    local = poly_interp_test();
    if (!local) overall = false;

    local = spline_interp_test();
    if (!local) overall = false;

    local = llsq_mimo_test();
    if (!local) overall = false;

    local = llsq_miso_test();
    if (!local) overall = false;

    local = peak_detect_test();
    if (!local) overall = false;

    // End
    if (overall) {
        printf("C API MEASUREMENT INTERPOLATION TESTS COMPLETED SUCCESSFULLY.\n\n");
    }
    else {
        printf("C API MEASUREMENT INTERPOLATION TESTS FAILED.\n\n");
    }
}


bool poly_interp_test() {
    // Local Variables
    bool rst;
    const int npts = 21;
    const int ni = 15;
    const double tol = 1.0e-8;
    int i, flag;
    double yi[15], delta;
    double x[] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 
        0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 
        0.9, 0.95, 1.0};
    double y[] = {0.0, 0.30901699437494700, 0.58778525229247300, 
        0.80901699437494700, 0.95105651629515400, 1.0, 
        0.95105651629515400, 0.80901699437494700, 0.58778525229247300, 
        0.30901699437494800, 0.0, -0.30901699437494700, 
        -0.58778525229247300, -0.80901699437494700, 
        -0.95105651629515400, -1.0, -0.95105651629515300, 
        -0.80901699437494700, -0.58778525229247200, 
        -0.30901699437494600, 0.0};
    double xi[] = {0.0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 
        0.6, 0.675, 0.75, 0.825, 0.9, 0.975, 1.0};
    double ans[] = {0.0, 0.4484011230, 0.8090169940, 0.9755282580, 
        0.9510565160, 0.6984011230, 0.3090169940, -0.1545084970, 
        -0.5877852520, -0.8800367550, -1.0, -0.88003675500, 
        -0.58778525200, -0.154508496999999, 0.0};

    // Process
    rst = true;
    flag = c_interpolate(1, npts, x, y, ni, xi, yi);
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("POLY_INTERP_TEST FAILED\nOutput Flag: %i\n", flag);
    }

    // Test
    for (i = 0; i < ni; ++i) {
        delta = ans[i] - yi[i];
        if (fabs(delta) > tol) {
            rst = false;
            printf("POLY_INTERP_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\nIndex: %i\n",
                ans[i], yi[i], delta, i);
        }
    }

    // End
    return rst;
}



bool spline_interp_test() {
    // Local Variables
    bool rst;
    const int npts = 21;
    const int ni = 15;
    const double tol = 2.0e-3;
    int i, flag;
    double yi[15], delta;
    double x[] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 
        0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 
        0.9, 0.95, 1.0};
    double y[] = {0.0, 0.30901699437494700, 0.58778525229247300, 
        0.80901699437494700, 0.95105651629515400, 1.0, 
        0.95105651629515400, 0.80901699437494700, 0.58778525229247300, 
        0.30901699437494800, 0.0, -0.30901699437494700, 
        -0.58778525229247300, -0.80901699437494700, 
        -0.95105651629515400, -1.0, -0.95105651629515300, 
        -0.80901699437494700, -0.58778525229247200, 
        -0.30901699437494600, 0.0};
    double xi[] = {0.0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 
        0.6, 0.675, 0.75, 0.825, 0.9, 0.975, 1.0};
    double ans[] = {0.0, 0.45395574322220, 0.80901699400000, 0.98766310413990,
        0.95105651600000, 0.70708838900300, 0.30901699400000,
        -0.15643039704573, -0.58778525200000, -0.89098339112519, -1.0,
        -0.890981711444173, -0.587785252000000, -0.156516060777801, 0.0};

    // Process
    rst = true;
    flag = c_spline(npts, x, y, ni, xi, yi, SPLINE_QUADRATIC_OVER_INTERVAL, 
        0.0, SPLINE_QUADRATIC_OVER_INTERVAL, 0.0);
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("SPLINE_INTERP_TEST FAILED\nOutput Flag: %i\n", flag);
    }

    // Test
    for (i = 0; i < ni; ++i) {
        delta = ans[i] - yi[i];
        if (fabs(delta) > tol) {
            rst = false;
            printf("SPLINE_INTERP_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\nIndex: %i\n",
                ans[i], yi[i], delta, i);
        }
    }

    // End
    return rst;
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
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("LLSQ_MIMO_TEST FAILED\nOutput Flag: %i\n", flag);
    }

    // Test
    rst = true;
    for (j = 0; j < dof; ++j) {
        for (i = 0; i < dof; ++i) {
            delta = ans[INDEX(i,j,dof)] - a[INDEX(i,j,dof)];
            if (fabs(delta) > tol) {
                rst = false;
                printf("LLSQ_MIMO_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\nRow: %i\nColumn: %i\n",
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
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("LLSQ_MISO_TEST FAILED\nOutput Flag: %i\n", flag);
    }

    // Test
    for (i = 0; i < dof; ++i) {
        delta = ans[i] - a[i];
        if (fabs(delta) > tol) {
            rst = false;
            printf("LLSQ_MISO_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\nRow: %i\n",
                ans[i], a[i], delta, i);
        }
    }

    // End
    return rst;
}



bool peak_detect_test() {
    // Local Variables
    bool rst;
    const double dx = 0.01;
    const int npts = 100;
    const int nMax = 2;
    const int nMin = 2;
    const int maxInd1 = 12;
    const int maxInd2 = 62;
    const int minInd1 = 37;
    const int minInd2 = 87;
    const int bufferSize = 100;
    double pi, x, y[100];
    int i, mxBuffer[100], mnBuffer[100], nmx, nmn, flag;

    // Initialization
    rst = true;
    pi = 2.0 * acos(0.0);
    x = y[0] = 0.0;
    for (i = 1; i < npts; ++i) {
        x += dx;
        y[i] = exp(-0.7 * x) * sin(4.0 * pi * x);
    }

    // Locate the peaks
    flag = c_peak_detect(npts, y, 0.1, bufferSize, mxBuffer, &nmx, bufferSize, 
        mnBuffer, &nmn);
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("PEAK_DETECT_TEST FAILED\nOutput Flag: %i\n", flag);
    }

    // Locate the peaks
    if (nmx != nMax) {
        rst = false;
        printf("PEAK_DETECT_TEST FAILED.\n");
        printf("Expected to find %i peak values, but found %i.\n", nMax, nmx);
    }
    if (nmn != nMin) {
        rst = false;
        printf("PEAK_DETECT_TEST FAILED.\n");
        printf("Expected to find %i valley values, but found %i.\n", nMin, nmn);
    }

    // Ensure additional testing is OK without overstepping our bounds
    if (nmn < nMin || nmx < nMax) return rst;

    // Check the location of each peak
    if (mxBuffer[0] != maxInd1) {
        rst = false;
        printf("PEAK_DETECT_TEST FAILED.\n");
        printf("Expected a peak at %i, but found %i.\n", maxInd1, mxBuffer[0]);
    }
    if (mxBuffer[1] != maxInd2) {
        rst = false;
        printf("PEAK_DETECT_TEST FAILED.\n");
        printf("Expected a peak at %i, but found %i.\n", maxInd2, mxBuffer[1]);
    }

    if (mnBuffer[0] != minInd1) {
        rst = false;
        printf("PEAK_DETECT_TEST FAILED.\n");
        printf("Expected a valley at %i, but found %i.\n", minInd1, mnBuffer[0]);
    }
    if (mnBuffer[1] != minInd2) {
        rst = false;
        printf("PEAK_DETECT_TEST FAILED.\n");
        printf("Expected a valley at %i, but found %i.\n", minInd2, mnBuffer[1]);
    }

    // End
    return rst;
}

