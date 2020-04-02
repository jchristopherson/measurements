// c_stats_tests.c

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "measurements.h"
#include "c_array.h"

// Global Variables
const double tol = 1.0e-8;

// Macros
#define SQR(x) ((x) * (x))

// Function Prototypes
bool mean_test();
bool median_test_odd();
bool median_test_even();
bool variance_test();
bool std_dev_test();

// Program
int main() {
    // Local Variables
    bool local, overall;

    // Initialization
    overall = true;

    // Introduction
    printf("C API MEASUREMENT STATS TESTS UNDERWAY...\n");

    // Testing
    local = mean_test();
    if (!local) overall = false;

    local = median_test_even();
    if (!local) overall = false;

    local = median_test_odd();
    if (!local) overall = false;

    local = variance_test();
    if (!local) overall = false;

    local = std_dev_test();
    if (!local) overall = false;

    // End
    if (overall) {
        printf("C API MEASUREMENT STATS TESTS COMPLETED SUCCESSFULLY.\n");
    }
    else {
        printf("C API MEASUREMENT STATS TESTS FAILED.");
    }
}


bool mean_test() {
    // Local Variables
    const int npts = 10000;
    double x[10000], ans, computed, delta;
    bool rst;

    // Initialization
    rst = true;
    create_random_array(npts, x);

    // Compute the actual solution the old-school way
    ans = array_sum(npts, x) / ((double)npts);

    // Utilize the actual method
    computed = c_mean(npts, x);

    // Test
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("MEAN TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // End
    return rst;
}


bool median_test_odd() {
    // Local Variables
    const int npts = 10001;
    double x[10001], ans, computed, delta;
    bool rst;

    // Initialization
    rst = true;
    create_random_array(npts, x);
    array_sort(npts, x);

    // Compute the answer
    ans = x[npts / 2 + 1];

    // Process
    computed = c_median(npts, x);

    // Test
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("MEDIAN_TEST_ODD FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // End
    return rst;
}



bool median_test_even() {
    // Local Variables
    const int npts = 10000;
    double x[10000], x1, x2, ans, computed, delta;
    bool rst;

    // Initialization
    rst = true;
    create_random_array(npts, x);
    array_sort(npts, x);

    // Compute the answer
    x1 = x[npts / 2];
    x2 = x[npts / 2 + 1];
    ans = 0.5 * (x1 + x2);

    // Process
    computed = c_median(npts, x);

    // Test
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("MEDIAN_TEST_EVEN FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // End
    return rst;
}



bool variance_test() {
    // Local Variables
    const int npts = 10000;
    double x[10000], avg, ans, computed, delta;
    bool rst;
    int i;

    // Initialization
    rst = true;
    create_random_array(npts, x);

    // Process
    computed = c_variance(npts, x);

    // Compute the actual solution
    avg = c_mean(npts, x);
    ans = 0.0;
    for (i = 0; i < npts; ++i) {
        ans += SQR(x[i] - avg);
    }
    ans /= (npts - 1.0);

    // Test
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("VARIANCE_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // End
    return rst;
}




bool std_dev_test() {
    // Local Variables
    const int npts = 10000;
    double x[10000], ans, computed, delta;
    bool rst;

    // Initialization
    rst = true;
    create_random_array(npts, x);

    // Process
    computed = c_standard_deviation(npts, x);

    // Compute the actual solution
    ans = sqrt(c_variance(npts, x));

    // Test
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("STD_DEV_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // End
    return rst;
}
