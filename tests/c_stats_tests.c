// c_stats_tests.c

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "measurements.h"
#include "c_array.h"
#include "c_test_macros.h"

// Global Variables
const double tol = 1.0e-8;

// Function Prototypes
bool mean_test();
bool median_test_odd();
bool median_test_even();
bool variance_test();
bool std_dev_test();

bool range_test();
bool z_score_test();
bool t_score_test();
bool confidence_interval_test();
bool normal_distribution_test();
bool t_distribution_test();
bool beta_distribution_test();
bool f_distribution_test();
bool t_test_test();
bool f_test_test();

// Program
int main() {
    // Local Variables
    bool local, overall;

    // Initialization
    overall = true;

    // Introduction
    printf("\nC API MEASUREMENT STATS TESTS UNDERWAY...\n");

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

    local = range_test();
    if (!local) overall = false;

    local = z_score_test();
    if (!local) overall = false;

    local = t_score_test();
    if (!local) overall = false;

    local = confidence_interval_test();
    if (!local) overall = false;

    // local = normal_distribution_test();
    // if (!local) overall = false;

    // local = t_distribution_test();
    // if (!local) overall = false;

    // local = beta_distribution_test();
    // if (!local) overall = false;

    // local = f_distribution_test();
    // if (!local) overall = false;

    // local = t_test_test();
    // if (!local) overall = false;

    // local = f_test_test();
    // if (!local) overall = false;

    // End
    if (overall) {
        printf("C API MEASUREMENT STATS TESTS COMPLETED SUCCESSFULLY.\n\n");
    }
    else {
        printf("C API MEASUREMENT STATS TESTS FAILED.\n\n");
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



bool range_test() {
    // Local Variables
    bool rst;
    const int npts = 10000;
    double x[10000], ans, computed, delta;

    // Initialization
    rst = true;
    create_random_array(npts, x);

    // Process
    computed = c_data_range(npts, x);

    // Compute the actual solution
    ans = array_max(npts, x) - array_min(npts, x);

    // Test
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("RANGE_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // End
    return rst;
}



bool z_score_test() {
    // Local Variables
    bool rst;
    const double tol = 1.0e-8;
    const double z80 = 1.281551565545;
    const double z90 = 1.644853626951;
    const double z95 = 1.959963984540;
    const double z99 = 2.575829303549;
    const double z999 = 3.290526731492;
    double ans, computed, delta;
    int flag;

    // Initialization
    rst = true;

    // 80%
    ans = z80;
    flag = c_z_score(0.8, &computed);
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("Z_SCORE_TEST 80% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 90%
    ans = z90;
    flag = c_z_score(0.9, &computed);
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("Z_SCORE_TEST 90% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 95%
    ans = z95;
    flag = c_z_score(0.95, &computed);
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("Z_SCORE_TEST 95% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 99%
    ans = z99;
    flag = c_z_score(0.99, &computed);
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("Z_SCORE_TEST 99% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 99.9%
    ans = z999;
    flag = c_z_score(0.999, &computed);
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("Z_SCORE_TEST 99.9% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // End
    return rst;
}



bool t_score_test() {
    // Local Variables
    bool rst;
    const int npts = 21;
    const double tol = 1.0e-3;
    const double t80 = 0.860;
    const double t85 = 1.064;
    const double t90 = 1.325;
    const double t95 = 1.725;
    double ans, computed, delta;
    int flag;

    // Initialization
    rst = true;

    // 80%
    ans = t80;
    flag = c_t_score(0.8, npts, &computed);
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("T_SCORE_TEST 80% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 85%
    ans = t85;
    flag = c_t_score(0.85, npts, &computed);
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("T_SCORE_TEST 85% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 90%
    ans = t90;
    flag = c_t_score(0.9, npts, &computed);
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("T_SCORE_TEST 90% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 95%
    ans = t95;
    flag = c_t_score(0.95, npts, &computed);
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("T_SCORE_TEST 95% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // End
    return rst;
}



bool confidence_interval_test() {
    // Local Variables
    bool rst;
    const double alpha = 0.95;
    const double ans = 0.196943590756784;
    const double tol = 1.0e-8;
    const int npts = 10;
    int flag;
    double zval, delta, computed;
    double x[] = {0.1266904993777120, 0.3827711768054340, 0.1370953805850570, 
        0.5852213531153070, 0.2267533281658030, 0.0999861358308985, 
        0.5851003510284570, 0.8136628645855180, 0.7400357894369070, 
        0.978777475520868};

    // Initialization
    rst = true;

    // Process
    flag = c_z_score(alpha, &zval);
    computed = c_confidence_interval(npts, x, zval);

    // Test
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("CONFIDENCE_INTERVAL_TEST 85% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // End
    return rst;
}



bool normal_distribution_test() {
    // Local Variables
    bool rst;

    // End
    return rst;
}



bool t_distribution_test() {
    // Local Variables
    bool rst;

    // End
    return rst;
}



bool beta_distribution_test() {
    // Local Variables
    bool rst;

    // End
    return rst;
}



bool f_distribution_test() {
    // Local Variables
    bool rst;

    // End
    return rst;
}



bool t_test_test() {
    // Local Variables
    bool rst;

    // End
    return rst;
}



bool f_test_test() {
    // Local Variables
    bool rst;

    // End
    return rst;
}



