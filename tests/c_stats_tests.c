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

    local = normal_distribution_test();
    if (!local) overall = false;

    local = t_distribution_test();
    if (!local) overall = false;

    local = beta_distribution_test();
    if (!local) overall = false;

    local = f_distribution_test();
    if (!local) overall = false;

    local = t_test_test();
    if (!local) overall = false;

    local = f_test_test();
    if (!local) overall = false;

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

    // Compute the actual solution the ole-school way
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
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("Z_SCORE_TEST 80% FAILED\nOutput Flag: %i\n", flag);
    }
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("Z_SCORE_TEST 80% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 90%
    ans = z90;
    flag = c_z_score(0.9, &computed);
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("Z_SCORE_TEST 90% FAILED\nOutput Flag: %i\n", flag);
    }
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("Z_SCORE_TEST 90% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 95%
    ans = z95;
    flag = c_z_score(0.95, &computed);
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("Z_SCORE_TEST 95% FAILED\nOutput Flag: %i\n", flag);
    }
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("Z_SCORE_TEST 95% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 99%
    ans = z99;
    flag = c_z_score(0.99, &computed);
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("Z_SCORE_TEST 99% FAILED\nOutput Flag: %i\n", flag);
    }
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("Z_SCORE_TEST 99% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 99.9%
    ans = z999;
    flag = c_z_score(0.999, &computed);
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("Z_SCORE_TEST 99.9% FAILED\nOutput Flag: %i\n", flag);
    }
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
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("T_SCORE_TEST 80% FAILED\nOutput Flag: %i\n", flag);
    }
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("T_SCORE_TEST 80% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 85%
    ans = t85;
    flag = c_t_score(0.85, npts, &computed);
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("T_SCORE_TEST 85% FAILED\nOutput Flag: %i\n", flag);
    }
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("T_SCORE_TEST 85% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 90%
    ans = t90;
    flag = c_t_score(0.9, npts, &computed);
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("T_SCORE_TEST 90% FAILED\nOutput Flag: %i\n", flag);
    }
    delta = ans - computed;
    if (fabs(delta) > tol) {
        rst = false;
        printf("T_SCORE_TEST 90% FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, computed, delta);
    }

    // 95%
    ans = t95;
    flag = c_t_score(0.95, npts, &computed);
    if (flag != M_NO_ERROR) {
        rst = false;
        printf("T_SCORE_TEST 95% FAILED\nOutput Flag: %i\n", flag);
    }
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
    const double tol = 1.0e-8;
    const int npts = 21;
    const double avg = 5.26404971932070e-1;
    const double sigma = 2.11034089338365e-1;
    double delta, f[21];
    int i;
    double x[] = {-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, 
        -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 
        3.5, 4.0, 4.5, 5.0};
    double ans[] = {2.306272166314290e-149, 1.229690115816890e-123, 
        2.392009083889280e-100, 1.697508606279160e-79, 4.394840953499330e-61, 
        4.151034651360500e-45, 1.430380472227950e-31, 1.798161951816050e-20, 
        8.246849202359420e-12, 1.379841874725850e-5, 8.422724695133530e-02, 
        1.875676375640870, 1.523860560727160e-1, 4.516630549640910e-5, 
        4.883890565618130e-11, 1.926634346379860e-19, 2.772775386838060e-30, 
        1.455834812310320e-43, 2.788634062219750e-59, 1.948735843727060e-77, 
        4.968169870317840e-98};
    
    // Process
    rst = true;
    c_normal_distribution(avg, sigma, npts, x, f);

    // Test
    for (i = 0; i < npts; ++i) {
        delta = ans[i] - f[i];
        if (fabs(delta) > tol) {
            rst = false;
            printf("NORMAL_DISTRIBUTION_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\nIndex: %i\n",
                ans[i], f[i], delta, i);
        }
    }

    // End
    return rst;
}



bool t_distribution_test() {
    // Local Variables
    bool rst;
    const double tol = 1.0e-8;
    const int npts = 21;
    const double dof = 20.0;
    double delta, f[21];
    int i;
    double x[] = {-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, 
        -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 
        3.5, 4.0, 4.5, 5.0};
    double ans[] = {0.0000789891062440353, 0.0002548336678335860, 
        0.0008224743001331390, 0.0026105772275963500, 
        0.0079637866461806600, 0.0226694437191449000, 
        0.0580872152473570000, 0.1286273829721460000, 
        0.2360456491267010000, 0.3458086123837420000, 
        0.3939885857114330000, 0.3458086123837420000, 
        0.2360456491267010000, 0.1286273829721460000, 
        0.0580872152473570000, 0.0226694437191449000, 
        0.0079637866461806600, 0.0026105772275963500, 
        0.0008224743001331390, 0.0002548336678335860, 
        0.0000789891062440353};
    
    // Process
    rst = true;
    c_t_distribution(dof, npts, x, f);

    // Test
    for (i = 0; i < npts; ++i) {
        delta = ans[i] - f[i];
        if (fabs(delta) > tol) {
            rst = false;
            printf("T_DISTRIBUTION_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\nIndex: %i\n",
                ans[i], f[i], delta, i);
        }
    }

    // End
    return rst;
}



bool beta_distribution_test() {
    // Local Variables
    bool rst;
    const double tol = 1.0e-8;
    const int npts = 18;
    const double a = 1.0;
    const double b = 3.0;
    int i;
    double delta, f[18];
    double x[] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 
        0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95};
    double ans[] = {2.43000, 2.16750, 1.92000, 1.68750, 1.47000, 1.26750, 
        1.08000, 0.90750, 0.75000, 0.60750, 0.48000, 0.36750, 
        0.27000, 0.18750, 0.12000, 0.06750, 0.03000, 0.00750};

    // Process
    rst = true;
    c_beta_distribution(a, b, npts, x, f);

    // Test
    for (i = 0; i < npts; ++i) {
        delta = ans[i] - f[i];
        if (fabs(delta) > tol) {
            rst = false;
            printf("BETA_DISTRIBUTION_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\nIndex: %i\n",
                ans[i], f[i], delta, i);
        }
    }

    // End
    return rst;
}



bool f_distribution_test() {
    // Local Variables
    bool rst;
    const double tol = 1.0e-8;
    const int npts = 18;
    const double d1 = 10.0;
    const double d2 = 1.0;
    int i;
    double delta, f[18];
    double x[] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 
        0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95};
    double ans[] = {0.27189774911348000, 0.40342757249598100, 
        0.46776063476011400, 0.48916613056974600, 0.48666000366211000, 
        0.47170875848615800, 0.45079130426395800, 0.42748989305298300, 
        0.40375575782592800, 0.38062550389189700, 0.35862153897298900, 
        0.33797693751338100, 0.31876293731516800, 0.30096192124722400, 
        0.28450947518162900, 0.26931867071908400, 0.25529401072011300, 
        0.24233930873721800};

    // Process
    rst = true;
    c_f_distribution(d1, d2, npts, x, f);

    // Test
    for (i = 0; i < npts; ++i) {
        delta = ans[i] - f[i];
        if (fabs(delta) > tol) {
            rst = false;
            printf("F_DISTRIBUTION_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\nIndex: %i\n",
                ans[i], f[i], delta, i);
        }
    }

    // End
    return rst;
}



bool t_test_test() {
    // Local Variables
    bool rst;
    const int n1 = 15;
    const int n2 = 12;
    const double tol = 1.e-8;
    const double ans1 = 0.6574262172637550;
    const double ans2 = 0.6621143019183570;
    const double ans3 = 0.2286303186286610;
    statistic t;
    double delta;
    double x1[] = {0.5875701638177520, 
        0.4925038626867790, 
        0.0997847738770978, 
        0.5540924574002610, 
        0.0833121626712929, 
        0.1738451549308330, 
        0.3521655274264620, 
        0.2239625528107020, 
        0.6871828620071030, 
        0.3248518075223050, 
        0.0551977898473518, 
        0.8648552295498370, 
        0.9239272586628300, 
        0.4917627939852090, 
        0.3508690262031490};
    double x2[] = {0.7557955531972870, 
        0.0482975843398515, 
        0.7609889442453010, 
        0.4898203045069780, 
        0.4382872343657070, 
        0.9676872466629530, 
        0.1167483190258670, 
        0.0399776180777329, 
        0.2528774837460510, 
        0.6824976673552180, 
        0.6602062072459940, 
        0.4015093296585650};
    double x3[] = {0.3201877837239090, 
        0.0980256595288177, 
        0.6897988918691660, 
        0.1785484851694640, 
        0.0991062800273234, 
        0.1195800744029930, 
        0.8476199670433790, 
        0.8536320559829150, 
        0.6394323340044970, 
        0.8848230532535040, 
        0.9300526849294520, 
        0.6703901525053320, 
        0.7168448453351630, 
        0.9870657922150660, 
        0.2874068518452400};

    // Equal Variance Assumption
    c_t_test(n1, x1, n2, x2, EQUAL_VARIANCE_ASSUMPTION, &t);
    delta = ans1 - t.probability;
    if (fabs(delta) > tol) {
        rst = false;
        printf("T_TEST_TEST - EQUAL - FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans1, t.probability, delta);
    }

    // Unequal Variance Assumption
    c_t_test(n1, x1, n2, x2, UNEQUAL_VARIANCE_ASSUMPTION, &t);
    delta = ans2 - t.probability;
    if (fabs(delta) > tol) {
        rst = false;
        printf("T_TEST_TEST - UNEQUAL - FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans2, t.probability, delta);
    }

    // Paired Data Assumption
    c_t_test(n1, x1, n1, x3, PAIRED_DATA_SET_ASSUMPTION, &t);
    delta = ans3 - t.probability;
    if (fabs(delta) > tol) {
        rst = false;
        printf("T_TEST_TEST - PAIRED - FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans3, t.probability, delta);
    }

    // End
    return rst;
}



bool f_test_test() {
    // Local Variables
    bool rst;
    const int n1 = 15;
    const int n2 = 15;
    const double tol = 1.e-8;
    const double ans = 0.46459183414;
    statistic f;
    double delta;
    double x1[] = {0.583876981, 
        0.948726432, 
        0.170087492, 
        0.054000359, 
        0.490074142, 
        0.914983043, 
        0.277306622, 
        0.146725403, 
        0.753100664, 
        0.503987212, 
        0.59636224, 
        0.478197551, 
        0.967097051, 
        0.863679029, 
        0.610842794};
    double x2[] = {0.21863735, 
        0.78811899, 
        0.89805827, 
        0.20362421, 
        0.20809033, 
        0.01437123, 
        0.74745013, 
        0.49003164, 
        0.62212261, 
        0.46324487, 
        0.61387136, 
        0.99010175, 
        0.28220738, 
        0.39286567, 
        0.04010292};

    // Process
    c_f_test(n1, x1, n2, x2, &f);
    delta = ans - f.probability;
    if (fabs(delta) > tol) {
        rst = false;
        printf("F_TEST_TEST FAILED\nExpected: %f\nComputed: %f\nDifference: %f\n",
            ans, f.probability, delta);
    }

    // End
    return rst;
}



