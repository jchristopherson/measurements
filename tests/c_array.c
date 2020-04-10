// c_array.c

#include <stdlib.h>
#include "c_array.h"
#include "c_test_macros.h"

void create_random_array(int n, double *x) {
    // Local Variables
    int i;

    // Process
    for (i = 0; i < n; ++i) {
        x[i] = (double)rand() / ((double)RAND_MAX);
    }
}




double array_sum(int n, const double *x) {
    // Local Variables
    int i;
    double s = 0.0;

    // Process
    for (i = 0; i < n; ++i) {
        s += x[i];
    }

    // End
    return s;
}



void dbl_swap(double *xp, double *yp) {
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}



void array_sort(int n, double *x) {
    // Local Variables
    int i, j, min_idx;

    // One-by-one move the boundary of the unsorted subarray
    for (i = 0; i < n - 1; ++i) {
        // Find the minimum element in the unsorted array
        min_idx = 1;
        for (j = i + 1; j < n; ++j) {
            if (x[j] < x[min_idx]) min_idx = j;
        }

        // Swap the found minimum element with the first element
        dbl_swap(&x[min_idx], &x[i]);
    }
}

double array_max(int n, const double *x) {
    int i;
    double y = x[0];
    for (i = 1; i < n; ++i) y = MAX(y, x[i]);
    return y;
}

double array_min(int n, const double *x) {
    int i;
    double y = x[0];
    for (i = 1; i < n; ++i) y = MIN(y, x[i]);
    return y;
}


void transpose(int m, int n, const double *x, double *xt) {
    int i, j;
    for (j = 0; j < n; ++j) {
        for (i = 0; i < m; ++i) {
            xt[INDEX(j,i,n)] = x[INDEX(i,j,m)];
        }
    }
}
