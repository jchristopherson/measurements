#ifndef MEASUREMENTS_H_
#define MEASUREMENTS_H_

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Tests to see if an array is monotonically increasing or 
 * decreasing.
 *
 * @param n The number of data points in the array.
 * @param x The N-element array to test.
 * 
 * @return Returns true if the array is monotonically increasing or
 * monotonically decreasing; else, returns false.
 */
bool c_is_monotonic(int n, const double *x);

/**
 * Computes the mean of a data set.
 *
 * @param n The number of data points.
 * @param x An N-element array containing the data set.
 *
 * @return The mean of @p x.
 */
double c_mean(int n, const double *x);

/**
 * Computes the median of a data set.
 *
 * @param x The data set.
 *
 * @return The median of @p x.
 */
double c_median(int n, const double *x);

/**
 * Computes the sample variance of a data set.
!!
!! @param x The data set.
!!
!! @param n The number of data points.
!! @param x An N-element array containing the data set.
!!
!! @par Remarks
!! To avoid overflow-type issues, Welford's algorithm is employed.  A 
!! simple illustration of this algorithm can be found
!! [here](https://www.johndcook.com/blog/standard_deviation/).
 */
double c_variance(int n, const double *x);

/**
 * Computes the standard deviation of a data set.
!!
!! @param n The number of data points.
!! @param x An N-element array containing the data set.
!!
!! @return The standard deviation of @p x.
 */
double c_standard_deviation(int n, const double *x);

#ifdef __cplusplus
}
#endif  // __cplusplus
#endif  // MEASUREMENTS_H_
