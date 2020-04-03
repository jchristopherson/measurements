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
 *
 * @param x The data set.
 *
 * @param n The number of data points.
 * @param x An N-element array containing the data set.
 *
 * @par Remarks
 * To avoid overflow-type issues, Welford's algorithm is employed.  A 
 * simple illustration of this algorithm can be found
 * [here](https://www.johndcook.com/blog/standard_deviation/).
 */
double c_variance(int n, const double *x);

/**
 * Computes the standard deviation of a data set.
 *
 * @param n The number of data points.
 * @param x An N-element array containing the data set.
 *
 * @return The standard deviation of @p x.
 */
double c_standard_deviation(int n, const double *x);

/**
 * Computes the range of a data set.
 *
 * @param n The number of data points.
 * @param x An N-element array containing the data set.
 *
 * @return The range of the data set.
 */
double c_data_range(int n, const double *x);

/**
 * Computes the z-score typically used in confidence interval
 * calculations.
 *
 * @param c The confidence level.  This value must be between 0
 *  and 1 such that: 0 < c < 1.
 * @param z The computed z-score.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_INVALID_INPUT_ERROR: Occurs if @p c is not within its allowed
 *      range.
 */
int c_z_score(double c, double *z);

/**
 * Evaluates the probability distribution function of a normal
!! distribution.
!!
!! @param mu The population mean.
!! @param sigma The population standard deviation.
!! @param n The number of values at which to evaluat the function.
!! @param x An N-element array containing the values at which to 
!!  evaluate the distrubition funciton.
!! @param f An N-element array where the function output will be
!!  written.
 */
void c_normal_distribution(double mu, double sigma, int n, const double *x,
    double *f);

/**
 * Evalautes the probability distribution function of Student's
!! t-distribution.
!!
!! @param dof The number of degrees of freedom of the data set.
!! @param n The number of values at which to evaluat the function.
!! @param t An N-element array containing the values at which to 
!!  evaluate the distrubition funciton.
!! @param f An N-element array where the function output will be
!!  written.
 */
void c_t_distribution(int dof, int n, const double *t, double *f);

#ifdef __cplusplus
}
#endif  // __cplusplus
#endif  // MEASUREMENTS_H_
