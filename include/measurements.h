#ifndef MEASUREMENTS_H_
#define MEASUREMENTS_H_

#include <stdbool.h>
#include <complex.h>

/** A flag denoting normal operation - no error. */
#define M_NO_ERROR                          0
/** A flag denoting an invalid input error state. */
#define M_INVALID_INPUT_ERROR               10000
/** A flag denoting an array size error. */
#define M_ARRAY_SIZE_ERROR                  10001
/** A flag denoting an out-of-memory condition. */
#define M_OUT_OF_MEMORY_ERROR               10002
/** A flag denoting a non-monotonic array error. */
#define M_NONMONOTONIC_ARRAY_ERROR          10003
/** A flag denoting that no data has been defined. */
#define M_NO_DATA_DEFINED_ERROR             10004
/** A flag denoting an underdefined problem error. */
#define M_UNDERDEFINED_PROBLEM              10005

/** 
 * Indicates that the spline is quadratic over the interval under
 * consideration (beginning or ending interval).  This is equivalent to
 * allowing a "natural" boundary condition at either the initial or final
 * point. 
 */
#define SPLINE_QUADRATIC_OVER_INTERVAL      1000
/** 
 * Indicates a known first derivative at either the beginning or ending
 * point. 
 */
#define SPLINE_KNOWN_FIRST_DERIVATIVE       1001
/** 
 * Indicates a known second derivative at either the beginning or ending
 * point. 
 */
#define SPLINE_KNOWN_SECOND_DERIVATIVE      1002
/** 
 * Indicates a continuous third derivative at either the beginning or ending
 * point. 
 */
#define SPLINE_CONTINUOUS_THIRD_DERIVATIVE  1003

/** Defines an equal variance assumption. */
#define EQUAL_VARIANCE_ASSUMPTION           2000
/** Defines an unequal variance assumption. */
#define UNEQUAL_VARIANCE_ASSUMPTION         2001
/** Defines a paired data set assumption. */
#define PAIRED_DATA_SET_ASSUMPTION          2002

/** 
 * A type containing information regarding a given statistic and its
 * associated probability.
 */
typedef struct {
    /** The value of the statistic. */
    double value;
    /** The probability defining the significance of the statistic value. */
    double probability;
} statistic;

/** A single entry (row) in an ANOVA table. */
typedef struct {
    /** The number of degrees of freedom. */
    int dof;
    /** The sum of the squares. */
    double sum_of_squares;
    /** The mean of the squares. */
    double mean_of_squares;
    /** The F statistic. */
    double f_stat;
    /** The mean value. */
    double mean;
    /** The probability (p-value) that the null hypothesis is true. */
    double probability;
} anova_table_entry;

/** A basic ANOVA table. */
typedef struct {
    /** The results from comparison of variances between the data sets. */
    anova_table_entry between;
    /** The residual variation. */
    anova_table_entry residual;
    /** The total variation information. */
    anova_table_entry total;
} anova_table;

/** A Gage R&R table. */
typedef struct {
    /** The combined operator information. */
    anova_table_entry operators;
    /** The combined part information. */
    anova_table_entry parts;
    /** The operator-by-part information. */
    anova_table_entry operator_by_part;
    /** The measurement equipment information. */
    anova_table_entry equipment;
    /** The total variability information. */
    anova_table_entry total;
} gage_table;

/**
 * Defines a window function.
 *
 * @param bin The index or bin number (0 <= @p bin <= @p n).
 * @param n The transform length.
 * @return The window function value.
 */
typedef double (*c_window_function)(int bin, int n);

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
 * Computes the t-score typically used in confidence interval
 * calculations when the population size is limited.
 *
 * @param c The confidence level.  This value must be between 0
 *  and 1 such that: 0 < c < 1.
 * @param n The population size.
 * @param t The computed t-score.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_INVALID_INPUT_ERROR: Occurs if @p c is not within its allowed
 *      range.
 */
int c_t_score(double c, int n, double *t);

/**
 * Computes the confidence interval of a data set.
 *
 * @param n The number of data points.
 * @param x An N-element array containing the data set.
 * @param zval The critical value (z or t value).
 *
 * @return The confidence interval as referenced from the population
 *  mean.
 */
double c_confidence_interval(int n, const double *x, double zval);

/**
 * Evaluates the probability distribution function of the
 * F distribution.
 *
 * @param d1 A model parameter.
 * @param d2 A model parameter.
 * @param n The number of values at which to evaluat the function.
 * @param x An N-element array containing the values at which to 
 *  evaluate the distrubition funciton.
 * @param f An N-element array where the function output will be
 *  written.
 */
void c_f_distribution(double d1, double d2, int n, const double *x, double *f);

/**
 * Applies Student's t-test to compute the t-statistic.  A two-tailed 
 * distribution is assumed.
 *
 * @param n1 The number of data points in the first data set.
 * @param x1 An @p n1 element array containing the first data set.
 * @param n2 The number of data points in the second data set.
 * @param x2 An @p n2 element array containing the second data set.
 * @param method An input flag defining which method to utilize.
 *  - EQUAL_VARIANCE_ASSUMPTION: This flag enforces an assumption that
 *      the variances of both populations are equivalent.
 *  - UNEQUAL_VARIANCE_ASSUMPTION: This flag enforces an assumption 
 *      that the varainces of both populations are not necessarily
 *      equivalent.
 *  - PAIRED_DATA_SET_ASSUMPTION: This flag enforces an assumption that
 *      the data sets are paired.  This requires that both data sets
 *      are the same size.  If this flag is defined, and the supplied
 *      data sets are different sized, the routine switches to
 *      EQUAL_VARIANCE_ASSUMPTION.
 * If no value is specified, the EQUAL_VARIANCE_ASSUMPTION is utilized.
 * @param rst Student's t-statistic and the associated probability term
 * that establishes the significance of the t-statistic.
 */
void c_t_test(int n1, const double *x1, int n2, const double *x2, int method, 
    statistic *rst);

/**
 * Applies the F-test for different variances.
 *
 * @param n1 The number of data points in the first data set.
 * @param x1 An @p n1 element array containing the first data set.
 * @param n2 The number of data points in the second data set.
 * @param x2 An @p n2 element array containing the second data set.
 * @param rst The F-test statistic and associated probability term that
 * establishes the signficance of the f-statistic.
 */
void c_f_test(int n1, const double *x1, int n2, const double *x2, 
    statistic *rst);

/**
 * Removes all NaN's from an array.
 * 
 * @param nx The array length.
 * @param x The array of length @p nx on which to operate.
 * @param y An array of length @p nx where the results will be written.
 * @param ny The actual number of items written to @p y.
 */
void c_remove_nans(int nx, const double *x, double *y, int *ny);

/**
 * Removes all zeros from an array.
 * 
 * @param nx The array length.
 * @param x The array of length @p nx on which to operate.
 * @param y An array of length @p nx where the results will be written.
 * @param ny The actual number of items written to @p y.
 */
void c_remove_zeros(int nx, const double *x, double *y, int *ny);

/**
 * Computes the R-squared value of a data set and a model of
 * the data.
 *
 * @param n The array length.
 * @param y An N-element array containing the dependent variables 
 *  from the data set.
 * @param ym An N-element array containing the corresponding modeled
 *  values.
 *
 * @return The R-squared value.
 */
double c_r_squared(int n, const double *y, const double *ym);

/**
 * Computes the adjusted R-squared value of a data set and a model of
 * the data.
 *
 * @param p The number of model parameters.
 * @param n The array length.
 * @param y An N-element array containing the dependent variables 
 *  from the data set.
 * @param ym An N-element array containing the corresponding modeled
 *  values.
 *
 * @return The R-squared value.
 */
double c_adjusted_r_squared(int p, int n, const double *y, const double *ym);

/**
 * Computes the value of the regularized beta function.
 *
 * @param x The upper limit of the integration.
 * @param a The first argument of the function.
 * @param b The second argument of the function.
 *
 * @return The value of the regularized beta function at @p a and @p b.
 */
double c_regularized_beta(double x, double a, double b);

/**
 * Computes the beta function.
 *
 * @param a The first argument of the function.
 * @param b The second argument of the function.
 *
 * @return The value of the beta function at @p a and @p b.
 */
double c_beta(double a, double b);

/**
 * Computes the incomplete beta function.
 *
 * @param x The upper limit of the integration.
 * @param a The first argument of the function.
 * @param b The second argument of the function.
 *
 * @return The value of the incomplete beta function.
 */
double c_incomplete_beta(double x, double a, double b);

/**
 * Computes the digamma function \f$ \psi(x) = 
 * \frac{d}{dx}\left( \ln \left( \Gamma \left( x \right) \right) 
 * \right) \f$.
 *
 * @param x The value at which to evaluate the function.
 * @return The function value.
 */
double c_digamma(double x);

/**
 * Performs a polynomial interpolation.
 *
 * @param order The order of the polynomial.  This value must be at 
 *  least 1, but not exceed @p npts - 1.
 * @param npts The number of data points.
 * @param x An @p npts element array containing the independent variable
 *  data.  This array must be monotonic.
 * @param y An @p npts element array containing the dependent variable
 *  data.
 * @param ni The number of points to interpolate.
 * @param xi An @p ni element array containing the independent variable
 *  values at which to interpolate.
 * @param yi An @p ni element array corresponding to @p xi where the
 *  interpolated values will be written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - M_INVALID_INPUT_ERROR: Occurs if @p order is less than 1.
 *  - M_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
 *      increasing or decreasing.
 */
int c_interpolate(int order, int npts, const double *x, const double *y,
    int ni, const double *xi, double *yi);

/**
 * Performs a spline interpolation.
 *
 * @param npts The number of data points.
 * @param x An @p npts element array containing the independent variable
 *  data.  This array must be monotonic.
 * @param y An @p npts element array containing the dependent variable
 *  data.
 * @param ni The number of points to interpolate.
 * @param xi An @p ni element array containing the independent variable
 *  values at which to interpolate.
 * @param yi An @p ni element array corresponding to @p xi where the
 *  interpolated values will be written.
 * @param ibcbeg An input that defines the nature of the
 *  boundary condition at the beginning of the spline.
 *  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
 *      initial interval.  No value is required for @p ybcbeg.  This is
 *      often considered a natural boundary condition.
 *  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at 
 *      its initial point is provided in @p ybcbeg.
 *  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at 
 *      its initial point is provided in @p ybcbeg.
 *  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is 
 *      continuous at x(2).  No value is required for @p ybcbeg.
 * @param ybcbeg If needed, the value of the initial point boundary
 *  condition.
 * @param ibcend An input that defines the nature of the
 *  boundary condition at the end of the spline.
 *  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
 *      final interval.  No value is required for @p ybcend.  This is
 *      often considered a natural boundary condition.
 *  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at 
 *      its initial point is provided in @p ybcend.
 *  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at 
 *      its initial point is provided in @p ybcend.
 *  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is 
 *      continuous at x(n-1).  No value is required for @p ybcend.
 * @param ybcend If needed, the value of the final point boundary
 *  condition.  If needed, but not supplied, a default value of zero 
 *  will be used.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - M_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
 *      increasing or decreasing.
 */
int c_spline(int npts, const double *x, const double *y, int ni, 
    const double *xi, double *yi, int ibcbeg, double ybcbeg, int ibcend,
    double ybcend);

/**
 * Smooths a data set using a robust locally weighted scatterplot 
 * smoothing (LOWESS) algorithm. 
 *
 * @param npts The number of data points.
 * @param x An @p npts element array containing the independent variable
 *  data.  This array must be monotonic.
 * @param y An @p npts element array containing the dependent variable
 *  data.
 * @param factor Specifies the amount of smoothing.  More specifically, 
 *  this value is the fraction of points used to compute each value.  
 *  As this value increases, the output becomes smoother.  Choosing a 
 *  value in the range of 0.2 to 0.8 usually results in a good fit.  As 
 *  such, a reasonable starting point, in the absence of better 
 *  information, is a value of 0.5.
 * @param ys An @p npts element array where the smoothed data will
 *  be written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - M_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
 *      increasing or decreasing.
 */
int c_lowess_smoothing(int npts, const double *x, const double *y, 
    double factor, double *ys);

/**
 * Applies a moving average to smooth a data set.
 *
 * @param npts The number of data points.
 * @param x An @p npts element array that on input contains the 
 *  signal to smooth.  On output, the smoothed signal.
 * @param navg The size of the averaging window.  This value must be
 *  at least 2, but no more than the number of elements in @p x.
 *
 * @return  An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_INVALID_INPUT_ERROR: Occurs if @p navg is less than 2, or 
 *      greater than @p npts.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 */
int c_moving_average(int npts, double *x, int navg);

/**
 * Fits the multiple input, multiple output linear model
 * A * X = Y by solving for matrix A in a least-squares sense.
 *
 * @param m The number of rows in matrix A.
 * @param n The number of columns in matrix A.
 * @param k The number of data points to fit (number of columns in
 *  either X or Y).
 * @param x An N-by-K matrix of known independent variables.  K
 *  must be greater than or equal to N.
 * @param ldx The leading dimension of matrix X.
 * @param y An M-by-K matrix of known dependent variables.  Notice,
 *  M must be less than or equal to N, and K must be greater than or
 *  equal to M.
 * @param ldy The leading dimension of matrix Y.
 * @param The M-by-N coefficient matrix A.
 * @param lda The leading dimension of matrix A.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - M_ARRAY_SIZE_ERROR: Occurs if there is a size-mismatch in the
 *      matrix equation.
 *  - M_UNDERDEFINED_PROBLEM: Occurs if there is insufficient data 
 *      (e.g. k < n).
 */
int c_linear_least_squares_mimo(int m, int n, int k, const double *x, int ldx,
    const double *y, int ldy, double *a, int lda);

/**
 * Fits the multiple input, single output linear model
 * A * X = Y by solving for array A in a least-squares sense.
 *
 * @param n The number of coefficients to find.
 * @param k The number of data points to fit.
 * @param x An N-by-K matrix of known independent variables.  K
 *  must be greater than or equal to N.
 * @param ldx The leading dimension of matrix X.
 * @param y A K-element array containing the known dependent 
 *  variables.
 * @param a The N-element coefficient array A.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - M_ARRAY_SIZE_ERROR: Occurs if there is a size-mismatch in the
 *      matrix equation.
 *  - M_UNDERDEFINED_PROBLEM: Occurs if there is insufficient data 
 *      (e.g. k < n).
 */
int c_linear_least_squares_miso(int n, int k, const double *x, int ldx,
    const double *y, double *a);

/**
 * Attempts to locate all peaks and valleys within the given
 * data set bound by the constraints specified.
 *
 * @param n The number of data points in the array to search.
 * @param x The N-element array to search.
 * @param delta A threshold level used to denote the minimum change
 *  acceptable in determining a peak or valley from neighboring points.
 * @param szmxind The size of the @p mxind buffer.
 * @param mxind A @p szmxind element array where the indices of the
 *  peak values will be written.  The indices are zero-based.
 * @param nmxind The actual number of peak value indices written to
 *  @p mxind.
 * @param szmnind The size of the @p mnind buffer.
 * @param mnind A @p szmnind element array where the indices of the
 *  valley values will be written.  The indices are zero-based.
 * @param nmnind The actual number of valley value indices written to
 *  @p mnind.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 */
int c_peak_detect(int n, const double *x, double delta, int szmxind,
    int *mxind, int *nmxind, int szmnind, int *mnind, int *nmnind);

/**
 * Tests to see if a number is an integer power of two.
 *
 * @param n The integer to test.
 * @return Returns true if @p n is a power of two; else, false.
 */
bool c_is_power_of_two(int n);

/**
 * Provides the next higher integer power of two.
 *
 * @param x The value to test.
 * @return The next power of two higher than @p x.  If @p x is already
 * a power of two, it's value is simply returned.  For instance, if @p
 * is set to 128, then a value of 7 is returned.  However, if a value
 * of 129 is supplied, then a value of 8 is returned.
 */
int c_next_power_of_two(int x);

/**
 * Computes the Fourier transform of a discretely sampled signal.
 *
 * @param n The number of data points.
 * @param x The N-element array containing the data to transform.
 * @param nf The number of elements in @p f.  Ideally, this should be
 *  N / 2 + 1 if N is even; else, (N + 1) / 2 if N is odd.
 * @param f The positive half and DC component of the transformed data.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 */
int c_fourier_transform(int n, const double *x, int nf, double complex *f);

/**
 * Computes the periodogram power spectral density (PSD) estimate
 * of a signal.
 *
 * @param n The number of data points.
 * @param x The N-element array containing the data to transform.
 * @param winfun The window function to apply.
 * @param nfft The length Fourier transform to apply.  This must
 *  be an integer power of two, even if @p x is not an integer power
 *  of two in length.  If this parameter is less than the length of
 *  @p x, @p x will be overlapped, windowed, and averaged to achieve
 *  an estimate of the power spectrum.  If this parameter is larger
 *  than the length of @p x, @p x will be padded with zeros prior
 *  to windowing.
 * @param np
 * @param p The NP-element array containing the periodogram.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - M_INVALID_INPUT_ERROR: Occurs if @p nfft is not an integer power
 *      of two.
 */
int c_periodogram(int n, const double *x, c_window_function winfun, int nfft,
    int np, double *p);

/**
 * Computes a suitable frequency for a Fourier transformed data set.
 *
 * @param fs The rate at which the data was sampled.  Notice, the
 *  returned frequency value will be expressed in the same units as this
 *  value.
 * @param i The frequency bin such that 0 <= i < m where m is
 *  @p nxfrm / 2 + 1 if @p nxfrm is even; else, (@p nxfrm + 1) / 2 if
 *  @p nxfrm is odd.
 * @param nxfrm The length of the signal that was transformed.
 *
 * @return The frequency value.
 */
double c_fourier_frequency(double fs, int i, int nxfrm);

/**
 * Defines a rectangular window.
 *
 * @param j The index or bin number (0 <= @p bin <= @p n).
 * @param n The transform length.
 *
 * @return The value of the window function at index @p j.
 */
double c_rectangular_window(int j, int n);

/**
 * Defines a Hann window.
 *
 * @param j The index or bin number (0 <= @p bin <= @p n).
 * @param n The transform length.
 *
 * @return The value of the window function at index @p j.
 */
double c_hann_window(int j, int n);

/**
 * Defines a Hamming window.
 *
 * @param j The index or bin number (0 <= @p bin <= @p n).
 * @param n The transform length.
 *
 * @return The value of the window function at index @p j.
 */
double c_hamming_window(int j, int n);

/**
 * Defines a Welch window.
 *
 * @param j The index or bin number (0 <= @p bin <= @p n).
 * @param n The transform length.
 *
 * @return The value of the window function at index @p j.
 */
double c_welch_window(int j, int n);

/**
 * Defines a Blackman-Harris window.
 *
 * @param j The index or bin number (0 <= @p bin <= @p n).
 * @param n The transform length.
 *
 * @return The value of the window function at index @p j.
 */
double c_blackman_harris_window(int j, int n);

/**
 * Applies a low-pass filter to a signal.
 *
 * @param n The length of the input array.
 * @param x An N-element array containing the signal to filter.
 * @param fs The frequency at which @p x was sampled.
 * @param cutoff The cut-off frequency.  This value must be
 *  positive-valued, and must be less than the Nyquist frequency.
 * @param y An N-element array where the filtered signal will be
 *  written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff is either not positive
 *      valued, or is greater than or equal to the Nyquist frequency.
 */
int c_low_pass_filter(int n, const double *x, double fs, double cutoff,
    double *y);

/**
 * Applies a high-pass filter to a signal.
 *
 * @param n The length of the input array.
 * @param x An N-element array containing the signal to filter.
 * @param fs The frequency at which @p x was sampled.
 * @param cutoff The cut-off frequency.  This value must be
 *  positive-valued, and must be less than the Nyquist frequency.
 * @param y An N-element array where the filtered signal will be
 *  written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff is either not positive
 *      valued, or is greater than or equal to the Nyquist frequency.
 */
int c_high_pass_filter(int n, const double *x, double fs, double cutoff,
    double *y);

/**
 * Applies a band-pass filter to a signal.
 *
 * @param n The length of the input array.
 * @param x An N-element array containing the signal to filter.
 * @param fs The frequency at which @p x was sampled.
 * @param cutoff1 The lower cut-off frequency.  This value must be
 *  positive-valued, and must be less than the Nyquist frequency.
 * @param cutoff2 The upper cut-off frequency.  This value must be
 *  positive-valued, and must be less than the Nyquist frequency.
 * @param y An N-element array where the filtered signal will be
 *  written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff1 or @p cutoff2 is 
 *      either not positive-valued, or is greater than or equal to the 
 *      Nyquist frequency.
 */
int c_band_pass_filter(int n, const double *x, double fs, double cutoff1,
    double cutoff2, double *y);

/**
 * Applies a band-stop filter to a signal.
 *
 * @param n The length of the input array.
 * @param x An N-element array containing the signal to filter.
 * @param fs The frequency at which @p x was sampled.
 * @param cutoff1 The lower cut-off frequency.  This value must be
 *  positive-valued, and must be less than the Nyquist frequency.
 * @param cutoff2 The upper cut-off frequency.  This value must be
 *  positive-valued, and must be less than the Nyquist frequency.
 * @param y An N-element array where the filtered signal will be
 *  written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff1 or @p cutoff2 is 
 *      either not positive-valued, or is greater than or equal to the 
 *      Nyquist frequency.
 */
int c_band_stop_filter(int n, const double *x, double fs, double cutoff1,
    double cutoff2, double *y);

/**
 * Computes an FFT of a data set.  The results of the transform
 * are normalized such that an inverse transform will result in the 
 * original signal.
 *
 * @param n The length of the array.
 * @param x An N-element array containing the data to transform.
 * @param f An N-element array where the transformed data will be 
 *  written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 */
int c_fft(int n, const double *x, double complex *f);

/**
 * Computes the inverse FFT of a data set.
 *
 * @param n The length of the array.
 * @param x An N-element array containing the data to transform.
 * @param f An N-element array where the transformed data will be 
 *  written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 */
int c_ifft(int n, const double complex *x, double complex *f);

/**
 * Computes the cumulative sum of an array.
 *
 * @param n The length of the array.
 * @param x The N-element array on which to operate.
 * @param y The N-element array where the solution will be written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 */
int c_cumulative_sum(int n, const double *x, double *y);

/**
 * Computes the differences between each element in the array.
 *
 * @param n The length of the array.
 * @param x The N-element array on which to operate.
 * @param dx An N-1 element array where the solution will be written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 */
int c_difference(int n, const double *x, double *dx);

/**
 * Unwraps an array of phase angles.
 *
 * @param n The array length.
 * @param x The N-element array to unwrap.
 * @param cutoff An input that specifies the threshold value to use 
 *  when unwrapping steps in @p x.
 * @param y The N-element array where the results will be written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 */
int c_unwrap(int n, const double *x, double cutoff, double *y);

/**
 * Estimates the derivative of a discretely sampled signal by means
 * of finite differences.
 *
 * @param n The length of the arrays.
 * @param x The N-element array containing the independent variable 
 *  data.
 * @parma y The N-element array containing the dependent variable
 *  data.
 * @param dydx An N-element array where the results will be written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 */
int c_finite_difference(int n, const double *x, const double *y, double *dydx);

/**
 * Estimates the indefinite integral of a discretely sampled signal.
 *
 * @param n The length of the arrays.
 * @param x The N-element array containing the independent variable 
 *  data.
 * @parma y The N-element array containing the dependent variable
 *  data.
 * @param c The initial condition.
 * @param f An N-element array where the results will be written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 */
int c_integrate(int n, const double *x, const double *y, double c, double *f);

/**
 * Estimates the definite integral of a discretely sampled signal
 * using a trapezoidal approach.
 *
 * @param n The length of the arrays.
 * @param x The N-element array containing the independent variable 
 *  data.
 * @parma y The N-element array containing the dependent variable
 *  data.
 *
 * @return The result of the integration.
 */
double c_trapz_integrate(int n, const double *x, const double *y);

/**
 * Removes the DC offset from a signal.
 *
 * @param n The array length.
 * @param x The N-element array on which to operate.
 * @param y An N-element array where the results will be written.
 *
 * @return An error flag with the following possible values.
 *  - M_NO_ERROR: No error occurred.  Normal operation.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 */
int c_remove_dc_offset(int n, const double *x, double *y);

/**
 * Computes a one-way analysis of variance (ANOVA).
 *
 * @param npts The number of data points in each data set.
 * @param nsets The number of data sets.
 * @param x An NPTS-by-NSETS matrix of data to analyze.
 * @param ldx The leading dimension of @p x.  This value must be
 *  greater than or equal to @p npts.
 * @param tbl The ANOVA table to populate with results.
 */
void c_anova(int npts, int nsets, const double *x, int ldx, anova_table *tbl);

/**
 * Computes a crossed-effects analysis of variance (ANOVA) for a set of 
 * measurement data in order to better understand variance contributions of
 * each part of a measurement system or gage analysis.
 * 
 * @param nparts The number of parts tested (must be greater than 1).
 * @param ntests The number of times each part is tested (must be greater 
 *  than 1).
 * @param nops The number of operators performing the tests (must be greater 
 *  than 1).
 * @param x An NPARTS-by-NTESTS-by-NOPS column-major 3D array containing the 
 *  data.
 * @param rst The gage_table to populate with results.
 * 
 * @return An error flag with the following possible values.
 *  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - M_INSUFFICIENT_DATA_ERROR: Occurs if there is only 1 row or 1
 *      column in @p x.
 */
int c_gage_anova(int nparts, int ntests, int nops, const double *x, 
    gage_table *rst);

/**
 * Evaluates the probability distribution function of a normal
 * distribution.
 *
 * @param mu The population mean.
 * @param sigma The population standard deviation.
 * @param n The number of values at which to evaluate the function.
 * @param x An N-element array containing the values at which to 
 *  evaluate the distrubition funciton.
 * @param f An N-element array where the function output will be
 *  written.
 */
void c_normal_distribution(double mu, double sigma, int n, const double *x,
    double *f);

/**
 * Evalautes the probability distribution function of Student's
 * t-distribution.
 *
 * @param dof The number of degrees of freedom of the data set.
 * @param n The number of values at which to evaluat the function.
 * @param t An N-element array containing the values at which to 
 *  evaluate the distrubition funciton.
 * @param f An N-element array where the function output will be
 *  written.
 */
void c_t_distribution(double dof, int n, const double *t, double *f);

/**
 * Evaluates the probability distribution function of a beta
 * distribution.
 *
 * @param a The first argument of the function.
 * @param b The second argument of the function.
 * @param n The number of values at which to evaluat the function.
 * @param x An N-element array containing the values at which to 
 *  evaluate the distrubition funciton.
 * @param f An N-element array where the function output will be
 *  written.
 */
void c_beta_distribution(double a, double b, int n, const double *x, double *f);

/** Evaluates the probability distribution function of the F-distribution.
 * 
 * @param d1 A model parameter.
 * @param d2 A model parameter.
 * @param n The number of values at which to evaluat the function.
 * @param x An N-element array containing the values at which to 
 *  evaluate the distrubition funciton.
 * @param f An N-element array where the function output will be
 *  written.
 */
void c_f_distribution(double d1, double d2, int n, const double *x, double *f);


#ifdef __cplusplus
}
#endif  // __cplusplus
#endif  // MEASUREMENTS_H_
