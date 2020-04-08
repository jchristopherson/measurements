#ifndef MEASUREMENTS_H_
#define MEASUREMENTS_H_

#include <stdbool.h>

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

/** A type containing variance components describing a measurement process. */
typedef struct {
    /** The measurement variation component.  In a gauge analysis this
     * is referred to as the gauge R&R.  It is the sum of the repeatability
     * and reproducibility variance components. */
    double measurement_variance;
    /** The part variance component. */
    double part_variance;
    /** The total process variance.  This is the sum of the measurement 
     * variance and part variance. */
    double total_variance;
    /** The equipment variance component.  This is often referred to
     * as the repeatability component. */
    double equipment_variance;
    /** The operator variance component. */
    double operator_variance;
    /** The operator by part variance component. */
    double operator_by_part_variance;
} process_variance;

/** A type containing gauge R&R results. */
typedef struct {
    /** The precision to tolerance ratio (P/T ratio).  This ratio
     * is simply the ratio of the measurement standard deviation to
     * the tolerance range. */
    double pt_ratio;
    /** The precision to total variation ratio (P/TV ratio).  This 
     * ratio is simply the ratio of the measurement standard deviation to
     * the total process standard deviation. */
    double ptv_ratio;
    /** The tolerance range. */
    double tolerance_range;
} grr_results;

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
 * distribution.
 *
 * @param mu The population mean.
 * @param sigma The population standard deviation.
 * @param n The number of values at which to evaluat the function.
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
 * Utilizes an analysis of variance (ANOVA) to determine the
 * variance components of a measurement process represented by the
 * supplied data set.
 *
 * @param nparts The number of parts being analyzed (must be greater 
 *  than 1).
 * @param ntests The number of tests performed on each part (must be
 *  greater than 1).
 * @param nops The number of operators performing tests (must be at
 *  least 1).
 * @param x An NPARTS-by-NTESTS-by-NOPS data set from the measurement 
 *  process to analyze.
 * @param alpha A parameter used to determine the appropriate 
 *  calculation path.  The recommended value is 0.05 (95% confidence level).
 * @param rst The resulting variance components.
 *
 * @remarks It is possible for this routine to return zero-valued
 *  variance components.  In such an event it is recommended that
 *  another variance estimator is utilized.
 */
void c_anova(int nparts, int ntests, int nops, const double *x, double alpha,
    process_variance *rst);

/**
 * Utilizes a control chart type approach to evaluate the 
 * measurement process utilized to collect the supplied data set.
 *
 * @param nparts The number of parts being analyzed (must be greater 
 *  than 1).
 * @param ntests The number of tests performed on each part (must be
 *  greater than 1).
 * @param nops The number of operators performing tests (must be at
 *  least 1).
 * @param x An NPARTS-by-NTESTS-by-NOPS data set from the measurement 
 *  process to analyze.
 * @param rst The resulting variance components.
 *
 * @remarks This approach is best suited for smaller data sets whose
 *  dimension doesn't exceed 25-30.  Anything over this size is better
 *  served by another technique.
 */
void c_control_chart_variance(int nparts, int ntests, int nops, 
    const double *x, process_variance *rst);

/**
 * Computes the gauge R&R statistics given a supplied set of
 * process variance data.
 *
 * @param k The multiplier to use when computing the P/T ratio.
 *  Typically, a value of 6 is used for this factor.
 * @param x The process variance data.
 * @param usl The upper specification limit.
 * @param lsl The lower specification limit.
 * @param grr The gauge R&R statistics.
 */
void c_compute_grr(double k, const process_variance *x, double usl, double lsl,
    grr_results *grr);

/**
 * Computes the discrimination ratio.
 *
 * @param tv The total variance.
 * @param mv The measurement system variance.
 *
 * @return The results of the operation.
 */
double c_discrimination_ratio(double tv, double mv);

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

#ifdef __cplusplus
}
#endif  // __cplusplus
#endif  // MEASUREMENTS_H_
