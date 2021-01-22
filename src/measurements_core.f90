! measurements_core.f90

!> @brief A module containing types and routines to support measurement-related 
!! calculations.
module measurements_core
    use iso_fortran_env
    use ferror
    implicit none

    !> @brief A flag denoting normal operation - no error.
    integer(int32), parameter :: M_NO_ERROR = 0
    !> @brief A flag denoting an invalid input error state.
    integer(int32), parameter :: M_INVALID_INPUT_ERROR = 10000
    !> @brief A flag denoting an array size error.
    integer(int32), parameter :: M_ARRAY_SIZE_ERROR = 10001
    !> @brief A flag denoting an out-of-memory condition.
    integer(int32), parameter :: M_OUT_OF_MEMORY_ERROR = 10002
    !> @brief A flag denoting a non-monotonic array error.
    integer(int32), parameter :: M_NONMONOTONIC_ARRAY_ERROR = 10003
    !> @brief A flag denoting that no data has been defined.
    integer(int32), parameter :: M_NO_DATA_DEFINED_ERROR = 10004
    !> @brief A flag denoting an underdefined problem error.
    integer(int32), parameter :: M_UNDERDEFINED_PROBLEM = 10005
    !> @brief A flag denoting an insufficient data error.
    integer(int32), parameter :: M_INSUFFICIENT_DATA_ERROR = 10006

    !> Indicates that the spline is quadratic over the interval under
    !! consideration (beginning or ending interval).  This is equivalent to
    !! allowing a "natural" boundary condition at either the initial or final
    !! point.
    integer(int32), parameter :: SPLINE_QUADRATIC_OVER_INTERVAL = 1000
    !> Indicates a known first derivative at either the beginning or ending
    !! point.
    integer(int32), parameter :: SPLINE_KNOWN_FIRST_DERIVATIVE = 1001
    !> Indicates a known second derivative at either the beginning or ending
    !! point.
    integer(int32), parameter :: SPLINE_KNOWN_SECOND_DERIVATIVE = 1002
    !> Indicates a continuous third derivative at either the beginning or ending
    !! point.
    integer(int32), parameter :: SPLINE_CONTINUOUS_THIRD_DERIVATIVE = 1003

    !> @brief Defines an equal variance assumption.
    integer(int32), parameter :: EQUAL_VARIANCE_ASSUMPTION = 2000
    !> @brief Defines an unequal variance assumption.
    integer(int32), parameter :: UNEQUAL_VARIANCE_ASSUMPTION = 2001
    !> @brief Defines a paired data set assumption.
    integer(int32), parameter :: PAIRED_DATA_SET_ASSUMPTION = 2002

    !> @brief A type containing information regarding a given statistic and its
    !! associated probability.
    type, bind(C) :: statistic
        !> @brief The value of the statistic.
        real(real64) :: value
        !> @brief The probability defining the significance of the statistic
        !! value.
        real(real64) :: probability
    end type

! ******************************************************************************
! MEASUREMENT_STATS.F90
! ------------------------------------------------------------------------------
    interface
        !> @brief Tests to see if an array is monotonically increasing or 
        !! decreasing.
        !!
        !! @param[in] x The array to test.
        !! 
        !! @return Returns true if the array is monotonically increasing or
        !! monotonically decreasing; else, returns false.
        pure module function is_monotonic(x) result(rst)
            real(real64), intent(in), dimension(:) :: x
            logical :: rst
        end function

        !> @brief Computes the mean of a data set.
        !!
        !! @param[in] x The data set.
        !!
        !! @return The mean of @p x.
        pure module function mean(x) result(z)
            real(real64), intent(in), dimension(:) :: x
            real(real64) :: z
        end function

        !> @brief Computes the median of a data set.
        !!
        !! @param[in] x The data set.
        !!
        !! @return The median of @p x.
        module function median(x) result(z)
            real(real64), intent(in), dimension(:) :: x
            real(real64) :: z
        end function

        !> @brief Computes the sample variance of a data set.
        !!
        !! @param[in] x The data set.
        !!
        !! @return The variance of @p x.
        !!
        !! @par Remarks
        !! To avoid overflow-type issues, Welford's algorithm is employed.  A 
        !! simple illustration of this algorithm can be found
        !! [here](https://www.johndcook.com/blog/standard_deviation/).
        pure module function variance(x) result(v)
            real(real64), intent(in), dimension(:) :: x
            real(real64) :: v
        end function

        !> @brief Computes the standard deviation of a data set.
        !!
        !! @param[in] x The data set.
        !!
        !! @return The standard deviation of @p x.
        !!
        !! @remarks
        !! The standard deviation is computed as \f$ \sigma = 
        !! \sqrt{\frac{\sum_{i=1}^{n}(x_{i} - \mu)^{2}}{n - 1}} \f$
        pure module function standard_deviation(x) result(s)
            real(real64), intent(in), dimension(:) :: x
            real(real64) :: s
        end function

        !> @brief Computes the range of a data set.
        !!
        !! @param[in] x The data set.
        !!
        !! @return The range of the data set.
        pure module function data_range(x) result(r)
            real(real64), intent(in), dimension(:) :: x
            real(real64) :: r
        end function

        !> @brief Computes the z-score typically used in confidence interval
        !! calculations.
        !!
        !! @param[in] c The confidence level.  This value must be between 0
        !!  and 1 such that: 0 < c < 1.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_INVALID_INPUT_ERROR: Occurs if @p c is not within its allowed
        !!      range.
        !!
        !! @return The z-score corresponding to @p c.
        module function z_score(c, err) result(z)
            real(real64), intent(in) :: c
            class(errors), intent(inout), optional, target :: err
            real(real64) :: z
        end function

        !> @brief Computes the t-score typically used in confidence interval
        !! calculations when the population size is limited.
        !!
        !! @param[in] c The confidence level.  This value must be between 0
        !!  and 1 such that: 0 < c < 1.
        !! @param[in] n The population size.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_INVALID_INPUT_ERROR: Occurs if @p c is not within its allowed
        !!      range.
        !!
        !! @return The t-score corresponding to @p c.
        module function t_score(c, n, err) result(t)
            real(real64), intent(in) :: c
            integer(int32), intent(in) :: n
            class(errors), intent(inout), optional, target :: err
            real(real64) :: t
        end function

        !> @brief Computes the confidence interval of a data set.
        !!
        !! @param[in] x The data set.
        !! @param[in] zval The critical value (z or t value).
        !!
        !! @return The confidence interval as referenced from the population
        !!  mean.
        !!
        !! @remarks
        !! This value is computed as follows.
        !! \f$ CI = z^{*} \frac{\sigma}{\sqrt{n}} \f$
        pure module function confidence_interval(x, zval) result(ci)
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(in) :: zval
            real(real64) :: ci
        end function

        !> @brief Evaluates the probability distribution function of the 
        !! normal distribution.
        !!
        !! @param[in] mu The population mean.
        !! @param[in] sigma The population standard deviation.
        !! @param[in] x The value at which to evaluate the distrubition 
        !!  funciton.
        !! @param[in] comp An optional input, that if set to true, allows
        !!  evaluation of the cumulative distribution function; else, if set
        !!  to false, the probability density function is evaluated.  The 
        !!  default is false such that the probabidlity density function is
        !!  evaluated.
        !!
        !! @return The value of the distribution function at @p x.
        !!
        !! @remarks
        !! The normal distribution has the form 
        !! \f$ f(x) = \frac{1}{\sigma \sqrt{2 \pi}} 
        !! \exp(-\frac{1}{2}(\frac{x-\mu}{\sigma})^{2}) \f$.
        pure elemental module function normal_distribution(mu, sigma, x) &
                result(f)
            real(real64), intent(in) :: mu, sigma, x
            real(real64) :: f
        end function

        !> @brief Evalautes the probability distribution function of 
        !! Student's t-distribution.
        !!
        !! @param[in] dof The number of degrees of freedom of the data set.
        !! @param[in] t The value at which to evaluate the distribution.
        !! @param[in] comp An optional input, that if set to true, allows
        !!  evaluation of the cumulative distribution function; else, if set
        !!  to false, the probability density function is evaluated.  The 
        !!  default is false such that the probabidlity density function is
        !!  evaluated.
        !!
        !! @return The value of the distribution function at @p t.
        !!
        !! @remarks
        !! Student's t-distribution has the form \f$ f(t) = 
        !! \frac{\Gamma(\frac{\nu + 1}{2})}{\sqrt{\nu \pi} 
        !! \Gamma(\frac{\nu}{2})} (1 + \frac{t^{2}}{\nu})^{-\frac{\nu + 1}{2}} 
        !! \f$.
        pure elemental module function t_distribution(dof, t) result(f)
            real(real64), intent(in) :: dof
            real(real64), intent(in) :: t
            real(real64) :: f
        end function

        !> @brief Evaluates the probability distribution function of the 
        !! beta distribution.
        !!
        !! @param[in] a The first argument of the function.
        !! @param[in] b The second argument of the function.
        !! @param[in] x The value at which to evaluate the distrubition 
        !!  funciton.
        !! @param[in] comp An optional input, that if set to true, allows
        !!  evaluation of the cumulative distribution function; else, if set
        !!  to false, the probability density function is evaluated.  The 
        !!  default is false such that the probabidlity density function is
        !!  evaluated.
        !!
        !! @return The value of the distribution function at @p x.
        !!
        !! @remarks The beta distribution has the form \f$ f(x) = 
        !! \frac{x^{a-1} (1 - x)^{b-1}}{\beta(a,b)} \f$.
        pure elemental module function beta_distribution(a, b, x) &
                result(z)
            real(real64), intent(in) :: a, b, x
            real(real64) :: z
        end function

        !> @brief Evaluates the probability distribution function of the 
        !! F-distribution.
        !!
        !! @param[in] d1 A model parameter.
        !! @param[in] d2 A model parameter.
        !! @param[in] x The value at which to evaluate the distrubition 
        !!  funciton.
        !!
        !! @return The value of the distribution function at @p x.
        !!
        !! @remarks The F distribution has the form 
        !! @par
        !! \f$ f(x) = 
        !! \frac{\sqrt{\alpha}}{x \beta(\frac{d_1}{2}, \frac{d_2}{2})} \f$
        !! @par
        !! where
        !! @par
        !! \f$ \alpha = 
        !! \frac{(d_1 x)^{d_1} d_{2}^{d_2}}{(d_1 x + d_2)^{d_1 + d_2}} \f$.
        pure elemental module function f_distribution(d1, d2, x) result(z)
            real(real64), intent(in) :: d1, d2, x
            real(real64) :: z
        end function

        !> @brief Applies Student's t-test to compute the t-statistic.  A 
        !! two-tailed distribution is assumed.
        !!
        !! @param[in] x1 The first data set.
        !! @param[in] x2 The second data set.
        !! @param[in] method An optional input defining which method to utilize.
        !!  - EQUAL_VARIANCE_ASSUMPTION: This flag enforces an assumption that
        !!      the variances of both populations are equivalent.
        !!  - UNEQUAL_VARIANCE_ASSUMPTION: This flag enforces an assumption 
        !!      that the varainces of both populations are not necessarily
        !!      equivalent.
        !!  - PAIRED_DATA_SET_ASSUMPTION: This flag enforces an assumption that
        !!      the data sets are paired.  This requires that both data sets
        !!      are the same size.  If this flag is defined, and the supplied
        !!      data sets are different sized, the routine switches to
        !!      EQUAL_VARIANCE_ASSUMPTION.
        !! If no value is specified, the EQUAL_VARIANCE_ASSUMPTION is utilized.
        !!
        !! @return Student's t-statistic and the associated probability term
        !! that establishes the significance of the t-statistic.
        !!
        !! @remarks
        !! Student's t-test can be used to understand the differences in means
        !! of two populations that have equivalent variances.
        pure module function t_test(x1, x2, method) result(rst)
            real(real64), intent(in), dimension(:) :: x1, x2
            integer(int32), intent(in), optional :: method
            type(statistic) :: rst
        end function

        !> @brief Applies the F-test for different variances.
        !!
        !! @param[in] x1 The first data set.
        !! @param[in] x2 The second data set.
        !!
        !! @return The F-test statistic and associated probability term that
        !! establishes the signficance of the f-statistic.
        pure module function f_test(x1, x2) result(rst)
            real(real64), intent(in), dimension(:) :: x1, x2
            type(statistic) :: rst
        end function

        !> @brief Removes all NaN's from an array.
        !!
        !! @param[in] x The array.
        !! 
        !! @return The array @p x without any NaN's.
        pure module function remove_nans(x) result(rst)
            real(real64), intent(in), dimension(:) :: x
            real(real64), allocatable, dimension(:) :: rst
        end function

        !> @brief Removes all values sufficiently close to zero from an array.
        !!
        !! @param[in] x The array.
        !! @param[in] tol An optional input that specifies the tolerance
        !!  defining acceptable closeness to zero to consider zero.  The 
        !!  default value is twice machine precision.
        !!
        !! @return The array @p x without zero values.
        pure module function remove_zeros(x, tol) result(rst)
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(in), optional :: tol
            real(real64), allocatable, dimension(:) :: rst
        end function

        !> @brief Computes the R-squared value of a data set and a model of
        !! the data.
        !!
        !! @param[in] y An N-element array containing the dependent variables 
        !!  from the data set.
        !! @param[in] ym An N-element array containing the corresponding modeled
        !!  values.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_ARRAY_SIZE_ERROR: Occurs if @p y and @p ym are not the same
        !!      size.
        !!
        !! @return The R-squared value.
        module function r_squared(y, ym, err) result(rst)
            ! Arguments
            real(real64), intent(in), dimension(:) :: y, ym
            class(errors), intent(inout), optional, target :: err
            real(real64) :: rst
        end function

        !> @brief Computes the probability of the null hypothesis being valid
        !! given the results of an f-test.
        !!
        !! @param[in] f The F-statistic.
        !! @param[in] dof1 The number of degrees of freedom in the first data
        !!  set.
        !! @param[in] dof2 The number of degrees of freedom in the second data
        !!  set.
        !!
        !! @return The probability value.  Subtract this value from 1 to 
        !! determine the validity of the null hypothesis.  If the result is
        !! less than the desired alpha, the null hypothesis is invalid.
        module function ftest_probability(f, dof1, dof2) result(rst)
            real(real64), intent(in) :: f
            integer(int32), intent(in) :: dof1, dof2
        end function
    end interface

! ******************************************************************************
! MEASUREMENTS_SPECIAL.F90
! ------------------------------------------------------------------------------
    interface
        !> @brief Computes the beta function.
        !!
        !! @param[in] a The first argument of the function.
        !! @param[in] b The second argument of the function.
        !!
        !! @return The value of the beta function at @p a and @p b.
        !!
        !! @remarks The beta function is related to the gamma function
        !! by the following relationship \f$ \beta(a,b) = 
        !! \frac{\Gamma(a) \Gamma(b)}{\Gamma(a + b)} \f$.
        pure elemental module function beta(a, b) result(z)
            real(real64), intent(in) :: a, b
            real(real64) :: z
        end function

        !> @brief Computes the value of the regularized incomplete beta 
        !! function.
        !!
        !! @param[in] x The upper limit of the integration for the incomplete
        !!  beta function.
        !! @param[in] a The first argument of the function.
        !! @param[in] b The second argument of the function.
        !!
        !! @return The value of the regularized beta function.
        !!
        !! @remarks The regularized beta function is defined as \f$ 
        !! I_{x}(a, b) = \frac{\beta(x; a, b)}{\beta(a, b)} \f$.
        pure elemental module function regularized_beta(x, a, b) &
                result(z)
            real(real64), intent(in) :: x, a, b
            real(real64) :: z
        end function

        !> @brief Computes the incomplete beta function.
        !!
        !! @param[in] x The upper limit of the integration.
        !! @param[in] a The first argument of the function.
        !! @param[in] b The second argument of the function.
        !!
        !! @return The value of the incomplete beta function.
        !!
        !! @remarks The incomplete beta function is defind as \f$ \beta(x;a,b) =
        !! \int_{0}^{x} t^{a-1} (1 - t)^{b-1} dt \f$.
        pure elemental module function incomplete_beta(x, a, b) result(z)
            real(real64), intent(in) :: x, a, b
            real(real64) :: z
        end function
    end interface

! ******************************************************************************
! MEASUREMENTS_INTERP.F90
! ------------------------------------------------------------------------------
    !> @brief Describes an abstract base class allowing for interpolation of X-Y
    !! type data sets.
    !!
    !! @par Notes
    !! This interpolation object is conceptually based upon the interpolation
    !! scheme utilized by the Numerical Recipes in C++ text.
    type, abstract :: interp_manager
    private
        integer(int32) :: m_order
        integer(int32) :: m_savedIndex
        integer(int32) :: m_indexCheck
        logical :: m_correlated
        real(real64), allocatable, dimension(:) :: m_x
        real(real64), allocatable, dimension(:) :: m_y
    contains
        !> @brief Initializes the interp_manager instance.
        procedure, public :: initialize => im_init
        !> @brief Attempts to locate the index in the array providing a lower
        !! bounds to the specified interpolation point.
        procedure, non_overridable, public :: locate => im_locate
        !> @brief Attempts to locate the index in the array providing a lower
        !! bounds to the specified interpolation point.
        procedure, non_overridable, public :: hunt => im_hunt
        !> @brief Interpolates to obtain the function value at the specified
        !!  independent variable.
        generic, public :: interpolate => im_perform, im_perform_array
        !> @brief Performs the actual interpolation.
        procedure(interp_xy), deferred :: raw_interp
        !> @brief Gets the number of stored data points.
        procedure, public :: get_count => im_get_num_pts
        !> @brief Gets the x component of the requested data point.
        procedure, public :: get_x => im_get_x
        !> @brief Gets the y component of the requested data point.
        procedure, public :: get_y => im_get_y

        procedure, non_overridable :: im_perform
        procedure, non_overridable :: im_perform_array
    end type

    !> @brief Extends the interp_manager class allowing for linear, piecewise
    !! interpolation of a data set.
    type, extends(interp_manager) :: linear_interp
    contains
        !> @brief Performs the actual interpolation.
        procedure :: raw_interp => li_raw_interp
    end type

    !> @brief Extends the interp_manager class allowing for polynomial
    !! interpolation of a data set.
    type, extends(interp_manager) :: polynomial_interp
    private
        real(real64), allocatable, dimension(:) :: m_c
        real(real64), allocatable, dimension(:) :: m_d
        real(real64) :: m_dy
    contains
        !> @brief Initializes the polynomial_interp instance.
        procedure, public :: initialize => pi_init
        !> @brief Performs the actual interpolation.
        procedure :: raw_interp => pi_raw_interp
    end type

    !> @brief Extends the interp_manager class allowing for cubic spline
    !! interpolation of a data set.
    type, extends(interp_manager) :: spline_interp
    private
        real(real64), allocatable, dimension(:) :: m_ypp
    contains
        !> @brief Performs the actual interpolation.
        procedure :: raw_interp => si_raw_interp
        !> @brief Computes the second derivative terms for the cubic-spline
        !! model.
        procedure :: compute_diff2 => si_second_deriv
        !> @brief Initializes the spline_interp instance.
        procedure, public :: initialize => si_init_1
        !> @brief Initializes the spline_interp instance while allowing
        !! definition of boundary conditions.
        procedure, public :: initialize_spline => si_init_2
        !> @brief Interpolates to obtain the first derivative value at the
        !! specified independent variable.
        generic, public :: first_derivative => si_diff1, si_diff1_array
        !> @brief Interpolates to obtain the second derivative value at the
        !! specified independent variable.
        generic, public :: second_derivative => si_diff2, si_diff2_array

        procedure :: si_diff1
        procedure :: si_diff1_array
        procedure :: si_diff2
        procedure :: si_diff2_array
    end type

! ------------------------------------------------------------------------------
    interface
        !> @brief Defines the signature of a method used to interpolate a single
        !!  value in an X-Y data set.
        !!
        !! @param[in,out] this The interp_manager based instance.
        !! @param[in] jlo The array index below which @p pt is found in x.
        !! @param[in] pt The independent variable value to interpolate.
        !!
        !! @return The interpolated value.
        function interp_xy(this, jlo, pt) result(yy)
            use iso_fortran_env
            import interp_manager
            class(interp_manager), intent(inout) :: this
            integer(int32), intent(in) :: jlo
            real(real64), intent(in) :: pt
            real(real64) :: yy
        end function
    end interface

! ------------------------------------------------------------------------------
    ! INTERP_MANAGER ROUTINES
    interface
        !> @brief Initializes the specified interp_manager instance.
        !!
        !! @param[in,out] this The interp_manager instance.
        !! @param[in] x An N-element array containing the independent variable 
        !!  data.  The data in this array must be either monotonically 
        !!  increasing or decreasing.
        !! @param[in] y An N-element array containing the dependent variable 
        !!  data.
        !! @param[in] order The order of the interpolating polynomial.  Notice, 
        !!  this parameter is optional; however, if not specified, a default of 
        !!  1 is used.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
        !!      increasing or decreasing.
        module subroutine im_init(this, x, y, order, err)
            class(interp_manager), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x, y
            integer(int32), intent(in), optional :: order
            class(errors), intent(inout), optional, target :: err
        end subroutine

        !> @brief Attempts to locate the index in the array providing a lower 
        !!  bounds to the specified interpolation point.
        !!
        !! @param[in,out] this The interp_manager instance.
        !! @param[in] pt The interpolation point.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
        !!
        !! @return The array index below @p pt.
        module function im_locate(this, pt, err) result(j)
            class(interp_manager), intent(inout) :: this
            real(real64), intent(in) :: pt
            class(errors), intent(inout), optional, target :: err
            integer :: j
        end function

        !> @brief Attempts to locate the index in the array providing a lower 
        !!  bounds to the specified interpolation point.  This method is 
        !!  typically more efficient than locate when the current index does 
        !!  not stray too far from the previous.
        !!
        !! @param[in,out] this The interp_manager instance.
        !! @param[in] pt The interpolation point.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
        !!
        !! @return The array index below @p pt.
        module function im_hunt(this, pt, err) result(j)
            class(interp_manager), intent(inout) :: this
            real(real64), intent(in) :: pt
            class(errors), intent(inout), optional, target :: err
            integer(int32) :: j
        end function

        !> @brief Interpolates to obtain the function value at the specified
        !!  independent variable.
        !!
        !! @param[in,out] this The interp_manager instance.
        !! @param[in] pt The independent variable value to interpolate.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
        !!
        !! @return The interpolated value.
        module function im_perform(this, pt, err) result(yy)
            class(interp_manager), intent(inout) :: this
            real(real64), intent(in) :: pt
            class(errors), intent(inout), optional, target :: err
            real(real64) :: yy
        end function

        !> @brief Interpolates to obtain the function value at the specified
        !!  independent variables.
        !!
        !! @param[in,out] this The interp_manager instance.
        !! @param[in] pts An M-element array containing the independent variable
        !!  values to interpolate.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
        !!
        !! @return An M-element array containing the interpolated values.
        module function im_perform_array(this, pts, err) result(yy)
            class(interp_manager), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: pts
            class(errors), intent(inout), optional, target :: err
            real(real64), dimension(size(pts)) :: yy
        end function

        !> @brief Gets the number of stored data points.
        !!
        !! @param[in] this The interp_manager object.
        !!
        !! @return The number of data points.
        pure module function im_get_num_pts(this) result(n)
            class(interp_manager), intent(in) :: this
            integer(int32) :: n
        end function

        !> @brief Gets the x component of the requested data point.
        !!
        !! @param[in] this The interp_manager object.
        !! @param[in] ind The one-based index of the data point to retrieve.
        !!
        !! @return The x component of the requested data point.
        pure module function im_get_x(this, ind) result(x)
            class(interp_manager), intent(in) :: this
            integer(int32), intent(in) :: ind
            real(real64) :: x
        end function

        !> @brief Gets the y component of the requested data point.
        !!
        !! @param[in] this The interp_manager object.
        !! @param[in] ind The one-based index of the data point to retrieve.
        !!
        !! @return The y component of the requested data point.
        pure module function im_get_y(this, ind) result(y)
            class(interp_manager), intent(in) :: this
            integer(int32), intent(in) :: ind
            real(real64) :: y
        end function
    end interface

! ------------------------------------------------------------------------------
    ! LINAER_INTERP ROUTINES
    interface
        !> @brief Performs the actual linear interpolation.
        !!
        !! @param[in,out] this The linear_interp_mgr instance.
        !! @param[in] jlo The array index below which @p pt is found in x.
        !! @param[in] pt The independent variable value to interpolate.
        !!
        !! @return The interpolated value.
        module function li_raw_interp(this, jlo, pt) result(yy)
            class(linear_interp), intent(inout) :: this
            integer(int32), intent(in) :: jlo
            real(real64), intent(in) :: pt
            real(real64) :: yy
        end function
    end interface

! ------------------------------------------------------------------------------
    ! POLYNOMIAL_INTERP ROUTINES
    interface
        !> @brief Initializes the specified polynomial_interp instance.
        !!
        !! @param[in,out] this The polynomial_interp instance.
        !! @param[in] x An N-element array containing the independent variable 
        !!  data.  The data in this array must be either monotonically 
        !!  increasing or decreasing.
        !! @param[in] y An N-element array containing the dependent variable 
        !!  data.
        !! @param[in] order The order of the interpolating polynomial.  Notice, 
        !!  this parameter is optional; however, if not specified, a default 
        !!  of 1 is used.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_INVALID_INPUT_ERROR: Occurs if @p order is less than 1.
        !!  - M_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
        !!      increasing or decreasing.
        module subroutine pi_init(this, x, y, order, err)
            class(polynomial_interp), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x, y
            integer(int32), intent(in), optional :: order
            class(errors), intent(inout), optional, target :: err
        end subroutine

        !> @brief Performs the actual interpolation.
        !!
        !! @param[in,out] this The polynomial_interp instance.
        !! @param[in] jlo The array index below which @p pt is found in x.
        !! @param[in] pt The independent variable value to interpolate.
        !!
        !! @return The interpolated value.
        module function pi_raw_interp(this, jlo, pt) result(yy)
            class(polynomial_interp), intent(inout) :: this
            integer(int32), intent(in) :: jlo
            real(real64), intent(in) :: pt
            real(real64) :: yy
        end function
    end interface

! ------------------------------------------------------------------------------
    ! SPLINE_INTERP MEMBERS
    interface
        !> @brief Solves a pentadiagonal system of linear equations.  A
        !!  pentadiagonal matrix is all zeros with the exception of the 
        !!  diagonal, and the two immediate sub and super-diagonals.  The 
        !!  entries of row I are stored as follows:
        !!      A(I,I-2) -> A1(I)
        !!      A(I,I-1) -> A2(I)
        !!      A(I,I) -> A3(I)
        !!      A(I,I+1) -> A4(I)
        !!      A(I,I+2) -> A5(I)
        !!
        !! @param[in] a1 An N-element array as defined above.
        !! @param[in,out] a2 An N-element array as defined above.  This array is
        !!  overwritten by this routine during the solution process.
        !! @param[in,out] a3 An N-element array as defined above.  This array is
        !!  overwritten by this routine during the solution process.
        !! @param[in,out] a4 An N-element array as defined above.  This array is
        !!  overwritten by this routine during the solution process.
        !! @param[in] a5 An N-element array as defined above.
        !! @param[in,out] b An N-element array containing the right-hand-side.  
        !!  This array is overwritten by this routine during the solution 
        !!  process.
        !! @param[out] x An N-element array that, on output, contains the 
        !!  solution to the linear system.
        !!
        !! - [Spline Library](http://people.sc.fsu.edu/~jburkardt/f77_src/spline/spline.html)
        module subroutine penta_solve(a1, a2, a3, a4, a5, b, x)
            real(real64), intent(in), dimension(:) :: a1, a5
            real(real64), intent(inout), dimension(:) :: a2, a3, a4, b
            real(real64), intent(out), dimension(:) :: x
        end subroutine

        !> @brief Performs the actual interpolation.
        !!
        !! @param[in,out] this The spline_interp instance.
        !! @param[in] jlo The array index below which @p pt is found in x.
        !! @param[in] pt The independent variable value to interpolate.
        !!
        !! @return The interpolated value.
        module function si_raw_interp(this, jlo, pt) result(yy)
            class(spline_interp), intent(inout) :: this
            integer(int32), intent(in) :: jlo
            real(real64), intent(in) :: pt
            real(real64) :: yy
        end function

        !> @brief Computes the second derivative terms for the cubic-spline model.
        !!
        !! @param[in,out] this The spline_interp_mgr instance.
        !! @param[in] ibcbeg Defines the nature of the boundary condition at the
        !!  beginning of the spline.
        !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
        !!      initial interval.
        !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at 
        !!      its initial point is provided in @p ybcbeg.
        !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at 
        !!      its initial point is provided in @p ybcbeg.
        !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is 
        !!      continuous at x(2).
        !! @param[in] ybcbeg If needed, the value of the initial point boundary
        !!  condition.
        !! @param[in] ibcend Defines the nature of the boundary condition at the
        !!  end of the spline.
        !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
        !!      final interval.
        !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at 
        !!      its initial point is provided in @p ybcend.
        !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at 
        !!      its initial point is provided in @p ybcend.
        !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is 
        !!      continuous at x(n-1).
        !! @param[in] ybcend If needed, the value of the final point boundary
        !!  condition.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!
        !! @par Remarks
        !! This code is a slight modification of the SPLINE_CUBIC_SET routine 
        !! from the [SPLINE](http://people.sc.fsu.edu/~jburkardt/f77_src/spline/spline.html) 
        !! library.
        module subroutine si_second_deriv(this, ibcbeg, ybcbeg, ibcend, ybcend, err)
            class(spline_interp), intent(inout) :: this
            integer(int32), intent(in) :: ibcbeg, ibcend
            real(real64), intent(in) :: ybcbeg, ybcend
            class(errors), intent(inout), optional, target :: err
        end subroutine

        !> @brief Initializes the specified spline_interp instance.  The end 
        !!  points are considered free such that the interpolant is quadratic 
        !!  over both the initial and final intervals.
        !!
        !! @param[in,out] this The spline_interp instance.
        !! @param[in] x An N-element array containing the independent variable 
        !!  data.  The data in this array must be either monotonically 
        !!  increasing or decreasing.
        !! @param[in] y An N-element array containing the dependent variable 
        !!  data.
        !! @param[in] order The order of the interpolating polynomial.  This
        !!  parameter is ignored as the spline is a cubic approximation.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
        !!      increasing or decreasing.
        module subroutine si_init_1(this, x, y, order, err)
            class(spline_interp), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x, y
            integer(int32), intent(in), optional :: order
            class(errors), intent(inout), optional, target :: err
        end subroutine

        !> @brief Initializes the specified spline_interp instance.
        !!
        !! @param[in,out] this The spline_interp instance.
        !! @param[in] x An N-element array containing the independent variable 
        !!  data.  The data in this array must be either monotonically 
        !!  increasing or decreasing.
        !! @param[in] y An N-element array containing the dependent variable 
        !!  data.
        !! @param[in] ibcbeg An optional input that defines the nature of the
        !!  boundary condition at the beginning of the spline.  If no parameter,
        !!  or an invalid parameter, is specified, the default natural condition
        !!  (SPLINE_QUADRATIC_OVER_INTERVAL) is used.
        !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
        !!      initial interval.  No value is required for @p ybcbeg.
        !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at 
        !!      its initial point is provided in @p ybcbeg.
        !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at 
        !!      its initial point is provided in @p ybcbeg.
        !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is 
        !!      continuous at x(2).  No value is required for @p ybcbeg.
        !! @param[in] ybcbeg If needed, the value of the initial point boundary
        !!  condition.  If needed, but not supplied, a default value of zero 
        !!  will be used.
        !! @param[in] ibcend An optional input that defines the nature of the
        !!  boundary condition at the end of the spline.  If no parameter, or an
        !!  invalid parameter, is specified, the default natural condition
        !!  (SPLINE_QUADRATIC_OVER_INTERVAL) is used.
        !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
        !!      final interval.  No value is required for @p ybcend.
        !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at 
        !!      its initial point is provided in @p ybcend.
        !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at 
        !!      its initial point is provided in @p ybcend.
        !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is 
        !!      continuous at x(n-1).  No value is required for @p ybcend.
        !! @param[in] ybcend If needed, the value of the final point boundary
        !!  condition.  If needed, but not supplied, a default value of zero 
        !!  will be used.
        !!
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
        !!      increasing or decreasing.
        module subroutine si_init_2(this, x, y, ibcbeg, ybcbeg, ibcend, ybcend, err)
            class(spline_interp), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x, y
            integer(int32), intent(in), optional :: ibcbeg, ibcend
            real(real64), intent(in), optional :: ybcbeg, ybcend
            class(errors), intent(inout), optional, target :: err
        end subroutine

        !> @brief Interpolates to obtain the first derivative value at the 
        !! specified independent variable.
        !!
        !! @param[in,out] this The interp_manager instance.
        !! @param[in] pt The independent variable value to interpolate.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
        !!
        !! @return The interpolated value.
        module function si_diff1(this, pt, err) result(yy)
            ! Arguments
            class(spline_interp), intent(inout) :: this
            real(real64), intent(in) :: pt
            class(errors), intent(inout), optional, target :: err
            real(real64) :: yy
        end function

        !> @brief Interpolates to obtain the first derivative value at the 
        !!  specified independent variables.
        !!
        !! @param[in,out] this The interp_manager instance.
        !! @param[in] pts An M-element array containing the independent variable
        !!  values to interpolate.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
        !!
        !! @return An M-element array containing the interpolated values.
        module function si_diff1_array(this, pts, err) result(yy)
            ! Arguments
            class(spline_interp), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: pts
            class(errors), intent(inout), optional, target :: err
            real(real64), dimension(size(pts)) :: yy
        end function

        !> @brief Interpolates to obtain the second derivative value at the
        !! specified independent variable.
        !!
        !! @param[in,out] this The interp_manager instance.
        !! @param[in] pt The independent variable value to interpolate.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
        !!
        !! @return The interpolated value.
        module function si_diff2(this, pt, err) result(yy)
            ! Arguments
            class(spline_interp), intent(inout) :: this
            real(real64), intent(in) :: pt
            class(errors), intent(inout), optional, target :: err
            real(real64) :: yy
        end function

        !> @brief Interpolates to obtain the second derivative value at the
        !! specified independent variables.
        !!
        !! @param[in,out] this The interp_manager instance.
        !! @param[in] pts An M-element array containing the independent variable
        !!  values to interpolate.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
        !!
        !! @return An M-element array containing the interpolated values.
        module function si_diff2_array(this, pts, err) result(yy)
            ! Arguments
            class(spline_interp), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: pts
            class(errors), intent(inout), optional, target :: err
            real(real64), dimension(size(pts)) :: yy
        end function
    end interface

! ******************************************************************************
! MEASUREMENTS_SMOOTHING.F90
! ------------------------------------------------------------------------------
    !> @brief Defines a type for computing a smoothing of an X-Y data set using
    !! a robust locally weighted scatterplot smoothing (LOWESS) algorithm.
    type lowess_smoothing
        private
        !> N-element array of x data points - sorted into ascending order.
        real(real64), allocatable, dimension(:) :: m_x
        !> N-element array of y data points.
        real(real64), allocatable, dimension(:) :: m_y
        !> N-element array containing the robustness weights for each data
        !! point.
        real(real64), allocatable, dimension(:) :: m_weights
        !> N-element array containing the residuals (Y - YS)
        real(real64), allocatable, dimension(:) :: m_residuals
        !> Scaling parameter used to define the nature of the linear
        !! interpolations used by the algorithm.
        real(real64) :: m_delta
        !> Tracks whether or not ls_init has been called
        logical :: m_init = .false.
    contains
        !> @brief Initializes the lowess_smoothing object.
        procedure, public :: initialize => ls_init
        !> @brief Performs the actual smoothing operation.
        procedure, public :: smooth => ls_smooth
        !> @brief Gets the number of stored data points.
        procedure, public :: get_count => ls_get_num_pts
        !> @brief Gets the x component of the requested data point.
        procedure, public :: get_x => ls_get_x
        !> @brief Gets the y component of the requested data point.
        procedure, public :: get_y => ls_get_y
        !> @brief Gets the residuals from each data point.
        procedure, public :: get_residuals => ls_get_residual
    end type

! ------------------------------------------------------------------------------
    ! LOWESS_SMOOTHING ROUTINES
    interface
        !> @brief Initializes the lowess_smoothing object.
        !!
        !! @param[in,out] this The lowess_smoothing object.
        !! @param[in] x An N-element containing the independent variable values 
        !!  of the data set.  This array must be in a monotonically increasing 
        !!  order.  The routine is capable of sorting the array into ascending 
        !!  order, dependent upon the value of @p srt.  If sorting is performed,
        !!  this routine will also shuffle @p y to match.
        !! @param[in] y  An N-element array of the dependent variables 
        !!  corresponding to @p x.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
        !!      increasing or decreasing.
        module subroutine ls_init(this, x, y, err)
            class(lowess_smoothing), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x, y
            class(errors), intent(inout), optional, target :: err
        end subroutine

        !> @brief Performs the actual smoothing operation.
        !!
        !! @param[in,out] this The lowess_smoothing object.
        !! @param[in] f Specifies the amount of smoothing.  More specifically, 
        !!  this value is the fraction of points used to compute each value.  
        !!  As this value increases, the output becomes smoother.  Choosing a 
        !!  value in the range of 0.2 to 0.8 usually results in a good fit.  As 
        !!  such, a reasonable starting point, in the absence of better 
        !!  information, is a value of 0.5.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_NO_DATA_DEFINED_ERROR: Occurs if no data has been defined.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!
        !! @return The smoothed data points.
        module function ls_smooth(this, f, err) result(ys)
            class(lowess_smoothing), intent(inout) :: this
            real(real64), intent(in) :: f
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: ys
        end function

        !> @brief Gets the number of stored data points.
        !!
        !! @param[in] this The lowess_smoothing object.
        !!
        !! @return The number of data points.
        pure module function ls_get_num_pts(this) result(n)
            class(lowess_smoothing), intent(in) :: this
            integer(int32) :: n
        end function

        !> @brief Gets the x component of the requested data point.
        !!
        !! @param[in] this The lowess_smoothing object.
        !! @param[in] ind The one-based index of the data point to retrieve.
        !!
        !! @return The x component of the requested data point.
        pure module function ls_get_x(this, ind) result(x)
            class(lowess_smoothing), intent(in) :: this
            integer(int32), intent(in) :: ind
            real(real64) :: x
        end function

        !> @brief Gets the y component of the requested data point.
        !!
        !! @param[in] this The lowess_smoothing object.
        !! @param[in] ind The one-based index of the data point to retrieve.
        !!
        !! @return The y component of the requested data point.
        pure module function ls_get_y(this, ind) result(y)
            class(lowess_smoothing), intent(in) :: this
            integer(int32), intent(in) :: ind
            real(real64) :: y
        end function

        !> @brief Gets the residuals from each data point.
        !!
        !! @param[in] this The lowess_smoothing object.
        !! @param[out] x An N-element array where the residual data should be
        !!  written.
        module subroutine ls_get_residual(this, x)
            class(lowess_smoothing), intent(in) :: this
            real(real64), intent(out), dimension(:) :: x
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    ! ADDITIONAL SMOOTHING ROUTINES
    interface
        !> @brief Applies a moving average to smooth a data set.
        !!
        !! @param[in,out] x On input, the signal to smooth.  On output, the 
        !!  smoothed signal.
        !! @param[in] npts The size of the averaging window.  This value must be
        !!  at least 2, but no more than the number of elements in @p x.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_INVALID_INPUT_ERROR: Occurs if @p npts is less than 2, or 
        !!      greater than the length of @p x.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        module subroutine moving_average(x, npts, err)
            real(real64), intent(inout), dimension(:) :: x
            integer(int32), intent(in) :: npts
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ******************************************************************************
! MEASUREMENTS_REGRESSION.F90
! ------------------------------------------------------------------------------
    interface
        !> @brief Fits the multiple input, multiple output linear model
        !! A * X = Y by solving for matrix A in a least-squares sense.
        !!
        !! @param[in] x An N-by-K matrix of known independent variables.  K
        !!  must be greater than or equal to N.
        !! @param[in] y An M-by-K matrix of known dependent variables.  Notice,
        !!  M must be less than or equal to N, and K must be greater than or
        !!  equal to M.
        !! @param[in] thrsh  An optional input that is used to determine the
        !!  closeness of zero for the singular values of matrix @p X.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_ARRAY_SIZE_ERROR: Occurs if there is a size-mismatch in the
        !!      matrix equation.
        !!  - M_UNDERDEFINED_PROBLEM: Occurs if there is insufficient data 
        !!      (e.g. k < n), or the problem is not sized appropriately 
        !!      (e.g. m > n).
        !!
        !! @return The M-by-N coefficient matrix A.
        !!
        !! @remarks
        !! Solving the linear system is straight-forward when M <= N by means
        !! of singular value decomposition.  Specifically, the Moore-Penrose
        !! pseudo-inverse utilizing singular value decomposition.  The
        !! solution is obtained as follows.
        !! @par
        !! \f$ A X X^{+} = Y X^{+} \f$
        !! @par
        !! \f$ X X^{+} = I \f$ as X is underdetermined.  Then
        !! @par
        !! \f$ A = Y X^{+} \f$
        !! @par
        !! where \f$ X^{+} \f$ is the Moore-Penrose pseudo-inverse of X.
        module function linear_least_squares_mimo(x, y, thrsh, err) result(a)
            real(real64), intent(in), dimension(:,:) :: x, y
            real(real64), intent(in), optional :: thrsh
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:,:) :: a
        end function

        !> @brief Fits the multiple input, single output linear model
        !! A * X = Y by solving for array A in a least-squares sense.
        !!
        !! @param[in] x An N-by-K matrix of known independent variables.  K
        !!  must be greater than or equal to N.
        !! @param[in] y A K-element array containing the known dependent 
        !!  variables.
        !! @param[in] thrsh  An optional input that is used to determine the
        !!  closeness of zero for the singular values of matrix @p X.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_ARRAY_SIZE_ERROR: Occurs if there is a size-mismatch in the
        !!      matrix equation.
        !!  - M_UNDERDEFINED_PROBLEM: Occurs if there is insufficient data 
        !!      (e.g. k < n), or the problem is not sized appropriately 
        !!      (e.g. size(y) /= n).
        !!
        !! @return The N-element coefficient array A.
        !!
        !! @remarks
        !! Solving the linear system is straight-forward by means of
        !! singular value decomposition.  Specifically, the Moore-Penrose
        !! pseudo-inverse utilizing singular value decomposition.  The
        !! solution is obtained as follows.
        !! @par
        !! \f$ A X X^{+} = Y X^{+} \f$
        !! @par
        !! \f$ X X^{+} = I \f$ as X is underdetermined.  Then
        !! @par
        !! \f$ A = Y X^{+} \f$
        !! @par
        !! where \f$ X^{+} \f$ is the Moore-Penrose pseudo-inverse of X.
        module function linear_least_squares_miso(x, y, thrsh, err) result(a)
            real(real64), intent(in), dimension(:,:) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(in), optional :: thrsh
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: a
        end function
    end interface

! ******************************************************************************
! MEASUREMENTS_PEAK.F90
! ------------------------------------------------------------------------------
    !> @brief Describes peak and valley values and associated locations within 
    !! the data set.
    type peak_info
        !> @brief A list of the maximum (peak) values.
        real(real64), allocatable, dimension(:) :: max_values
        !> @brief A list of the indices of the maximum values.
        integer(int32), allocatable, dimension(:) :: max_value_indices
        !> @brief A list of the minimum (valley) values.
        real(real64), allocatable, dimension(:) :: min_values
        !> @brief A list of the indices of the minimum values.
        integer(int32), allocatable, dimension(:) :: min_value_indices
    end type

! ------------------------------------------------------------------------------
    interface
        !> @brief Attempts to locate all peaks and valleys within the given
        !! data set bound by the constraints specified.
        !!
        !! @param[in] v The data set.
        !! @param[in] delta A threshold level used to denote the minimum change
        !!  acceptable in determining a peak or valley from neighboring points.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!
        !! @return A peak_info type containing the results of the search.
        !!
        !! @par References
        !! - http://billauer.co.il/peakdet.html
        module function peak_detect(v, delta, err) result(rst)
            real(real64), intent(in), dimension(:) :: v
            real(real64), intent(in) :: delta
            class(errors), intent(inout), optional, target :: err
            type(peak_info) :: rst
        end function
    end interface

! ******************************************************************************
! MEASUREMENT_SPECTRUM.F90
! ------------------------------------------------------------------------------
    interface
        !> @brief Defines a window function.
        !!
        !! @param[in] bin The index or bin number (0 <= @p bin <= @p n).
        !! @param[in] n The transform length.
        !! @return The window function value.
        function window_function(bin, n) result(x)
            use iso_fortran_env
            integer(int32), intent(in) :: bin, n
            real(real64) :: x
        end function
    end interface

! ------------------------------------------------------------------------------
    interface
        !> @brief Tests to see if a number is an integer power of two.
        !!
        !! @param[in] n The integer to test.
        !! @return Returns true if @p n is a power of two; else, false.
        pure elemental module function is_power_of_two(n) result(rst)
            integer(int32), intent(in) :: n
            logical :: rst
        end function

        !> @brief Provides the next higher integer power of two.
        !!
        !! @param[in] x The value to test.
        !! @return The next power of two higher than @p x.  If @p x is already
        !! a power of two, it's value is simply returned.  For instance, if @p
        !! is set to 128, then a value of 7 is returned.  However, if a value
        !! of 129 is supplied, then a value of 8 is returned.
        pure elemental module function next_power_of_two(x) result(n)
            integer(int32), intent(in) :: x
            integer(int32) :: n
        end function

        !> @brief Pads an array with zeros to arrive at an array whose length
        !! is the next integer power of two higher.  If the array is already
        !! an integer power of two in length no action is taken.
        !!
        !! @param[in] x The array to pad.
        !! @return The padded array.
        pure module function pad_with_zeros(x) result(xp)
            real(real64), intent(in), dimension(:) :: x
            real(real64), allocatable, dimension(:) :: xp
        end function

        !> @brief Computes the Fourier transform of a discretely sampled signal.
        !!
        !! @param[in] x An N-element array containing the signal.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!
        !! @return The positive frequency side of the complex-valued Fourier 
        !!  transform of @p x.
        module function fourier_transform(x, err) result(f)
            real(real64), intent(in), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
            complex(real64), allocatable, dimension(:) :: f
        end function

        !> @brief Computes the periodogram power spectral density (PSD) estimate
        !! of a signal.                                                                                                                                                             
        !!
        !! @param[in] x An N-element array containing the signal.
        !! @param[in] winfun The window function to apply.
        !! @param[in] nfft The length Fourier transform to apply.  This must
        !!  be an integer power of two, even if @p x is not an integer power
        !!  of two in length.  If this parameter is less than the length of
        !!  @p x, @p x will be overlapped, windowed, and averaged to achieve
        !!  an estimate of the power spectrum.  If this parameter is larger
        !!  than the length of @p x, @p x will be padded with zeros prior
        !!  to windowing.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_INVALID_INPUT_ERROR: Occurs if @p nfft is not an integer power
        !!      of two.
        !!
        !! @return The periodogram (power spectrum) of @p x.
        module function periodogram(x, winfun, nfft, err) result(p)
            real(real64), intent(in), dimension(:) :: x
            procedure(window_function), pointer, intent(in) :: winfun
            integer(int32), intent(in) :: nfft
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: p
        end function

        !> @brief Computes a suitable frequency for a Fourier transformed data
        !! set.
        !!
        !! @param[in] fs The rate at which the data was sampled.  Notice, the
        !!  returned frequency value will be expressed in the same units as this
        !!  value.
        !! @param[in] i The frequency bin such that 0 <= i < m where m is
        !!  @p nxfrm / 2 + 1 if @p nxfrm is even; else, (@p nxfrm + 1) / 2 if
        !!  @p nxfrm is odd.
        !! @param[in] nxfrm The length of the signal that was transformed.
        !!
        !! @return The frequency value.
        pure elemental module function fourier_frequency(fs, i, nxfrm) result(f)
            real(real64), intent(in) :: fs
            integer(int32), intent(in) :: i, nxfrm
            real(real64) :: f
        end function

        !> @brief Defines a rectangular window.
        !!
        !! @param[in] j The index or bin number (0 <= @p bin <= @p n).
        !! @param[in] n The transform length.
        !!
        !! @return The value of the window function at index @p j.
        pure module function rectangular_window(j, n) result(x)
            integer(int32), intent(in) :: j, n
            real(real64) :: x
        end function

        !> @brief Defines a Hann window.
        !!
        !! @param[in] j The index or bin number (0 <= @p bin <= @p n).
        !! @param[in] n The transform length.
        !!
        !! @return The value of the window function at index @p j.
        pure module function hann_window(j, n) result(x)
            integer(int32), intent(in) :: j, n
            real(real64) :: x
        end function

        !> @brief Defines a Hamming window.
        !!
        !! @param[in] j The index or bin number (0 <= @p bin <= @p n).
        !! @param[in] n The transform length.
        !!
        !! @return The value of the window function at index @p j.
        pure module function hamming_window(j, n) result(x)
            integer(int32), intent(in) :: j, n
            real(real64) :: x
        end function

        !> @brief Defines a Welch window.
        !!
        !! @param[in] j The index or bin number (0 <= @p bin <= @p n).
        !! @param[in] n The transform length.
        !!
        !! @return The value of the window function at index @p j.
        pure module function welch_window(j, n) result(x)
            integer(int32), intent(in) :: j, n
            real(real64) :: x
        end function

        !> @brief Defines a Blackman-Harris window.
        !!
        !! @param[in] j The index or bin number (0 <= @p bin <= @p n).
        !! @param[in] n The transform length.
        !!
        !! @return The value of the window function at index @p j.
        pure module function blackman_harris_window(j, n) result(x)
            integer(int32), intent(in) :: j, n
            real(real64) :: x
        end function

        !> @brief Computes an FFT of a data set.  The results of the transform
        !! are normalized such that an inverse transform will result in the 
        !! original signal.
        !!
        !! @param[in] x An N-element array containing the data to transform.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!
        !! @return The complex-valued Fourier transform of @p x.
        module function fft(x, err) result(f)
            real(real64), intent(in), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
            complex(real64), allocatable, dimension(:) :: f
        end function

        !> @brief Computes the inverse FFT of a data set.
        !!
        !! @param[in] x An N-element array containing the data to transform.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!
        !! @return The resulting inverse Fourier transform of @p x.
        module function ifft(x, err) result(f)
            complex(real64), intent(in), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
            complex(real64), allocatable, dimension(:) :: f
        end function
    end interface

! ******************************************************************************
! MEASUREMENTS_FILTER.F90
! ------------------------------------------------------------------------------
    interface
        !> @brief Applies a low-pass filter to a signal.
        !!
        !! @param[in] x The array containing the signal to filter.
        !! @param[in] fs The frequency at which @p x was sampled.
        !! @param[in] cutoff The cut-off frequency.  This value must be
        !!  positive-valued, and must be less than the Nyquist frequency.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff is either not positive
        !!      valued, or is greater than or equal to the Nyquist frequency.
        !!
        !! @return The filtered signal.
        module function low_pass_filter(x, fs, cutoff, err) result(y)
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(in) :: fs, cutoff
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: y
        end function

        !> @brief Applies a high-pass filter to a signal.
        !!
        !! @param[in] x The array containing the signal to filter.
        !! @param[in] fs The frequency at which @p x was sampled.
        !! @param[in] cutoff The cut-off frequency.  This value must be
        !!  positive-valued, and must be less than the Nyquist frequency.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff is either not positive
        !!      valued, or is greater than or equal to the Nyquist frequency.
        !!
        !! @return The filtered signal.
        module function high_pass_filter(x, fs, cutoff, err) result(y)
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(in) :: fs, cutoff
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: y
        end function

        !> @brief Applies a band-pass filter to a signal.
        !!
        !! @param[in] x The array containing the signal to filter.
        !! @param[in] fs The frequency at which @p x was sampled.
        !! @param[in] cutoff1 The lower cut-off frequency.  This value must be
        !!  positive-valued, and must be less than the Nyquist frequency.
        !! @param[in] cutoff2 The upper cut-off frequency.  This value must be
        !!  positive-valued, and must be less than the Nyquist frequency.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff1 or @p cutoff2 is 
        !!      either not positive-valued, or is greater than or equal to the 
        !!      Nyquist frequency.
        !!
        !! @return The filtered signal.
        module function band_pass_filter(x, fs, cutoff1, cutoff2, err) result(y)
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(in) :: fs, cutoff1, cutoff2
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: y
        end function

        !> @brief Applies a band-stop filter to a signal.
        !!
        !! @param[in] x The array containing the signal to filter.
        !! @param[in] fs The frequency at which @p x was sampled.
        !! @param[in] cutoff1 The lower cut-off frequency.  This value must be
        !!  positive-valued, and must be less than the Nyquist frequency.
        !! @param[in] cutoff2 The upper cut-off frequency.  This value must be
        !!  positive-valued, and must be less than the Nyquist frequency.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff1 or @p cutoff2 is 
        !!      either not positive-valued, or is greater than or equal to the 
        !!      Nyquist frequency.
        !!
        !! @return The filtered signal.
        module function band_stop_filter(x, fs, cutoff1, cutoff2, err) result(y)
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(in) :: fs, cutoff1, cutoff2
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: y
        end function
    end interface

! ******************************************************************************
! MEASUREMENTS_ARRAY.F90
! ------------------------------------------------------------------------------
    interface
        !> @brief Computes the cumulative sum of an array.
        !!
        !! @param[in] x The N-element array on which to operate.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!
        !! @return The resulting N-element array.
        module function cumulative_sum(x, err) result(y)
            real(real64), intent(in), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: y
        end function

        !> @brief Computes the differences between each element in the array.
        !!
        !! @param[in] x The N-element array on which to operate.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!
        !! @return The N-1 element array containing the differences between 
        !!  each element in @p x.
        module function difference(x, err) result(dx)
            real(real64), intent(in), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: dx
        end function

        !> @brief Unwraps an array of phase angles.
        !!
        !! @param[in] x The array to unwrap.
        !! @param[in] cutoff  An optional input that specifies the threshold
        !!  value to use when unwrapping steps in @p x.  The default is pi.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!
        !! @return The unwraped phase array.
        module function unwrap(x, cutoff, err) result(p)
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(in), optional :: cutoff
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: p
        end function

        !> @brief Estimates the derivative of a discretely sampled signal by
        !! means of finite differences.
        !!
        !! @param[in] x The N-element array containing the independent variable
        !!  data.
        !! @parma[in] y The N-element array containing the dependent variable
        !!  data.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
        !!
        !! @return The N-element derivative estimate.
        module function finite_difference(x, y, err) result(dydx)
            real(real64), intent(in), dimension(:) :: x, y
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: dydx
        end function

        !> @brief Estimates the indefinite integral of a discretely sampled
        !! signal.
        !!
        !! @param[in] x The N-element array containing the independent variable
        !!  data.
        !! @parma[in] y The N-element array containing the dependent variable
        !!  data.
        !! @param[in] c An optional input defining the initial condition.  The
        !!  default is 0.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
        !!
        !! @return The N-element integral estimate.
        module function integrate(x, y, c, err) result(f)
            real(real64), intent(in), dimension(:) :: x, y
            real(real64), intent(in), optional :: c
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: f
        end function

        !> @brief Estimates the definite integral of a discretely sampled signal
        !! using a trapezoidal approach.
        !!
        !! @param[in] The N-element array containing the independent variable
        !!  data.
        !! @parma[in] y The N-element array containing the dependent variable
        !!  data.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
        !!
        !! @return The result of the integration.
        module function trapz_integrate(x, y, err) result(f)
            real(real64), intent(in), dimension(:) :: x, y
            class(errors), intent(inout), optional, target :: err
            real(real64) :: f
        end function

        !> @brief Removes the DC offset from a signal.
        !!
        !! @param[in] x The signal on which to operate.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!
        !! @return The signal with the DC offset removed.
        module function remove_dc_offset(x, err) result(y)
            real(real64), intent(in), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:) :: y
        end function
    end interface

! ******************************************************************************
! MEASUREMENTS_ANOVA.F90
! ------------------------------------------------------------------------------
    !> @brief A single entry in an ANOVA table.
    type anova_table_entry
        !> @brief The number of degrees of freedom.
        integer(int32) :: dof
        !> @brief The sum of the squares.
        real(real64) :: sum_of_squares
        !> @brief The mean of the squares.
        real(real64) :: mean_of_squares
        !> @brief The F statistic.
        real(real64) :: f_stat
        !> @brief The mean value.
        real(real64) :: mean
        !> @brief The probability (p-value) that the null hypothesis is true.
        real(real64) :: probability
    end type

! ------------------------------------------------------------------------------
    !> @brief A gage repeatablilty and reproducibility (GR&R) ANOVA table.
    type gage_anova_table
        !> @brief The individual operator information.
        type(anova_table_entry), allocatable, dimension(:) :: operator
        !> @brief The combined operator information.
        type(anova_table_entry) :: operators
        !> @brief The individual part information.
        type(anova_table_entry), allocatable, dimension(:) :: part
        !> @brief The combined part information
        type(anova_table_entry) :: parts
        !> @brief The operator-by-part information
        type(anova_table_entry) :: operator_by_part
        !> @brief The measurement equipment information.
        type(anova_table_entry) :: equipment
        !> @brief The total variability information.
        type(anova_table_entry) :: total
    end type

! ------------------------------------------------------------------------------
    interface
        !> @brief Computes a crossed-effects analysis of variance (ANOVA) for
        !! a set of measurement data in order to better understand variance
        !! contributions of each part of a measurement system or gage analysis.
        !!
        !! @param[in] x An M-by-N-by-P matrix containing the measurement results
        !! taken by each operator where M is the number of parts tested by each
        !! operator, N is the number of times each operator tested each part, 
        !! and P is the number of test operators or inspectors.  Each value
        !! must be greater than one.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution.  If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - M_INSUFFICIENT_DATA_ERROR:
        !!
        !! @return A @p gage_anova_table object containing the variance results.
        !!
        !! @par Example
        !! The following example illustrates an analysis of variance example
        !! worked from https://www.spcforexcel.com/knowledge/measurement-systems-analysis/anova-gage-rr-part-1.
        !! @code{.f90}
        !! program main
        !!     use iso_fortran_env
        !!     use measurements_core
        !!     implicit none
        !!
        !!     ! Parameters
        !!     integer(int32), parameter :: nops = 3
        !!     integer(int32), parameter :: nparts = 5
        !!     integer(int32), parameter :: ntrials = 3
        !!
        !!     ! Local Variables
        !!     real(real64) :: x(nparts, ntrials, nops)
        !!     type(gage_anova_table) :: rst
        !!
        !!     ! Operator 1 Results
        !!     x(:,:,1) = reshape([&
        !!         3.29d0, 2.44d0, 4.34d0, 3.47d0, 2.2d0, &
        !!         3.41d0, 2.32d0, 4.17d0, 3.5d0, 2.08d0, &
        !!         3.64d0, 2.42d0, 4.27d0, 3.64d0, 2.16d0], &
        !!         [nparts, ntrials])
        !!
        !!     ! Operator 2 Results
        !!     x(:,:,2) = reshape([ &
        !!         3.08d0, 2.53d0, 4.19d0, 3.01d0, 2.44d0, &
        !!         3.25d0, 1.78d0, 3.94d0, 4.03d0, 1.8d0, &
        !!         3.07d0, 2.32d0, 4.34d0, 3.2d0, 1.72d0], &
        !!         [nparts, ntrials])
        !!
        !!     ! Operator 3 Results
        !!     x(:,:,3) = reshape([ &
        !!         3.04d0, 1.62d0, 3.88d0, 3.14d0, 1.54d0, &
        !!         2.89d0, 1.87d0, 4.09d0, 3.2d0, 1.93d0, &
        !!         2.85d0, 2.04d0, 3.67d0, 3.11d0, 1.55d0], &
        !!         [nparts, ntrials])
        !!
        !!     ! Perform the ANOVA
        !!     rst = gage_anova(x)
        !!
        !!     ! Display the results
        !!     print '(A)', "Operator Results:"
        !!     print '(AI0)', achar(9) // "DOF: ", rst%operators%dof
        !!     print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%operators%sum_of_squares
        !!     print '(AF0.3)', achar(9) // "Mean of Squares (Variance): ", rst%operators%mean_of_squares
        !!     print '(AF0.3)', achar(9) // "F Statistic: ", rst%operators%f_stat
        !!     print '(AF0.5)', achar(9) // "Probability: ", rst%operators%probability
        !!
        !!     print '(A)', new_line('a') // "Part Results:"
        !!     print '(AI0)', achar(9) // "DOF: ", rst%parts%dof
        !!     print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%parts%sum_of_squares
        !!     print '(AF0.3)', achar(9) // "Mean of Squares (Variance): ", rst%parts%mean_of_squares
        !!     print '(AF0.3)', achar(9) // "F Statistic: ", rst%parts%f_stat
        !!     print '(AF0.5)', achar(9) // "Probability: ", rst%parts%probability
        !!
        !!     print '(A)', new_line('a') // "Equipment Results:"
        !!     print '(AI0)', achar(9) // "DOF: ", rst%equipment%dof
        !!     print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%equipment%sum_of_squares
        !!     print '(AF0.3)', achar(9) // "Mean of Squares (Variance): ", rst%equipment%mean_of_squares
        !!
        !!     print '(A)', new_line('a') // "Operator-Part Interaction Results:"
        !!     print '(AI0)', achar(9) // "DOF: ", rst%operator_by_part%dof
        !!     print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%operator_by_part%sum_of_squares
        !!     print '(AF0.3)', achar(9) // "Mean of Squares (Variance): ", rst%operator_by_part%mean_of_squares
        !!     print '(AF0.3)', achar(9) // "F Statistic: ", rst%operator_by_part%f_stat
        !!     print '(AF0.5)', achar(9) // "Probability: ", rst%operator_by_part%probability
        !!
        !!     print '(A)', new_line('a') // "Total Results:"
        !!     print '(AI0)', achar(9) // "DOF: ", rst%total%dof
        !!     print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%total%sum_of_squares
        !!     print '(AF0.3)', achar(9) // "Mean of Squares (Variance): ", rst%total%mean_of_squares
        !!     print '(AF0.3)', achar(9) // "Overall Mean: ", rst%total%mean
        !! end program
        !! @endcode
        !! @code{.txt}
        !! Operator Results:
        !!         DOF: 2
        !!         Sum of Squares: 1.630
        !!         Mean of Squares (Variance): .815
        !!         F Statistic: 100.322
        !!         Probability: .00000
        !!
        !! Part Results:
        !!         DOF: 4
        !!         Sum of Squares: 28.909
        !!         Mean of Squares (Variance): 7.227
        !!         F Statistic: 889.458
        !!         Probability: -.00000
        !!
        !! Equipment Results:
        !!         DOF: 30
        !!         Sum of Squares: 1.712
        !!         Mean of Squares (Variance): .057
        !!
        !! Operator-Part Interaction Results:
        !!         DOF: 8
        !!         Sum of Squares: .065
        !!         Mean of Squares (Variance): .008
        !!         F Statistic: .142
        !!         Probability: .99637
        !!
        !! Total Results:
        !!         DOF: 44
        !!         Sum of Squares: 32.317
        !!         Mean of Squares (Variance): .734
        !!         Overall Mean: 2.944
        !! @endcode
        module function gage_anova(x, err) result(rst)
            real(real64), intent(in), target, contiguous, dimension(:,:,:) :: x
            class(errors), intent(inout), optional, target :: err
            type(gage_anova_table) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
end module
