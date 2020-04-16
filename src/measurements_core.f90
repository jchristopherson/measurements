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

    !> @brief A type containing variance components describing a measurement
    !! process.
    type, bind(C) :: process_variance
        !> @brief The measurement variation component.  In a gauge analysis this
        !! is referred to as the gauge R&R.  It is the sum of the repeatability
        !! and reproducibility variance components.
        real(real64) :: measurement_variance
        !> @brief The part variance component.
        real(real64) :: part_variance
        !> @brief The total process variance.  This is the sum of the 
        !! measurement variance and part variance.
        real(real64) :: total_variance
        !> @brief The equipment variance component.  This is often referred to
        !! as the repeatability component.
        real(real64) :: equipment_variance
        !> @brief The operator variance component.
        real(real64) :: operator_variance
        !> @brief The operator by part variance component.
        real(real64) :: operator_by_part_variance
    end type

    !> @brief A type containing gauge R&R results.
    type, bind(C) :: grr_results
        !> @brief The precision to tolerance ratio (P/T ratio).  This ratio
        !! is simply the ratio of the measurement standard deviation to
        !! the tolerance range.
        real(real64) :: pt_ratio
        !> @brief The precision to total variation ratio (P/TV ratio).  This 
        !! ratio is simply the ratio of the measurement standard deviation to
        !! the total process standard deviation.
        real(real64) :: ptv_ratio
        !> @brief The tolerance range.
        real(real64) :: tolerance_range
    end type

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

        !> @brief Utilizes an analysis of variance (ANOVA) to determine the
        !! variance components of a measurement process represented by the
        !! supplied data set.
        !!
        !! @param[in] x An M-by-N-by-P data set from the measurement process
        !!  to analyze where M is the number of parts tested (must be greater
        !!  than 1), N is the number of tests performed per part (must be
        !!  greater than 1), and P is the number of operators performing the
        !!  tests (must be at least 1).
        !! @param[in] alpha An optional parameter used to determine the 
        !!  appropriate calculation path.  The default value is 0.05.
        !!
        !! @return The resulting variance components of the process.
        !!
        !! @remarks It is possible for this routine to return zero-valued
        !!  variance components.  In such an event it is recommended that
        !!  another variance estimator is utilized.
        pure module function anova(x, alpha) result(rst)
            real(real64), intent(in), dimension(:,:,:) :: x
            real(real64), intent(in), optional :: alpha
            type(process_variance) :: rst
        end function

        !> @brief Utilizes a control chart type approach to evaluate the 
        !! measurement process utilized to collect the supplied data set.
        !!
        !! @param[in] x An M-by-N-by-P data set from the measurement process
        !!  to analyze where M is the number of parts tested (must be greater
        !!  than 1), N is the number of tests performed per part (must be
        !!  greater than 1), and P is the number of operators performing the
        !!  tests (must be at least 1).
        !!
        !! @return The resulting variance components of the process.
        !!
        !! @remarks This approach is best suited for smaller data sets whose
        !!  dimension doesn't exceed 25-30.  Anything over this size is better
        !!  served by another technique.
        pure module function control_chart_variance(x) result(rst)
            real(real64), intent(in), dimension(:,:,:) :: x
            type(process_variance) :: rst
        end function

        !> @brief Computes the gauge R&R statistics given a supplied set of
        !! process variance data.
        !!
        !! @param[in] k The multiplier to use when computing the P/T ratio.
        !!  Typically, a value of 6 is used for this factor.
        !! @param[in] x The process variance data.
        !! @param[in] usl The upper specification limit.
        !! @param[in] lsl The lower specification limit.
        !!
        !! @return The gauge R&R statistics.
        pure module function compute_grr(k, x, usl, lsl) result(rst)
            real(real64), intent(in) :: k, usl, lsl
            type(process_variance), intent(in) :: x
            type(grr_results) :: rst
        end function

        !> @brief Computes the discrimination ratio.
        !!
        !! @param[in] tv The total variance.
        !! @param[in] mv The measurement system variance.
        !!
        !! @return The results of the operation.
        !!
        !! @par 
        !! The discrimination ratio is computed as follows.
        !! @par
        !! /f$ DR = \sqrt{\frac{2 \sigma_{total}^2}{\sigma_{meas}^2} - 1} /f$
        !! @par
        !! An alternate means of computing this parameter (as used in JMP)
        !! is as follows.
        !! @par
        !! /f$ DR = 1.41 \frac{\sigma_{parts}}{\sigma_{meas}} /f$
        pure elemental module function discrimination_ratio(tv, mv) result(x)
            real(real64), intent(in) :: tv, mv
            real(real64) :: x
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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
