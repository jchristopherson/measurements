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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
