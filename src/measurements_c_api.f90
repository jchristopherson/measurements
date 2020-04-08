! measurements_c_api.f90

!> @brief A module providing C bindings to the MEASUREMENTS library.
module measurements_c_api
    use iso_c_binding
    use iso_fortran_env
    use measurements_core
    use ferror
    implicit none
contains
! ------------------------------------------------------------------------------
    !> @brief Tests to see if an array is monotonically increasing or 
    !! decreasing.
    !!
    !! @param[in] n The number of data points in the array.
    !! @param[in] x The N-element array to test.
    !! 
    !! @return Returns true if the array is monotonically increasing or
    !! monotonically decreasing; else, returns false.
    function c_is_monotonic(n, x) bind(C, name = "c_is_monotonic") result(rst)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        logical(c_bool) :: rst

        ! Process
        rst = logical(is_monotonic(x), c_bool)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the mean of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the data set.
    !!
    !! @return The mean of @p x.
    function c_mean(n, x) bind(C, name = "c_mean") result(z)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double) :: z

        ! Process
        z = mean(x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the median of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the data set.
    !!
    !! @return The median of @p x.
    function c_median(n, x) bind(C, name = "c_median") result(z)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double) :: z

        ! Process
        z = median(x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the sample variance of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the data set.
    !!
    !! @par Remarks
    !! To avoid overflow-type issues, Welford's algorithm is employed.  A 
    !! simple illustration of this algorithm can be found
    !! [here](https://www.johndcook.com/blog/standard_deviation/).
    function c_variance(n, x) bind(C, name = "c_variance") result(v)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double) :: v

        ! Process
        v = variance(x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the standard deviation of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the data set.
    !!
    !! @return The standard deviation of @p x.
    !!
    !! @remarks
    !! The standard deviation is computed as \f$ \sigma = 
    !! \sqrt{\frac{\sum_{i=1}^{n}(x_{i} - \mu)^{2}}{n - 1}} \f$
    function c_standard_deviation(n, x) bind(C, name = "c_standard_deviation") &
            result(s)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double) :: s

        ! Process
        s = standard_deviation(x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the range of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the data set.
    !!
    !! @return The range of the data set.
    function c_data_range(n, x) bind(C, name = "c_data_range") result(r)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double) :: r

        ! Process
        r = data_range(x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the z-score typically used in confidence interval
    !! calculations.
    !!
    !! @param[in] c The confidence level.  This value must be between 0
    !!  and 1 such that: 0 < c < 1.
    !! @param[out] z The computed z-score.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_INVALID_INPUT_ERROR: Occurs if @p c is not within its allowed
    !!      range.
    function c_z_score(c, z) bind(C, name = "c_z_score") result(flag)
        ! Arguments
        real(c_double), intent(in), value :: c
        real(c_double), intent(out) :: z
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        z = z_score(c, err)
        if (err%has_error_occurred()) flag = err%get_error_flag()
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the confidence interval of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the data set.
    !! @param[in] zval The critical value (z or t value).
    !!
    !! @return The confidence interval as referenced from the population
    !!  mean.
    !!
    !! @remarks
    !! This value is computed as follows.
    !! \f$ CI = z^{*} \frac{\sigma}{\sqrt{n}} \f$
    function c_confidence_interval(n, x, zval) &
            bind(C, name = "c_confidence_interval") result(ci)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n), zval
        real(c_double) :: ci

        ! Process
        ci = confidence_interval(x, zval)
    end function

! ------------------------------------------------------------------------------
    !> @brief Evaluates the probability distribution function of a normal
    !! distribution.
    !!
    !! @param[in] mu The population mean.
    !! @param[in] sigma The population standard deviation.
    !! @param[in] n The number of values at which to evaluat the function.
    !! @param[in] x An N-element array containing the values at which to 
    !!  evaluate the distrubition funciton.
    !! @param[out] f An N-element array where the function output will be
    !!  written.
    !!
    !! @remarks
    !! The normal distribution has the form 
    !! \f$ f(x) = \frac{1}{\sigma \sqrt{2 \pi}} 
    !! \exp(-\frac{1}{2}(\frac{x-\mu}{\sigma})^{2}) \f$.
    subroutine c_normal_distribution(mu, sigma, n, x, f) &
            bind(C, name = "c_normal_distribution")
        ! Arguments
        real(c_double), intent(in), value :: mu, sigma
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out) :: f(n)

        ! Process
        f = normal_distribution(mu, sigma, x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Evalautes the probability distribution function of Student's
    !! t-distribution.
    !!
    !! @param[in] dof The number of degrees of freedom of the data set.
    !! @param[in] n The number of values at which to evaluat the function.
    !! @param[in] t An N-element array containing the values at which to 
    !!  evaluate the distrubition funciton.
    !! @param[out] f An N-element array where the function output will be
    !!  written.
    !!
    !! @remarks
    !! Student's t-distribution has the form \f$ f(t) = 
    !! \frac{\Gamma(\frac{\nu + 1}{2})}{\sqrt{\nu \pi} 
    !! \Gamma(\frac{\nu}{2})} (1 + \frac{t^{2}}{\nu})^{-\frac{\nu + 1}{2}} 
    !! \f$.
    subroutine c_t_distribution(dof, n, t, f) bind(C, name = "c_t_distribution")
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in), value :: dof
        real(c_double), intent(in) :: t(n)
        real(c_double), intent(out) :: f(n)

        ! Process
        f = t_distribution(dof, t)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Evaluates the probability distribution function of a beta
    !! distribution.
    !!
    !! @param[in] a The first argument of the function.
    !! @param[in] b The second argument of the function.
    !! @param[in] n The number of values at which to evaluat the function.
    !! @param[in] x An N-element array containing the values at which to 
    !!  evaluate the distrubition funciton.
    !! @param[out] f An N-element array where the function output will be
    !!  written.
    !!
    !! @remarks The beta distribution has the form \f$ f(x) = 
    !! \frac{x^{a-1} (1 - x)^{b-1}}{\beta(a,b)} \f$.
    subroutine c_beta_distribution(a, b, n, x, f) &
            bind(C, name = "c_beta_distribution")
        ! Arguments
        real(c_double), intent(in), value :: a, b
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out) :: f(n)

        ! Process
        f = beta_distribution(a, b, x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Evaluates the probability distribution function of the
    !! F distribution.
    !!
    !! @param[in] d1 A model parameter.
    !! @param[in] d2 A model parameter.
    !! @param[in] n The number of values at which to evaluat the function.
    !! @param[in] x An N-element array containing the values at which to 
    !!  evaluate the distrubition funciton.
    !! @param[out] f An N-element array where the function output will be
    !!  written.
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
    subroutine c_f_distribution(d1, d2, n, x, f) &
            bind(C, name = "c_f_distribution")
        ! Arguments
        real(c_double), intent(in), value :: d1, d2
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out) :: f(n)

        ! Process
        f = f_distribution(d1, d2, x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Utilizes an analysis of variance (ANOVA) to determine the
    !! variance components of a measurement process represented by the
    !! supplied data set.
    !!
    !! @param[in] nparts The number of parts being analyzed (must be greater 
    !!  than 1).
    !! @param[in] ntests The number of tests performed on each part (must be
    !!  greater than 1).
    !! @param[in] nops The number of operators performing tests (must be at
    !!  least 1).
    !! @param[in] x An NPARTS-by-NTESTS-by-NOPS data set from the measurement 
    !!  process to analyze.
    !! @param[in] alpha A parameter used to determine the appropriate 
    !!  calculation path.  The recommended value is 0.05 (95% confidence level).
    !! @param[out] rst The resulting variance components.
    !!
    !! @remarks It is possible for this routine to return zero-valued
    !!  variance components.  In such an event it is recommended that
    !!  another variance estimator is utilized.
    subroutine c_anova(nparts, ntests, nops, x, alpha, rst) &
            bind(C, name = "c_anova")
        ! Arguments
        integer(c_int), intent(in), value :: nparts, ntests, nops
        real(c_double), intent(in) :: x(nparts,ntests, nops)
        real(c_double), intent(in), value :: alpha
        type(process_variance), intent(out) :: rst

        ! Process
        rst = anova(x, alpha)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Utilizes a control chart type approach to evaluate the 
    !! measurement process utilized to collect the supplied data set.
    !!
    !! @param[in] nparts The number of parts being analyzed (must be greater 
    !!  than 1).
    !! @param[in] ntests The number of tests performed on each part (must be
    !!  greater than 1).
    !! @param[in] nops The number of operators performing tests (must be at
    !!  least 1).
    !! @param[in] x An NPARTS-by-NTESTS-by-NOPS data set from the measurement 
    !!  process to analyze.
    !! @param[out] rst The resulting variance components.
    !!
    !! @remarks This approach is best suited for smaller data sets whose
    !!  dimension doesn't exceed 25-30.  Anything over this size is better
    !!  served by another technique.
    subroutine c_control_chart_variance(nparts, ntests, nops, x, rst) &
            bind(C, name = "c_control_chart_variance")
        ! Arguments
        integer(c_int), intent(in), value :: nparts, ntests, nops
        real(c_double), intent(in) :: x(nparts, ntests, nops)
        type(process_variance), intent(out) :: rst

        ! Process
        rst = control_chart_variance(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the gauge R&R statistics given a supplied set of
    !! process variance data.
    !!
    !! @param[in] k The multiplier to use when computing the P/T ratio.
    !!  Typically, a value of 6 is used for this factor.
    !! @param[in] x The process variance data.
    !! @param[in] usl The upper specification limit.
    !! @param[in] lsl The lower specification limit.
    !! @param[out] grr The gauge R&R statistics.
    subroutine c_compute_grr(k, x, usl, lsl, grr) &
            bind(C, name = "c_compute_grr")
        ! Arguments
        real(c_double), intent(in), value :: k, usl, lsl
        type(process_variance), intent(in) :: x
        type(grr_results), intent(out) :: grr

        ! Process
        grr = compute_grr(k, x, usl, lsl)
    end subroutine

! ------------------------------------------------------------------------------
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
    function c_discrimination_ratio(tv, mv) &
            bind(C, name = "c_discrimination_ratio") result(x)
        ! Arguments
        real(c_double), intent(in), value :: tv, mv
        real(c_double) :: x

        ! Process
        x = discrimination_ratio(tv, mv)
    end function
    
! ******************************************************************************
! SPECIAL FUNCTIONS
! ------------------------------------------------------------------------------
    !> @brief Computes the value of the regularized beta function.
    !!
    !! @param[in] x The upper limit of the integration.
    !! @param[in] a The first argument of the function.
    !! @param[in] b The second argument of the function.
    !!
    !! @return The value of the regularized beta function at @p a and @p b.
    !!
    !! @remarks The regularized beta function is defined as \f$ 
    !! I_{x}(a, b) = \frac{\beta(x; a, b)}{\beta(a, b)} \f$.
    function c_regularized_beta(x, a, b) &
            bind(C, name = "c_regularized_beta") result(z)
        ! Arguments
        real(c_double), intent(in), value :: x, a, b
        real(c_double) :: z

        ! Process
        z = regularized_beta(x, a, b)
    end function

! ------------------------------------------------------------------------------
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
    function c_beta(a, b) bind(C, name = "c_beta") result(z)
        ! Arguments
        real(c_double), intent(in), value :: a, b
        real(c_double) :: z

        ! Process
        z = beta(a, b)
    end function
! ------------------------------------------------------------------------------
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
    function c_incomplete_beta(x, a, b) bind(C, name = "incomplete_beta") &
            result(z)
        ! Arguments
        real(c_double), intent(in), value :: x, a, b
        real(c_double) :: z

        ! Process
        z = incomplete_beta(x, a, b)
    end function

! ******************************************************************************
! INTERPOLATION
! ------------------------------------------------------------------------------
    !> @brief Performs a polynomial interpolation.
    !!
    !! @param[in] order The order of the polynomial.  This value must be at 
    !!  least 1, but not exceed @p npts - 1.
    !! @param[in] npts The number of data points.
    !! @param[in] x An @p npts element array containing the independent variable
    !!  data.  This array must be monotonic.
    !! @param[in] y An @p npts element array containing the dependent variable
    !!  data.
    !! @param[in] ni The number of points to interpolate.
    !! @param[in] xi An @p ni element array containing the independent variable
    !!  values at which to interpolate.
    !! @param[out] yi An @p ni element array corresponding to @p xi where the
    !!  interpolated values will be written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - M_INVALID_INPUT_ERROR: Occurs if @p order is less than 1.
    !!  - M_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
    !!      increasing or decreasing.
    function c_interpolate(order, npts, x, y, ni, xi, yi) &
            bind(C, name = "c_interpolate") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: order, npts, ni
        real(c_double), intent(in) :: x(npts), y(npts), xi(ni)
        real(c_double), intent(out) :: yi(ni)
        integer(c_int) :: flag

        ! Local Variables
        type(polynomial_interp) :: interp
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        call interp%initialize(x, y, order, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
        yi = interp%interpolate(xi)
    end function

! ------------------------------------------------------------------------------
    !> @brief Performs a spline interpolation.
    !!
    !! @param[in] npts The number of data points.
    !! @param[in] x An @p npts element array containing the independent variable
    !!  data.  This array must be monotonic.
    !! @param[in] y An @p npts element array containing the dependent variable
    !!  data.
    !! @param[in] ni The number of points to interpolate.
    !! @param[in] xi An @p ni element array containing the independent variable
    !!  values at which to interpolate.
    !! @param[out] yi An @p ni element array corresponding to @p xi where the
    !!  interpolated values will be written.
    !! @param[in] ibcbeg An input that defines the nature of the
    !!  boundary condition at the beginning of the spline.
    !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
    !!      initial interval.  No value is required for @p ybcbeg.  This is
    !!      often considered a natural boundary condition.
    !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at 
    !!      its initial point is provided in @p ybcbeg.
    !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at 
    !!      its initial point is provided in @p ybcbeg.
    !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is 
    !!      continuous at x(2).  No value is required for @p ybcbeg.
    !! @param[in] ybcbeg If needed, the value of the initial point boundary
    !!  condition.
    !! @param[in] ibcend An input that defines the nature of the
    !!  boundary condition at the end of the spline.
    !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
    !!      final interval.  No value is required for @p ybcend.  This is
    !!      often considered a natural boundary condition.
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
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - M_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
    !!      increasing or decreasing.
    function c_spline(npts, x, y, ni, xi, yi, ibcbeg, ybcbeg, ibcend, ybcend) &
            bind(C, name = "c_spline") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: npts, ni, ibcbeg, ibcend
        real(c_double), intent(in), value :: ybcbeg, ybcend
        real(c_double), intent(in) :: x(npts), y(npts), xi(ni)
        real(c_double), intent(out) :: yi(ni)
        integer(c_int) :: flag

        ! Local Variables
        type(spline_interp) :: interp
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        call interp%initialize_spline(x, y, ibcbeg, ybcbeg, ibcend, ybcend, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
        yi = interp%interpolate(xi)
    end function

! ------------------------------------------------------------------------------
    !> @brief Smooths a data set using a robust locally weighted scatterplot 
    !! smoothing (LOWESS) algorithm. 
    !!
    !! @param[in] npts The number of data points.
    !! @param[in] x An @p npts element array containing the independent variable
    !!  data.  This array must be monotonic.
    !! @param[in] y An @p npts element array containing the dependent variable
    !!  data.
    !! @param[in] factor Specifies the amount of smoothing.  More specifically, 
    !!  this value is the fraction of points used to compute each value.  
    !!  As this value increases, the output becomes smoother.  Choosing a 
    !!  value in the range of 0.2 to 0.8 usually results in a good fit.  As 
    !!  such, a reasonable starting point, in the absence of better 
    !!  information, is a value of 0.5.
    !! @param[out] ys An @p npts element array where the smoothed data will
    !!  be written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - M_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
    !!      increasing or decreasing.
    function c_lowess_smoothing(npts, x, y, factor, ys) &
            bind(C, name = "c_lowess_smoothing") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: npts
        real(c_double), intent(in), value :: factor
        real(c_double), intent(in) :: x(npts), y(npts)
        real(c_double), intent(out) :: ys(npts)
        integer(c_int) :: flag

        ! Local Variables
        type(lowess_smoothing) :: fit
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        call fit%initialize(x, y, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if

        ys = fit%smooth(factor, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
