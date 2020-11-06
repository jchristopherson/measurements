! measurements_c_api.f90

!> @brief A module providing C bindings to the MEASUREMENTS library.
module measurements_c_api
    use iso_c_binding
    use iso_fortran_env
    use measurements_core
    use ferror
    implicit none

    interface
        !> @brief Defines a window function.
        !!
        !! @param[in] bin The index or bin number (0 <= @p bin <= @p n).
        !! @param[in] n The transform length.
        !! @return The window function value.
        function c_window_function(bin, n) result(x)
            use iso_c_binding
            integer(c_int), intent(in), value :: bin, n
            real(c_double) :: x
        end function
    end interface

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
    !> @brief Computes the t-score typically used in confidence interval
    !! calculations when the population size is limited.
    !!
    !! @param[in] c The confidence level.  This value must be between 0
    !!  and 1 such that: 0 < c < 1.
    !! @param[in] n The population size.
    !! @param[out] t The computed t-score.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_INVALID_INPUT_ERROR: Occurs if @p c is not within its allowed
    !!      range.
    function c_t_score(c, n, t) bind(C, name = "c_t_score") result(flag)
        ! Arguments
        real(c_double), intent(in), value :: c
        integer(c_int), intent(in), value :: n
        real(c_double), intent(out) :: t
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        t = t_score(c, n, err)
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
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in), value :: zval
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
    !> @brief Applies Student's t-test to compute the t-statistic.  A 
    !! two-tailed distribution is assumed.
    !!
    !! @param[in] n1 The number of data points in the first data set.
    !! @param[in] x1 An @p n1 element array containing the first data set.
    !! @param[in] n2 The number of data points in the second data set.
    !! @param[in] x2 An @p n2 element array containing the second data set.
    !! @param[in] method An input flag defining which method to utilize.
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
    !! @param[out] rst Student's t-statistic and the associated probability term
    !! that establishes the significance of the t-statistic.
    subroutine c_t_test(n1, x1, n2, x2, method, rst) bind(C, name = "c_t_test")
        ! Arguments
        integer(c_int), intent(in), value :: n1, n2, method
        real(c_double), intent(in) :: x1(n1), x2(n2)
        type(statistic), intent(out) :: rst

        ! Process
        rst = t_test(x1, x2, method)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Applies the F-test for different variances.
    !!
    !! @param[in] n1 The number of data points in the first data set.
    !! @param[in] x1 An @p n1 element array containing the first data set.
    !! @param[in] n2 The number of data points in the second data set.
    !! @param[in] x2 An @p n2 element array containing the second data set.
    !! @param[out] rst The F-test statistic and associated probability term that
    !! establishes the signficance of the f-statistic.
    subroutine c_f_test(n1, x1, n2, x2, rst) bind(C, name = "c_f_test")
        ! Arguments
        integer(c_int), intent(in), value :: n1, n2
        real(c_double), intent(in) :: x1(n1), x2(n2)
        type(statistic), intent(out) :: rst

        ! Process
        rst = f_test(x1, x2)
    end subroutine
    
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
    !> @brief Applies a moving average to smooth a data set.
    !!
    !! @param[in] npts The number of data points.
    !! @param[in,out] x An @p npts element array that on input contains the 
    !!  signal to smooth.  On output, the smoothed signal.
    !! @param[in] navg The size of the averaging window.  This value must be
    !!  at least 2, but no more than the number of elements in @p x.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_INVALID_INPUT_ERROR: Occurs if @p navg is less than 2, or 
    !!      greater than @p npts.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    function c_moving_average(npts, x, navg) &
            bind(C, name = "c_moving_average") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: npts, navg
        real(c_double), intent(inout) :: x(npts)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        call moving_average(x, navg, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ******************************************************************************
! REGRESSION
! ------------------------------------------------------------------------------
    !> @brief Fits the multiple input, multiple output linear model
    !! A * X = Y by solving for matrix A in a least-squares sense.
    !!
    !! @param[in] m The number of rows in matrix A.
    !! @param[in] n The number of columns in matrix A.
    !! @param[in] k The number of data points to fit (number of columns in
    !!  either X or Y).
    !! @param[in] x An N-by-K matrix of known independent variables.  K
    !!  must be greater than or equal to N.
    !! @param[in] ldx The leading dimension of matrix X.
    !! @param[in] y An M-by-K matrix of known dependent variables.  Notice,
    !!  M must be less than or equal to N, and K must be greater than or
    !!  equal to M.
    !! @param[in] ldy The leading dimension of matrix Y.
    !! @param[out] The M-by-N coefficient matrix A.
    !! @param[in] lda The leading dimension of matrix A.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - M_ARRAY_SIZE_ERROR: Occurs if there is a size-mismatch in the
    !!      matrix equation.
    !!  - M_UNDERDEFINED_PROBLEM: Occurs if there is insufficient data 
    !!      (e.g. k < n).
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
    function c_linear_least_squares_mimo(m, n, k, x, ldx, y, ldy, a, lda) &
            bind(C, name = "c_linear_least_squares_mimo") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, k, ldx, ldy, lda
        real(c_double), intent(in) :: x(ldx,*), y(ldy,*)
        real(c_double), intent(out) :: a(lda,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Input Check
        if (ldx < n) then
            flag = M_ARRAY_SIZE_ERROR
            return
        end if

        if (ldy < m) then
            flag = M_ARRAY_SIZE_ERROR
            return
        end if

        if (lda < m) then
            flag = M_ARRAY_SIZE_ERROR
            return
        end if

        ! Process
        a(1:m,1:n) = linear_least_squares_mimo(x(1:n,1:k), y(1:m,1:k), &
            err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Fits the multiple input, single output linear model
    !! A * X = Y by solving for array A in a least-squares sense.
    !!
    !! @param[in] n The number of coefficients to find.
    !! @param[in] k The number of data points to fit.
    !! @param[in] x An N-by-K matrix of known independent variables.  K
    !!  must be greater than or equal to N.
    !! @param[in] ldx The leading dimension of matrix X.
    !! @param[in] y A K-element array containing the known dependent 
    !!  variables.
    !! @param[out] The N element coefficient array A.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - M_ARRAY_SIZE_ERROR: Occurs if there is a size-mismatch in the
    !!      matrix equation.
    !!  - M_UNDERDEFINED_PROBLEM: Occurs if there is insufficient data 
    !!      (e.g. k < n).
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
    function c_linear_least_squares_miso(n, k, x, ldx, y, a) &
            bind(C, name = "c_linear_least_squares_miso") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, k, ldx
        real(c_double), intent(in) :: x(ldx,*), y(k)
        real(c_double), intent(out) :: a(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Input Check
        if (ldx < n) then
            flag = M_ARRAY_SIZE_ERROR
            return
        end if

        ! Process
        a = linear_least_squares_miso(x(1:n,1:k), y, err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ******************************************************************************
! PEAK DETECTION
! ------------------------------------------------------------------------------
    !> @brief Attempts to locate all peaks and valleys within the given
    !! data set bound by the constraints specified.
    !!
    !! @param[in] n The number of data points in the array to search.
    !! @param[in] x The N-element array to search.
    !! @param[in] delta A threshold level used to denote the minimum change
    !!  acceptable in determining a peak or valley from neighboring points.
    !! @param[in] szmxind The size of the @p mxind buffer.
    !! @param[out] mxind A @p szmxind element array where the indices of the
    !!  peak values will be written.  The indices are zero-based.
    !! @param[out] nmxind The actual number of peak value indices written to
    !!  @p mxind.
    !! @param[in] szmnind The size of the @p mnind buffer.
    !! @param[out] mnind A @p szmnind element array where the indices of the
    !!  valley values will be written.  The indices are zero-based.
    !! @param[out] nmnind The actual number of valley value indices written to
    !!  @p mnind.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    function c_peak_detect(n, x, delta, szmxind, mxind, nmxind, szmnind, &
            mnind, nmnind) bind(C, name = "c_peak_detect") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, szmxind, szmnind
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in), value :: delta
        integer(c_int), intent(out) :: mxind(szmxind), mnind(szmnind), nmxind, &
            nmnind
        integer(c_int) :: flag

        ! Local Variables
        integer(c_int) :: npts
        type(errors) :: err
        type(peak_info) :: pks

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        pks = peak_detect(x, delta, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if

        ! Copy over the results - subtract 1 to zero-base the indices
        npts = min(szmxind, size(pks%max_value_indices))
        mxind(1:npts) = pks%max_value_indices(1:npts) - 1
        nmxind = npts

        npts = min(szmnind, size(pks%min_value_indices))
        mnind(1:npts)  = pks%min_value_indices(1:npts) - 1
        nmnind = npts
    end function

! ------------------------------------------------------------------------------
    !> @brief Tests to see if a number is an integer power of two.
    !!
    !! @param[in] n The integer to test.
    !! @return Returns true if @p n is a power of two; else, false.
    function c_is_power_of_two(n) bind(C, name = "c_is_power_of_two") &
            result(rst)
        integer(c_int), intent(in), value :: n
        logical(c_bool) :: rst
        rst = is_power_of_two(n)
    end function

! ------------------------------------------------------------------------------
    !> @brief Provides the next higher integer power of two.
    !!
    !! @param[in] x The value to test.
    !! @return The next power of two higher than @p x.  If @p x is already
    !! a power of two, it's value is simply returned.  For instance, if @p
    !! is set to 128, then a value of 7 is returned.  However, if a value
    !! of 129 is supplied, then a value of 8 is returned.
    function c_next_power_of_two(x) bind(C, name = "c_next_power_of_two") &
            result(n)
        integer(c_int), intent(in), value :: x
        integer(c_int) :: n
        n = next_power_of_two(x);
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the Fourier transform of a discretely sampled signal.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x The N-element array containing the data to transform.
    !! @param[in] nf The number of elements in @p f.  Ideally, this should be
    !!  N / 2 + 1 if N is even; else, (N + 1) / 2 if N is odd.
    !! @param[out] f The positive half and DC component of the transformed data.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    function c_fourier_transform(n, x, nf, f) &
            bind(C, name = "c_fourier_transform") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, nf
        real(c_double), intent(in) :: x(n)
        complex(c_double), intent(out) :: f(nf)
        integer(c_int) :: flag

        ! Local Variables
        integer(int32) :: m
        complex(real64), allocatable, dimension(:) :: fc
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)
        if (mod(n, 2) == 0) then
            m = min(nf, n / 2 + 1)
        else
            m = min(nf, (n + 1) / 2)
        end if

        ! Process
        fc = fourier_transform(x, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
        f(1:m) = fc(1:m)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the periodogram power spectral density (PSD) estimate
    !! of a signal.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x The N-element array containing the data to transform.
    !! @param[in] winfun The window function to apply.
    !! @param[in] nfft The length Fourier transform to apply.  This must
    !!  be an integer power of two, even if @p x is not an integer power
    !!  of two in length.  If this parameter is less than the length of
    !!  @p x, @p x will be overlapped, windowed, and averaged to achieve
    !!  an estimate of the power spectrum.  If this parameter is larger
    !!  than the length of @p x, @p x will be padded with zeros prior
    !!  to windowing.
    !! @param[in] np
    !! @param[out] p The NP-element array containing the periodogram.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - M_INVALID_INPUT_ERROR: Occurs if @p nfft is not an integer power
    !!      of two.
    function c_periodogram(n, x, winfun, nfft, np, p) &
            bind(C, name = "c_periodogram") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, nfft, np
        real(c_double), intent(in) :: x(n)
        type(c_funptr), intent(in), value :: winfun
        real(c_double), intent(out) :: p(np)
        integer(c_int) :: flag

        ! Local Variables
        procedure(c_window_function), pointer :: cfptr
        procedure(window_function), pointer :: fptr
        type(errors) :: err
        real(real64), allocatable, dimension(:) :: pc
        integer(int32) :: m

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)
        call c_f_procpointer(winfun, cfptr)
        fptr => win
        m = min(np, nfft / 2 + 1)

        ! Process
        pc = periodogram(x, fptr, nfft, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
        p(1:m) = pc(1:m)
    contains
        function win(jj, nn) result(xx)
            integer(int32), intent(in) :: jj, nn
            real(real64) :: xx
            xx = cfptr(jj, nn)
        end function
    end function

! ------------------------------------------------------------------------------
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
    function c_fourier_frequency(fs, i, nxfrm) &
            bind(C, name = "c_fourier_frequency") result(f)
        real(c_double), intent(in), value :: fs
        integer(c_int), intent(in), value :: i, nxfrm
        real(c_double) :: f
        f = fourier_frequency(fs, i, nxfrm)
    end function

! ------------------------------------------------------------------------------
    !> @brief Defines a rectangular window.
    !!
    !! @param[in] j The index or bin number (0 <= @p bin <= @p n).
    !! @param[in] n The transform length.
    !!
    !! @return The value of the window function at index @p j.
    function c_rectangular_window(j, n) bind(C, name = "c_rectangular_window") &
            result(x)
        integer(c_int), intent(in), value :: j, n
        real(c_double) :: x
        x = rectangular_window(j, n)
    end function

! ------------------------------------------------------------------------------
    !> @brief Defines a Hann window.
    !!
    !! @param[in] j The index or bin number (0 <= @p bin <= @p n).
    !! @param[in] n The transform length.
    !!
    !! @return The value of the window function at index @p j.
    function c_hann_window(j, n) bind(C, name = "c_hann_window") result(x)
        integer(c_int), intent(in), value :: j, n
        real(c_double) :: x
        x = hann_window(j, n)
    end function

! ------------------------------------------------------------------------------
    !> @brief Defines a Hamming window.
    !!
    !! @param[in] j The index or bin number (0 <= @p bin <= @p n).
    !! @param[in] n The transform length.
    !!
    !! @return The value of the window function at index @p j.
    function c_hamming_window(j, n) bind(C, name = "c_hamming_window") result(x)
        integer(c_int), intent(in), value :: j, n
        real(c_double) :: x
        x = hamming_window(j, n)
    end function

! ------------------------------------------------------------------------------
    !> @brief Defines a Welch window.
    !!
    !! @param[in] j The index or bin number (0 <= @p bin <= @p n).
    !! @param[in] n The transform length.
    !!
    !! @return The value of the window function at index @p j.
    function c_welch_window(j, n) bind(C, name = "c_welch_window") result(x)
        integer(c_int), intent(in), value :: j, n
        real(c_double) :: x
        x = welch_window(j, n)
    end function

! ------------------------------------------------------------------------------
    !> @brief Defines a Blackman-Harris window.
    !!
    !! @param[in] j The index or bin number (0 <= @p bin <= @p n).
    !! @param[in] n The transform length.
    !!
    !! @return The value of the window function at index @p j.
    function c_blackman_harris_window(j, n) &
            bind(C, name = "c_blackman_harris_window") result(x)
        integer(c_int), intent(in), value :: j, n
        real(c_double) :: x
        x = blackman_harris_window(j, n)
    end function

! ------------------------------------------------------------------------------
    !> @brief Applies a low-pass filter to a signal.
    !!
    !! @param[in] n The length of the input array.
    !! @param[in] x An N-element array containing the signal to filter.
    !! @param[in] fs The frequency at which @p x was sampled.
    !! @param[in] cutoff The cut-off frequency.  This value must be
    !!  positive-valued, and must be less than the Nyquist frequency.
    !! @param[out] y An N-element array where the filtered signal will be
    !!  written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff is either not positive
    !!      valued, or is greater than or equal to the Nyquist frequency.
    function c_low_pass_filter(n, x, fs, cutoff, y) &
            bind(C, name = "c_low_pass_filter") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in), value :: fs, cutoff
        real(c_double), intent(out) :: y(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        y = low_pass_filter(x, fs, cutoff, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Applies a high-pass filter to a signal.
    !!
    !! @param[in] n The length of the input array.
    !! @param[in] x An N-element array containing the signal to filter.
    !! @param[in] fs The frequency at which @p x was sampled.
    !! @param[in] cutoff The cut-off frequency.  This value must be
    !!  positive-valued, and must be less than the Nyquist frequency.
    !! @param[out] y An N-element array where the filtered signal will be
    !!  written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff is either not positive
    !!      valued, or is greater than or equal to the Nyquist frequency.
    function c_high_pass_filter(n, x, fs, cutoff, y) &
            bind(C, name = "c_high_pass_filter") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in), value :: fs, cutoff
        real(c_double), intent(out) :: y(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        y = high_pass_filter(x, fs, cutoff, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Applies a band-pass filter to a signal.
    !!
    !! @param[in] n The length of the input array.
    !! @param[in] x An N-element array containing the signal to filter.
    !! @param[in] fs The frequency at which @p x was sampled.
    !! @param[in] cutoff1 The lower cut-off frequency.  This value must be
    !!  positive-valued, and must be less than the Nyquist frequency.
    !! @param[in] cutoff2 The upper cut-off frequency.  This value must be
    !!  positive-valued, and must be less than the Nyquist frequency.
    !! @param[out] y An N-element array where the filtered signal will be
    !!  written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff1 or @p cutoff2 is 
    !!      either not positive-valued, or is greater than or equal to the 
    !!      Nyquist frequency.
    function c_band_pass_filter(n, x, fs, cutoff1, cutoff2, y) &
            bind(C, name = "c_band_pass_filter") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in), value :: fs, cutoff1, cutoff2
        real(c_double), intent(out) :: y(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        y = band_pass_filter(x, fs, cutoff1, cutoff2, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Applies a band-stop filter to a signal.
    !!
    !! @param[in] n The length of the input array.
    !! @param[in] x An N-element array containing the signal to filter.
    !! @param[in] fs The frequency at which @p x was sampled.
    !! @param[in] cutoff1 The lower cut-off frequency.  This value must be
    !!  positive-valued, and must be less than the Nyquist frequency.
    !! @param[in] cutoff2 The upper cut-off frequency.  This value must be
    !!  positive-valued, and must be less than the Nyquist frequency.
    !! @param[out] y An N-element array where the filtered signal will be
    !!  written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - M_INVALID_INPUT_ERROR: Occurs if @p cutoff1 or @p cutoff2 is 
    !!      either not positive-valued, or is greater than or equal to the 
    !!      Nyquist frequency.
    function c_band_stop_filter(n, x, fs, cutoff1, cutoff2, y) &
            bind(C, name = "c_band_stop_filter") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in), value :: fs, cutoff1, cutoff2
        real(c_double), intent(out) :: y(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        y = band_stop_filter(x, fs, cutoff1, cutoff2, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes an FFT of a data set.  The results of the transform
    !! are normalized such that an inverse transform will result in the 
    !! original signal.
    !!
    !! @param[in] n The length of the array.
    !! @param[in] x An N-element array containing the data to transform.
    !! @param[out] f An N-element array where the transformed data will be 
    !!  written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    function c_fft(n, x, f) bind(C, name = "c_fft") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        complex(c_double), intent(out) :: f(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        f = fft(x, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the inverse FFT of a data set.
    !!
    !! @param[in] n The length of the array.
    !! @param[in] x An N-element array containing the data to transform.
    !! @param[out] f An N-element array where the transformed data will be 
    !!  written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    function c_ifft(n, x, f) bind(C, name = "c_ifft") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        complex(c_double), intent(in) :: x(n)
        complex(c_double), intent(out) :: f(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        f = ifft(x, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the cumulative sum of an array.
    !!
    !! @param[in] n The length of the array.
    !! @param[in] x The N-element array on which to operate.
    !! @param[out] y The N-element array where the solution will be written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    function c_cumulative_sum(n, x, y) bind(C, name = "c_cumulative_sum") &
            result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out) :: y(n)
        integer(c_int) :: flag

        ! Local Variabes
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        y = cumulative_sum(x, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the differences between each element in the array.
    !!
    !! @param[in] n The length of the array.
    !! @param[in] x The N-element array on which to operate.
    !! @param[out] dx An N-1 element array where the solution will be written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    function c_difference(n, x, dx) bind(C, name = "c_difference") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out) :: dx(n-1)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        dx = difference(x, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Unwraps an array of phase angles.
    !!
    !! @param[in] n The array length.
    !! @param[in] x The N-element array to unwrap.
    !! @param[in] cutoff An input that specifies the threshold value to use 
    !!  when unwrapping steps in @p x.
    !! @param[out] y The N-element array where the results will be written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    function c_unwrap(n, x, cutoff, y) bind(C, name = "c_unwrap") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in), value :: cutoff
        real(c_double), intent(out) :: y(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        y = unwrap(x, cutoff, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Estimates the derivative of a discretely sampled signal by means
    !! of finite differences.
    !!
    !! @param[in] n The length of the arrays.
    !! @param[in] x The N-element array containing the independent variable 
    !!  data.
    !! @parma[in] y The N-element array containing the dependent variable
    !!  data.
    !! @param[out] dydx An N-element array where the results will be written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    function c_finite_difference(n, x, y, dydx) &
            bind(C, name = "c_finite_difference") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n), y(n)
        real(c_double), intent(out) :: dydx(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        dydx = finite_difference(x, y, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Estimates the indefinite integral of a discretely sampled signal.
    !!
    !! @param[in] n The length of the arrays.
    !! @param[in] x The N-element array containing the independent variable 
    !!  data.
    !! @parma[in] y The N-element array containing the dependent variable
    !!  data.
    !! @param[in] c The initial condition.
    !! @param[out] f An N-element array where the results will be written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    function c_integrate(n, x, y, c, f) bind(C, name = "c_integrate") &
            result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n), y(n)
        real(c_double), intent(in), value :: c
        real(c_double), intent(out) :: f(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        f = integrate(x, y, c, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Estimates the definite integral of a discretely sampled signal
    !! using a trapezoidal approach.
    !!
    !! @param[in] n The length of the arrays.
    !! @param[in] x The N-element array containing the independent variable 
    !!  data.
    !! @parma[in] y The N-element array containing the dependent variable
    !!  data.
    !!
    !! @return The result of the integration.
    function c_trapz_integrate(n, x, y) bind(C, name = "c_trapz_integrate") &
            result(f)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n), y(n)
        real(c_double) :: f

        ! Process
        f = trapz_integrate(x, y)
    end function

! ------------------------------------------------------------------------------
    !> @brief Removes the DC offset from a signal.
    !!
    !! @param[in] n The array length.
    !! @param[in] x The N-element array on which to operate.
    !! @param[out] y An N-element array where the results will be written.
    !!
    !! @return An error flag with the following possible values.
    !!  - M_NO_ERROR: No error occurred.  Normal operation.
    !!  - M_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    function c_remove_dc_offset(n, x, y) bind(C, name = "c_remove_dc_offset") &
            result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out) :: y(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Initialization
        flag = M_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Process
        y = remove_dc_offset(x, err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
end module
