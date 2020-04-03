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
        integer(c_int), intent(in), value :: dof, n
        real(c_double), intent(in) :: t(n)
        real(c_double), intent(out) :: f(n)

        ! Process
        f = t_distribution(dof, t)
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
