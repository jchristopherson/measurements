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
    end interface

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
