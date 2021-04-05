! measurements_distribtion.f90

submodule (measurements_core) measurements_distribtion
contains
! ------------------------------------------------------------------------------
pure elemental module function normal_distribution_pdf(mu, sigma, x) result(f)
    ! Arguments
    real(real64), intent(in) :: mu, sigma, x
    real(real64) :: f

    ! Constants
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Process
    f = (1.0d0 / (sigma * sqrt(2.0d0 * pi))) * &
        exp(-0.5 * ((x - mu) / sigma)**2)
end function

! ------------------------------------------------------------------------------
pure elemental module function normal_distribution_cdf(mu, sigma, x) result(rst)
    ! Arguments
    real(real64), intent(in) :: mu, sigma, x
    real(real64) :: rst

    ! Process
    rst = 0.5d0 * (1.0d0 + erf((x - mu) / (sqrt(2.0d0) * sigma)))
end function

! ------------------------------------------------------------------------------
pure elemental module function log_normal_distribution_pdf(mu, sigma, x) &
        result(rst)
    ! Arguments
    real(real64), intent(in) :: mu, sigma, x
    real(real64) :: rst

    ! Constants
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Process
    rst = (1.0d0 / (x * sigma * sqrt(2.0d0 * pi))) * &
        exp(-(log(x) - mu)**2 / (2.0d0 * sigma**2))
end function

! ------------------------------------------------------------------------------
pure elemental module function log_normal_distribution_cdf(mu, sigma, x) &
        result(rst)
    ! Arguments
    real(real64), intent(in) :: mu, sigma, x
    real(real64) :: rst

    ! Process
    rst = 0.5d0 + 0.5d0 * erf((log(x) - mu) / (sigma * sqrt(2.0d0)))
end function

! ------------------------------------------------------------------------------
pure elemental module function t_distribution_pdf(dof, t) result(f)
    ! Arguments
    real(real64), intent(in) :: dof
    real(real64), intent(in) :: t
    real(real64) :: f

    ! Constants
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    real(real64) :: arg

    ! Process
    arg = 0.5d0 * (dof + 1)
    f = (gamma(arg) / (gamma(0.5d0 * dof) * sqrt(dof * pi))) * &
        (1.0d0 + t**2 / dof)**(-arg)
end function

! ------------------------------------------------------------------------------
pure elemental module function t_distribution_cdf(dof, x) result(rst)
    ! Arguments
    real(real64), intent(in) :: dof, x
    real(real64) :: rst

    ! Local Variables
    real(real64) :: z

    ! Process
    z = dof / (x**2 + dof)
    if (x < 0.0d0) then
        rst = 0.5d0 * regularized_beta(z, 0.5d0 * dof, 0.5d0)
    else
        rst = 1.0d0 - 0.5d0 * regularized_beta(z, 0.5d0 * dof, 0.5d0)
    end if
end function

! ------------------------------------------------------------------------------
pure elemental module function beta_distribution_pdf(a, b, x) result(z)
    ! Arguments
    real(real64), intent(in) :: a, b, x
    real(real64) :: z

    ! Process
    z = x**(a - 1.0d0) * (1.0d0 - x)**(b - 1.0d0) / beta(a, b)
end function

! ------------------------------------------------------------------------------
pure elemental module function beta_distribution_cdf(a, b, x) result(rst)
    ! Arguments
    real(real64), intent(in) :: a, b, x
    real(real64) :: rst

    ! Process
    rst = regularized_beta(x, a, b)
end function

! ------------------------------------------------------------------------------
pure elemental module function f_distribution_pdf(d1, d2, x) result(z)
    ! Arguments
    real(real64), intent(in) :: d1, d2, x
    real(real64) :: z

    ! Local Variables
    real(real64) :: arg

    ! Process
    arg = ((d1 * x)**d1) * (d2**d2) / ((d1 * x + d2)**(d1 + d2))
    z = sqrt(arg) / (x * beta(0.5d0 * d1, 0.5d0 * d2))
end function

! ------------------------------------------------------------------------------
pure elemental module function f_distribution_cdf(d1, d2, x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x, d1, d2
    real(real64) :: rst

    ! Local Variables
    real(real64) :: arg

    ! Process
    arg = d1 * x / (d1 * x + d2)
    rst = regularized_beta(arg, 0.5d0 * d1, 0.5d0 * d2)
end function

! ******************************************************************************
! NORMAL_DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental module function nrm_pdf(this, x) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = normal_distribution_pdf( &
        this%get_mean(), &
        this%get_standard_deviation(), &
        x &
    )
end function

! ------------------------------------------------------------------------------
pure elemental module function nrm_cdf(this, x) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = normal_distribution_cdf( &
        this%get_mean(), &
        this%get_standard_deviation(), &
        x &
    )
end function

! ------------------------------------------------------------------------------
pure module function nrm_mean(this) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%m_mean
end function

! ------------------------------------------------------------------------------
pure module function nrm_median(this) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%m_mean
end function

! ------------------------------------------------------------------------------
pure module function nrm_mode(this) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%m_mean
end function

! ------------------------------------------------------------------------------
pure module function nrm_variance(this) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%m_sigma**2
end function

! ------------------------------------------------------------------------------
module subroutine nrm_set_params(this, x)
    class(normal_distribution), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x
    this%m_mean = x(1)
    this%m_sigma = x(2)
end subroutine

! ------------------------------------------------------------------------------
pure module function nrm_get_mean(this) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%m_mean
end function

! --------------------
module subroutine nrm_set_mean(this, x)
    class(normal_distribution), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_mean = x
end subroutine

! ------------------------------------------------------------------------------
pure module function nrm_get_sigma(this) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%m_sigma
end function

! --------------------
module subroutine nrm_set_sigma(this, x)
    class(normal_distribution), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_sigma = x
end subroutine

! ******************************************************************************
! STUDENT'S T DISTRIBUTION
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule
