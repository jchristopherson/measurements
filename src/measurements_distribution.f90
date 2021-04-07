! measurements_distribtion.f90

submodule (measurements_core) measurements_distribtion
    use :: ieee_arithmetic
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
        this%get_mu(), &
        this%get_sigma(), &
        x &
    )
end function

! ------------------------------------------------------------------------------
pure elemental module function nrm_cdf(this, x) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = normal_distribution_cdf( &
        this%get_mu(), &
        this%get_sigma(), &
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
pure module function nrm_get_mu(this) result(rst)
    class(normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%m_mean
end function

! --------------------
module subroutine nrm_set_mu(this, x)
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
! LOG NORMAL DISTRIBUTION
! ------------------------------------------------------------------------------
pure elemental module function lnrm_pdf(this, x) result(rst)
    class(log_normal_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = log_normal_distribution_pdf( &
        this%get_mu(), &
        this%get_sigma(), &
        x &
    )
end function

! ------------------------------------------------------------------------------
pure elemental module function lnrm_cdf(this, x) result(rst)
    class(log_normal_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = log_normal_distribution_cdf( &
        this%get_mu(), &
        this%get_sigma(), &
        x &
    )
end function

! ------------------------------------------------------------------------------
pure module function lnrm_mean(this) result(rst)
    class(log_normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = exp(this%m_mean + 0.5d0 * this%m_sigma**2)
end function

! ------------------------------------------------------------------------------
pure module function lnrm_median(this) result(rst)
    class(log_normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = exp(this%m_mean)
end function

! ------------------------------------------------------------------------------
pure module function lnrm_mode(this) result(rst)
    class(log_normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = exp(this%m_mean - this%m_sigma**2)
end function

! ------------------------------------------------------------------------------
pure module function lnrm_variance(this) result(rst)
    class(log_normal_distribution), intent(in) :: this
    real(real64) :: rst
    rst = (exp(this%m_sigma**2) - 1.0d0) * &
        exp(2.0d0 * this%m_mean + this%m_sigma**2)
end function

! ******************************************************************************
! STUDENT'S T DISTRIBUTION
! ------------------------------------------------------------------------------
pure module function td_get_dof(this) result(rst)
    class(t_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%m_dof
end function

! --------------------
module subroutine td_set_dof(this, x)
    class(t_distribution), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_dof = x
end subroutine

! ------------------------------------------------------------------------------
pure elemental module function td_pdf(this, x) result(rst)
    class(t_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = t_distribution_pdf(this%get_dof(), x)
end function

! ------------------------------------------------------------------------------
pure elemental module function td_cdf(this, x) result(rst)
    class(t_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = t_distribution_cdf(this%get_dof(), x)
end function

! ------------------------------------------------------------------------------
pure module function td_mean(this) result(rst)
    class(t_distribution), intent(in) :: this
    real(real64) :: rst

    real(real64) :: nan
    nan = ieee_value(nan, IEEE_QUIET_NAN)

    if (this%get_dof() < 1.0d0) then
        rst = nan
    else
        rst = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
pure module function td_median(this) result(rst)
    class(t_distribution), intent(in) :: this
    real(real64) :: rst
    rst = 0.0d0
end function

! ------------------------------------------------------------------------------
pure module function td_mode(this) result(rst)
    class(t_distribution), intent(in) :: this
    real(real64) :: rst
    rst = 0.0d0
end function

! ------------------------------------------------------------------------------
pure module function td_variance(this) result(rst)
    class(t_distribution), intent(in) :: this
    real(real64) :: rst

    real(real64) :: inf
    inf = ieee_value(inf, IEEE_POSITIVE_INF)

    if (this%get_dof() > 2.0d0) then
        rst = this%get_dof() / (this%get_dof() - 2.0d0)
    else
        rst = inf
    end if
end function

! ------------------------------------------------------------------------------
module subroutine td_set_params(this, x)
    class(t_distribution), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x
    this%m_dof = x(1)
end subroutine

! ******************************************************************************
! BETA DISTRIBUTION
! ------------------------------------------------------------------------------
pure module function bd_get_alpha(this) result(rst)
    class(beta_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%m_alpha
end function

! --------------------
module subroutine bd_set_alpha(this, x)
    class(beta_distribution), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_alpha = x
end subroutine

! ------------------------------------------------------------------------------
pure module function bd_get_beta(this) result(rst)
    class(beta_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%m_beta
end function

! --------------------
module subroutine bd_set_beta(this, x)
    class(beta_distribution), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_beta = x
end subroutine

! ------------------------------------------------------------------------------
pure elemental module function bd_pdf(this, x) result(rst)
    class(beta_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = beta_distribution_pdf( &
        this%get_alpha(), &
        this%get_beta(), &
        x &
    )
end function

! ------------------------------------------------------------------------------
pure elemental module function bd_cdf(this, x) result(rst)
    class(beta_distribution), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64) :: rst
    rst = beta_distribution_cdf( &
        this%get_beta(), &
        this%get_beta(), &
        x &
    )
end function

! ------------------------------------------------------------------------------
pure module function bd_mean(this) result(rst)
    class(beta_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%get_alpha() / (this%get_alpha() + this%get_beta())
end function

! ------------------------------------------------------------------------------
pure module function bd_median(this) result(rst)
    class(beta_distribution), intent(in) :: this
    real(real64) :: rst
    ! rst = 1.0d0 / regularized_beta(0.5d0, this%get_alpha(), this%get_beta())
    ! Approximation
    rst = (this%get_alpha() - 1.0d0 / 3.0d0) / ( &
        this%get_alpha() + this%get_beta() - 2.0d0 / 3.0d0 &
    )
end function

! ------------------------------------------------------------------------------
pure module function bd_mode(this) result(rst)
    class(beta_distribution), intent(in) :: this
    real(real64) :: rst
    rst = (this%get_alpha() - 1.0d0) / &
        (this%get_alpha() + this%get_beta() - 2.0d0)
end function

! ------------------------------------------------------------------------------
pure module function bd_variance(this) result(rst)
    class(beta_distribution), intent(in) :: this
    real(real64) :: rst
    rst = this%get_alpha() * this%get_beta() / ( &
        (this%get_alpha() + this%get_beta())**2 * &
        (this%get_alpha() + this%get_beta() + 1.0d0) &
    )
end function

! ------------------------------------------------------------------------------
module subroutine bd_set_params(this, x)
    class(beta_distribution), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x
    this%m_alpha = x(1)
    this%m_beta = x(2)
end subroutine

! ------------------------------------------------------------------------------
pure module function bd_geometric_mean(this) result(rst)
    class(beta_distribution), intent(in) :: this
    real(real64) :: rst
    rst = exp(digamma(this%get_alpha()) - &
        digamma(this%get_alpha() + this%get_beta()) &
    )
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule
