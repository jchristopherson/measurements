! measurements_special.f90

submodule (measurements_core) measurements_special
contains
! ------------------------------------------------------------------------------
pure elemental module function beta(a, b) result(z)
    ! Arguments
    real(real64), intent(in) :: a, b
    real(real64) :: z

    ! Process
    ! REF: https://en.wikipedia.org/wiki/Beta_function
    z = gamma(a) * gamma(b) / gamma(a + b)
end function

! ------------------------------------------------------------------------------
pure elemental module function regularized_beta(x, a, b) result(z)
    ! Arguments
    real(real64), intent(in) :: x, a, b
    real(real64) :: z

    ! Process
    ! REF: https://en.wikipedia.org/wiki/Beta_function
    z = incomplete_beta(x, a, b) / beta(a, b)
end function

! ------------------------------------------------------------------------------
pure elemental module function incomplete_beta(x, a, b) result(z)
    ! Arguments
    real(real64), intent(in) :: x, a, b
    real(real64) :: z

    ! Process
    z = beta_distribution(a, b, x) * beta(a, b)
end function

! ------------------------------------------------------------------------------
end submodule
