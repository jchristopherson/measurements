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
    z = incomplete_beta(x, a, b) / beta(a, b)
end function

! ------------------------------------------------------------------------------
! REF:
! - https://github.com/dcwuser/metanumerics/blob/master/Numerics/Functions/AdvancedMath_Gamma.cs
pure elemental module function incomplete_beta(x, a, b) result(z)
    ! Arguments
    real(real64), intent(in) :: x, a, b
    real(real64) :: z

    ! Local Variables
    ! real(real64) :: xtp

    ! Process
    ! xtp = (a + 1.0d0) / (a + b + 2.0d0)
    ! if (x > xtp) then
    !     z = beta(a, b) - beta_distribution(b, a, 1.0d0 - x)
    ! else
    !     z = (x**a) * ((1.0d0 - x)**b) * inc_beta_partial_fraction(x, a, b)
    ! end if
    if (x == 0.0d0) then
        z = 0.0d0
    else if (x == 1.0d0) then
        ! The incomplete beta function simply reduces to the beta function
        z = beta(a, b)
    else
        z = (x**a) * ((1.0d0 - x)**b) * inc_beta_partial_fraction(x, a, b)
    end if
end function

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
pure elemental function inc_beta_partial_fraction(x, a, b) result(z)
    ! Arguments
    real(real64), intent(in) :: x, a, b
    real(real64) :: z

    ! Parameters
    integer(int32), parameter :: max_iter = 100000
    real(real64), parameter :: tol = 1.0d-12

    ! Local Variables
    integer(int32) :: m, k
    real(real64) :: ab, p, d, df, f, fold

    ! Initialization
    ab = a + b
    p = -x * ab / (a + 1.0d0)
    d = 1.0d0 / (1.0d0 + p)
    df = -p * d
    f = 1.0d0 + df

    ! Process
    do k = 2, max_iter
        fold = f
        m = k / 2
        p = x / ((a + (k - 1)) * (a + k))
        if (mod(k, 2) == 0) then
            p = p * m * (b - m)
        else
            p = p * (-(a + m) * (ab + m))
        end if
        d = 1.0d0 / (1.0d0 + p * d)
        df = (d - 1.0d0) * df
        f = f + df
        if (abs(f - fold) <= tol) then
            z = f / a
            exit
        end if
    end do
end function

! ------------------------------------------------------------------------------
end submodule
