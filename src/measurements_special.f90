! measurements_special.f90

submodule (measurements_core) measurements_special
    use ieee_arithmetic
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
    real(real64) :: c, tol

    ! Process
    tol = sqrt(epsilon(tol))
    c = (x**a) * ((1.0d0 - x)**b)
    if (x == 0.0d0) then
        z = 0.0d0
    else if (abs(1.0d0 - x) < tol) then
        ! The incomplete beta function simply reduces to the beta function
        z = beta(a, b)
    else
        z = c * inc_beta_partial_fraction(x, a, b)
    end if
end function

! ------------------------------------------------------------------------------
! REF:
! - https://people.sc.fsu.edu/~jburkardt/f_src/asa103/asa103.f90
pure elemental module function digamma(x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Parameters
    real(real64), parameter :: c = 8.5d0
    real(real64), parameter :: euler_mascheroni = 0.57721566490153286060d0
    
    ! Local Variables
    real(real64) :: r, x, x2, nan

    ! If x <= 0.0
    if (x <= 0.0) then
        nan = ieee_value(nan, IEEE_QUIET_NAN)
        rst = nan
        return
    end if

    ! Approximation for a small argument
    if (x <= 1.0d-6) then
        rst = -euler_mascheroni - 1.0d0 / x + 1.6449340668482264365d0 * x
        return
    end if

    ! Process
    rst = 0.0d0
    x2 = x
    do while (x2 < c)
        rst = rst - 1.0d0 / x2
        x2 = x2 + 1.0d0
    end do

    r = 1.0d0 / x2
    rst = rst + log(x2) - 0.5d0 * r
    r = r * r
    rst = rst &
        -r * (1.0d0 / 12.0d0 &
        - r * (1.0d0 / 120.0d0 &
        - r * (1.0d0 / 252.0d0 &
        - r * (1.0d0 / 240.0d0 &
        - r * (1.0d0 / 132.0d0) &
    ))))
end function




! !
! !  Reduce to DIGAMA(X + N).
! !
!   digamma = 0.0D+00
!   x2 = x

!   do while ( x2 < c )
!     digamma = digamma - 1.0D+00 / x2
!     x2 = x2 + 1.0D+00
!   end do
! !
! !  Use Stirling's (actually de Moivre's) expansion.
! !
!   r = 1.0D+00 / x2

!   digamma = digamma + log ( x2 ) - 0.5D+00 * r

!   r = r * r

!   digamma = digamma &
!     - r * ( 1.0D+00 / 12.0D+00 &
!     - r * ( 1.0D+00 / 120.0D+00 &
!     - r * ( 1.0D+00 / 252.0D+00 &
!     - r * ( 1.0D+00 / 240.0D+00 &
!     - r * ( 1.0D+00 / 132.0D+00 ) ) ) ) )

!   return
! end

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
