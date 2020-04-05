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
pure elemental module function incomplete_beta(x, a, b) result(z)
    ! Arguments
    real(real64), intent(in) :: x, a, b
    real(real64) :: z

    ! Process
    ! REF: https://en.wikipedia.org/wiki/Beta_function
    z = beta_distribution(a, b, x) * beta(a, b)
end function

! ------------------------------------------------------------------------------
! incomplete gamma function

! ------------------------------------------------------------------------------
! REF: ! https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
pure elemental module function psi(x) result(ps)
    ! Arguments
    real(real64), intent(in) :: x
    real(real64) :: ps

    ! Parameters
    real(real64), parameter :: a1 = -0.83333333333333333d-01
    real(real64), parameter :: a2 =  0.83333333333333333d-02
    real(real64), parameter :: a3 = -0.39682539682539683d-02
    real(real64), parameter :: a4 =  0.41666666666666667d-02
    real(real64), parameter :: a5 = -0.75757575757575758d-02
    real(real64), parameter :: a6 =  0.21092796092796093d-01
    real(real64), parameter :: a7 = -0.83333333333333333d-01
    real(real64), parameter :: a8 =  0.4432598039215686d+00
    real(real64), parameter :: el = 0.5772156649015329d+00
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: k, n
    real(real64) :: x2, s, xa

    ! Process
    xa = abs(x)
    s = 0.0d0

    ! Quick Return
    if (x == aint(x) .and. x <= 0.0d0) then
        ps = 1.0d300
        return
    else if (xa == aint(xa)) then
        n = int(xa)
        do k = 1, n - 1
            s = s + 1.0d0 / real(k, real64)
        end do
        ps = -el + s
    else if (xa + 0.5d0 == aint(xa + 0.5d0)) then
        n = int(xa - 0.5d0)
        do k = 1, n
            s = s + 1.0d0 / (2.0d0 * k - 1.0d0)
        end do
        ps = -el + 2.0d0 * s - 1.386294361119891d0
    else
        if (xa < 1.0d1) then
            n = 10 - int(xa)
            do k = 0, n - 1
                s = s + 1.0d0 / (xa + k)
            end do
            xa = xa + n
        end if
        x2 = 1.0d0 / (xa * xa)
        ps = log(xa) - 0.5d0 / xa + x2 * ((((((( &
             a8 &
            * x2 + a7) &
            * x2 + a6) &
            * x2 + a5) &
            * x2 + a4) &
            * x2 + a3) &
            * x2 + a2) &
            * x2 + a1)
        ps = ps - s
    end if

    if (x < 0.0d0) then
        ps = ps - pi * cos(pi * x) / sin(pi * x) - 1.0d0 / x
    end if
end function

! ------------------------------------------------------------------------------
! hypergeometric function
! Need PSI before HYGFX
! https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
! https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
pure elemental function hypergeom(a, b, c, x) result(hf)
    ! Arguments
    real(real64), intent(in) :: a, b, c, x
    real(real64) :: hf

    ! Parameters
    real(real64), parameter :: el = 0.5772156649015329d0
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    real(real64) :: a0, aa, bb, c0, c1, f0, f1, g0, g1, g2, g3, ga, gabc, gam, &
        gb, gbm, gc, gca, gcab, gcb, gm, hw, pa, pb, r, r0, r1, rm, rp, sm, &
        sp, sp0, x1
    integer(int32) :: j, k, m, nm
    logical :: l0, l1, l2, l3, l4, l5

    ! Initialization
    l0 = (c == aint(c)) .and. (c < 0.0d0)
    l1 = (1.0d0 - x < 1.0d-15) .and. (c - a - b <= 0.0d0)
    l2 = (a == aint(a)) .and. (a < 0.0d0)
    l3 = (b == aint(b)) .and. (b < 0.0d0)
    l4 = (c - a == aint(c - a)) .and. (c - a <= 0.0d0)
    l5 = (c - b == aint(c - b)) .and. (c - b <= 0.0d0)
end function

! ------------------------------------------------------------------------------
end submodule
