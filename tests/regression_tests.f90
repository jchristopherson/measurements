! regression_test.f90

module regression_test
    use iso_fortran_env
    use measurements_core
    implicit none
contains
! ------------------------------------------------------------------------------
function llsq_mimo_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8
    integer(int32), parameter :: dof = 2
    integer(int32), parameter :: npts = 34

    ! Local Variables
    integer(int32) :: i, j
    real(real64) :: xt(npts, dof), yt(npts, dof), x(dof, npts), y(dof, npts), &
        a(dof, dof), ans(dof, dof), delta

    ! Initialization
    rst = .true.
    xt = reshape([0.0d0, 0.38905d0, 0.77816d0, 0.97269d0, 1.16714d0, &
        1.556d0, 1.94484d0, 0.9726d0, -0.00001d0, 0.0d0, -0.388886d0, &
        -0.77775d0, -0.97215d0, -1.16654d0, -1.55533d0, -1.9441d0, -0.97171d0, &
        0.00004d0, 0.0d0, -0.00044d0, -0.0013d0, -0.0024d0, -0.00382d0, &
        -0.00528d0, -0.00257d0, 0.00015d0, 0.0d0, 0.00144d0, 0.00306d0, &
        0.00446d0, 0.00567d0, 0.00688d0, 0.00451d0, -0.00002d0, 0.0d0, &
        0.00122d0, 0.00259d0, 0.0029d0, 0.00314d0, 0.00338d0, 0.00356d0, &
        0.00477d0, -0.00001d0, 0.0d0, 0.00021d0, 0.00051d0, 0.00069d0, &
        0.00088d0, 0.0013d0, 0.00175d0, 0.00058d0, 0.00003d0, 0.0d0, &
        0.27156d0, 0.54329d0, 0.81507d0, 1.08682d0, 1.35881d0, 0.81553d0, &
        0.0001d0, 0.0d0, -0.27145d0, -0.54312d0, -0.81493d0, -1.0868d0, &
        -1.35879d0, -0.81548d0, 0.0d0], [npts, dof])
    yt = reshape([0.0d0, 3000.0d0, 6000.0d0, 7500.0d0, 9000.0d0, 12000.0d0, &
        15000.0d0, 7500.0d0, 0.0d0, 0.0d0, -3000.0d0, -6000.0d0, -7500.0d0, &
        -9000.0d0, -12000.0d0, -15000.0d0, -7500.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 67.7908728d0, 135.5817456d0, 203.3726184d0, &
        271.1634912d0, 338.954364d0, 203.3726184d0, 0.0d0, 0.0d0, &
        -67.7908728d0, -135.5817456d0, -203.3726184d0, -271.1634912d0, &
        -338.954364d0, -203.3726184d0, 0.0d0], [npts, dof])
    x = transpose(xt)
    y = transpose(yt)

    ! Define the solution
    ans = reshape([7713.710356038485d0, -0.2154564584637752d0, &
        33.5206925537885d0, 249.4768502933202d0], [dof, dof])

    ! Compute A
    a = linear_least_squares_mimo(x, y)

    ! Test
    do j = 1, dof
        do i = 1, dof
            delta = ans(i,j) - a(i,j)
            if (abs(delta) > tol) then
                rst = .false.
                print '(A)', "LLSQ_MIMO_TEST FAILED."
                print *, "Expected: ", ans(i,j)
                print *, "Computed: ", a(i,j)
                print *, "Difference (Expected - Computed): ", delta
                print *, "Row: ", i
                print *, "Column: ", j
            end if
        end do
    end do
end function

! ------------------------------------------------------------------------------
function llsq_miso_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8
    integer(int32), parameter :: dof = 2
    integer(int32), parameter :: npts = 18

    ! Local Variables
    integer(int32) :: i
    real(real64) :: xt(npts, dof), x(dof, npts), y(npts), a(dof), ans(dof), &
        delta

    ! Initialization
    rst = .true.
    xt = reshape([0.0d0, 0.38905d0, 0.77816d0, 0.97269d0, 1.16714d0, 1.556d0, &
        1.94484d0, 0.9726d0, -1.0d-05, 0.0d0, -0.388886d0, -0.77775d0, &
        -0.97215d0, -1.16654d0, -1.55533d0, -1.9441d0, -0.97171d0, 4.0d-05, &
        0.0d0, 0.0d0, 0.00122d0, 0.00259d0, 0.0029d0, 0.00314d0, 0.00338d0, &
        0.00356d0, 0.00477d0, -1.0d-05, 0.0d0, 0.00021d0, 0.00051d0, &
        0.00069d0, 0.00088d0, 0.0013d0, 0.00175d0, 0.00058d0, 3.0d-05], &
        [npts, dof])
    x = transpose(xt)
    y = [0.0d0, 3000.d0, 6000.d0, 7500.d0, 9000.d0, 12000.d0, 15000.d0, &
        7500.d0, 0.d0, 0.d0, -3000.d0, -6000.d0, -7500.d0, -9000.d0, &
        -12000.d0, -15000.d0, -7500.d0, 0.0d0]
    ans = [7714.3632680523087d0, -860.96082867646840d0]

    ! Compute A
    a = linear_least_squares_miso(x, y)

    ! Test
    do i = 1, dof
        delta = ans(i) - a(i)
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "LLSQ_MISO_TEST FAILED."
            print *, "Expected: ", ans(i)
            print *, "Computed: ", a(i)
            print *, "Difference (Expected - Computed): ", delta
            print *, "Index: ", i
        end if
    end do
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
