! measurements_regression.f90

submodule (measurements_core) measurements_regression
    use linalg_core
contains
! ------------------------------------------------------------------------------
module function linear_least_squares_mimo(x, y, thrsh, err) result(a)
    ! Arguments
    real(real64), intent(in), dimension(:,:) :: x, y
    real(real64), intent(in), optional :: thrsh
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:,:) :: a

    ! Local Variables
    integer(int32) :: m, n, k, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    real(real64), allocatable, dimension(:,:) :: xw, xinv
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = size(y, 1)
    k = size(y, 2)
    n = size(x, 1)

    ! Dimensionallity Check
    if (size(x, 2) /= k) then
        ! ERROR - Size Mismatch:
        call errmgr%report_error("linear_least_squares_mimo", &
            "The number of columns in matrix X must match that of matrix Y.", &
            M_ARRAY_SIZE_ERROR)
        return
    end if
    if (m > n .or. k < n) then
        ! ERROR - underdefined problem
        call errmgr%report_error("linear_least_squares_mimo", &
            "There is insufficient data to fit the problem as posed.", &
            M_UNDERDEFINED_PROBLEM)
        return
    end if

    ! Memory Allocation
    allocate(xinv(k, n), stat = flag)
    if (flag /= 0) then
        ! ERROR:
        call errmgr%report_error("linear_least_squares_mimo", &
            "There is insufficient memory to complete the operation.", &
            M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Compute A = Y * pinv(X)
    xw = x                  ! Copy X to avoid overwriting it
    call mtx_pinverse(xw, xinv, tol = thrsh, err = errmgr)
    if (errmgr%has_error_occurred()) return
    a = matmul(y, xinv)     ! Computes A = Y * pinv(X)
end function

! ------------------------------------------------------------------------------
module function linear_least_squares_miso(x, y, thrsh, err) result(a)
    ! Arguments
    real(real64), intent(in), dimension(:,:) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(in), optional :: thrsh
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:) :: a

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: n, k, flag
    real(real64), allocatable, dimension(:,:) :: xw, xinv
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x, 1)
    k = size(x, 2)

    ! Dimensionality Check
    if (size(y) /= k) then
        ! ERROR: Size mismatch
        call errmgr%report_error("linear_least_squares_miso", &
            "The number of columns in matrix X must match that of matrix Y.", &
            M_ARRAY_SIZE_ERROR)
        return
    end if

    if (k < n) then
        ! ERROR: Underdefined problem
        call errmgr%report_error("linear_least_squares_miso", &
            "There is insufficient data to fit the problem as posed.", &
            M_UNDERDEFINED_PROBLEM)
        return
    end if

    ! Memory Allocation
    allocate(xinv(k, n), stat = flag)
    if (flag /= 0) then
        ! ERROR:
        call errmgr%report_error("linear_least_squares_miso", &
            "There is insufficient memory to complete the operation.", &
            M_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Compute the pseudo-inverse of X
    xw = x      ! Prevent the inversion routine from overwriting X
    call mtx_pinverse(xw, xinv, tol = thrsh, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute A = Y * pinv(X)
    a = matmul(y, xinv)
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule
