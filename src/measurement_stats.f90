! measurement_stats.f90

submodule (measurements_core) measurement_stats
    use linalg_core
contains
! ------------------------------------------------------------------------------
pure module function is_monotonic(x) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    logical :: rst

    ! Local Variables
    integer(int32) :: i
    logical :: ascend

    ! Quick Return
    rst = .true.
    if (size(x) <= 1) return
    if (x(2) == x(1)) then
        rst = .false.
        return
    end if

    ! Process
    ascend = x(2) > x(1)
    if (ascend) then
        do i = 2, size(x)
            if (x(i) <= x(i-1)) then
                rst = .false.
                exit
            end if
        end do
    else
        do i = 2, size(x)
            if (x(i) >= x(i-1)) then
                rst = .false.
                exit
            end if
        end do
    end if
end function

! ------------------------------------------------------------------------------
pure module function mean(x) result(z)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64) :: z

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(x)
    if (n == 0) then
        z = zero
    else
        z = x(1)
        do i = 2, n
            z = z + (x(i) - z) / i
        end do
    end if
end function

! ------------------------------------------------------------------------------
module function median(x) result(z)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64) :: z

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: half = 0.5d0

    ! Local Variables
    integer(int32) :: n, nmid, nmidp1, iflag
    real(real64), allocatable, dimension(:) :: y

    ! Initialization
    n = size(x)
    nmid = n / 2
    nmidp1 = nmid + 1
    iflag = n - 2 * nmid

    ! Sort the data into ascending order
    y = x
    call sort(y, .true.)

    ! Compute the median
    if (iflag == 0) then
        z = half * (y(nmid) + y(nmidp1))
    else
        z = y(nmidp1)
    end if
end function

! ------------------------------------------------------------------------------
pure module function variance(x) result(v)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64) :: v

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: oldMean, newMean

    ! Process
    n = size(x)
    if (n <= 1) then
        v = zero
    else
        oldMean = x(1)
        v = zero
        do i = 2, n
            newMean = oldMean + (x(i) - oldMean) / i
            v = v + (x(i) - oldMean) * (x(i) - newMean)
            oldMean = newMean
        end do
        v = v / (n - one)
    end if
end function

! ------------------------------------------------------------------------------
pure module function standard_deviation(x) result(s)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    real(real64) :: s

    ! Process
    s = sqrt(variance(x))
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule
