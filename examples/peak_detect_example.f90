! peak_detect_example.f90

program main
    use iso_fortran_env
    use measurements_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: dt = 1.0d-3
    real(real64), parameter :: threshold = 1.0d-2

    ! Variables
    real(real64) :: t(npts), x(npts)
    type(peak_info) :: pks
    integer(int32) :: i

    ! Create a waveform
    do i = 1, npts
        if (i == 1) then
            t(i) = 0.0d0
        else
            t(i) = t(i-1) + dt
        end if
        x(i) = exp(-2.0d0 * t(i)) * sin(15.0d0 * t(i))
    end do

    ! Locate the peaks
    pks = peak_detect(x, threshold)

    ! Print the peak and valley information
    print '(A)', "Peaks (indices, t, x)"
    do i = 1, size(pks%max_values)
        print '(AI0AF7.5AF7.5)', achar(9), &
            pks%max_value_indices(i), achar(9), &
            t(pks%max_value_indices(i)), achar(9), &
            pks%max_values(i)
    end do

    print '(A)', "Valleys (indices, t, x)"
    do i = 1, size(pks%min_values)
        print '(AI0AF7.5AF8.5)', achar(9), &
            pks%min_value_indices(i), achar(9), &
            t(pks%min_value_indices(i)), achar(9), &
            pks%min_values(i)
    end do
end program
