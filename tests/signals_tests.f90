! signals_tests.f90

module signals_tests
    use iso_fortran_env
    use measurements_core
    implicit none
contains
! ------------------------------------------------------------------------------
    function fourier_transform_test() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter :: npts = 1024
        real(real64), parameter :: fs = 1024.0d0
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        real(real64), parameter :: f1 = 5.0d0
        real(real64), parameter :: f2 = 50.0d0
        real(real64), parameter :: x1 = 3.0d0
        real(real64), parameter :: x2 = 0.5d0
        integer(int32), parameter :: ind1 = 6
        integer(int32), parameter :: ind2 = 51
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        integer(int32) :: i
        real(real64) :: t, dt, x(npts), delta
        complex(real64), allocatable, dimension(:) :: f

        ! Initialization
        rst = .true.
        dt = 1.0d0 / fs
        t = 0.0d0
        x(1) = 0.0d0
        do i = 1, npts
            t = t + dt
            x(i) = x1 * sin(2.0d0 * pi * f1 * t) + x2 * sin(2.0d0 * pi * f2 * t)
        end do

        ! Compute the Fourier transform of X
        f = fourier_transform(x)

        ! Check
        if (size(f) /= npts / 2 + 1) then
            rst = .false.
            print '(A)', "FOURIER_TRANSFORM FAILED."
            print '(AI0AI0A)', "Expected an array of length ", npts / 2 + 1, &
                ", but found an array of length ", size(f), "."
        end if

        delta = abs(f(ind1)) - x1
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "FOURIER_TRANSFORM FAILED."
            print *, "Expected: ", x1
            print *, "Computed: ", abs(f(ind1))
            print *, "Difference (Expected - Computed): ", delta
        end if

        delta = abs(f(ind2)) - x2
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "FOURIER_TRANSFORM FAILED."
            print *, "Expected: ", x2
            print *, "Computed: ", abs(f(ind2))
            print *, "Difference (Expected - Computed): ", delta
        end if
    end function

! ------------------------------------------------------------------------------
    function periodogram_test() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: fs = 1.024d3
        real(real64), parameter :: x1 = 1.8d0
        real(real64), parameter :: f1 = 1.0d2
        real(real64), parameter :: max_t = 1.0d0
        real(real64), parameter :: ans = x1**2 / 2.0d0
        real(real64), parameter :: tol = 1.0d-8
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

        ! Local Variables
        integer(int32) :: i, n
        real(real64) :: t, dt, pmax, delta
        real(real64), allocatable, dimension(:) :: x, p
        procedure(window_function), pointer :: window

        ! Initialization
        rst = .true.
        dt = 1.0d0 / fs
        n = int(max_t / dt, int32) + 1
        allocate(x(n))
        t = 0.0d0
        x(1) = 0.0d0
        do i = 2, n
            t = t + dt
            x(i) = x1 * sin(2.0d0 * pi * f1 * t)
        end do
        
        ! Define the window function to utilize
        window => rectangular_window

        ! Compute the PSD - use a 256 point FFT
        p = periodogram(x, window, 256)
        
        ! Compute the power at the primary harmonic, and compare to the solution
        pmax = maxval(p)
        delta = ans - pmax
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "PERIODOGRAM_TEST FAILED."
            print *, "Expected: ", ans
            print *, "Computed: ", pmax
            print *, "Difference (Expected - Computed): ", delta
        end if
    end function

! ------------------------------------------------------------------------------
    function fourier_frequency_test() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter :: npts = 1024
        integer(int32), parameter :: nodd = 1025
        real(real64), parameter :: fs = 1024.0d0

        ! Local Variables
        integer(int32) :: m
        real(real64) :: f

        ! Initialization
        rst = .true.
        m = npts / 2 + 1

        ! Test 1 - Check DC - Even Valued
        f = fourier_frequency(fs, 0, npts)
        if (f /= 0.0d0) then
            rst = .false.
            print '(A)', "FOURIER_FREQUENCY_TEST FAILED."
            print *, "Expected: ", 0.0d0
            print *, "Computed: ", f
        end if

        ! Test 2 - Check Nyquist - Even Valued
        f = fourier_frequency(fs, m - 1, npts)
        if (f /= 0.5d0 * fs) then
            rst = .false.
            print '(A)', "FOURIER_FREQUENCY_TEST - EVEN FAILED."
            print *, "Expected: ", 0.5d0 * fs
            print *, "Computed: ", f
        end if

        ! Test 3 - Check Nyquist - Odd Valued
        m = (nodd + 1) / 2
        f = fourier_frequency(fs, m - 1, nodd)
        if (f /= 0.5d0 * fs) then
            rst = .false.
            print '(A)', "FOURIER_FREQUENCY_TEST - ODD FAILED."
            print *, "Expected: ", 0.5d0 * fs
            print *, "Computed: ", f
        end if
    end function

! ------------------------------------------------------------------------------
    function low_pass_filter_test() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: fs = 1024.0d0
        integer(int32), parameter :: npts = 2048
        real(real64), parameter :: cutoff = 100.0d0
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: x(npts), y(npts)
        integer(int32) :: i, ind, m
        complex(real64), allocatable, dimension(:) :: xfrm
        real(real64), allocatable, dimension(:) :: mag

        ! Initialization
        rst = .true.

        ! Construct a random noise signal, and then bias
        call random_number(x)
        x = x - 0.5d0

        ! Filter the data
        y = low_pass_filter(x, fs, cutoff)

        ! Find the index of the cutoff frequency.  Anything above this
        ! value should return a zero value magnitude in its Fourier
        ! transform
        if (mod(npts, 2) == 0) then
            m = npts / 2 + 1
        else
            m = (npts + 1) / 2
        end if
        ind = int(2.0d0 * cutoff * m / fs, int32) + 1 ! +1 accounts for the zero offset

        ! Compute the FFT of the filtered signal, and make sure any value
        ! after the filter frequency is sufficiently close to zero in 
        ! magnitude
        xfrm = fourier_transform(y)
        mag = abs(xfrm)

        ! do i = 1, size(mag)
        !     print *, mag(i), ",", fourier_frequency(fs, i, npts)
        ! end do

        ! Test
        do i = ind + 1, size(mag)
            if (mag(i) > tol) then
                rst = .false.
                print '(A)', "LOW_PASS_FILTER_TEST FAILED."
                print *, "|X| = ", mag(i)
                print *, "Frequency = ", fourier_frequency(fs, i - 1, npts)
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    function fft_test() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter:: npts = 1024
        real(real64), parameter :: fs = 1024.0d0
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        real(real64), parameter :: f1 = 5.0d0
        real(real64), parameter :: f2 = 50.0d0
        real(real64), parameter :: x1 = 3.0d0
        real(real64), parameter :: x2 = 0.5d0
        integer(int32), parameter :: ind1 = 6
        integer(int32), parameter :: ind2 = 51
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        integer(int32) :: i
        real(real64) :: t, dt, x(npts), delta
        complex(real64), allocatable, dimension(:) :: f, finv

        ! Initialization
        rst = .true.
        dt = 1.0d0 / fs
        t = 0.0d0
        x(1) = 0.0d0
        do i = 1, npts
            t = t + dt
            x(i) = x1 * sin(2.0d0 * pi * f1 * t) + x2 * sin(2.0d0 * pi * f2 * t)
        end do

        ! Compute the FFT of X
        f = fft(x)

        ! Test the FFT
        delta = 2.0d0 * abs(f(ind1)) - x1
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "FFT_TRANSFORM FAILED."
            print *, "Expected: ", x1
            print *, "Computed: ", 2.0d0 * abs(f(ind1))
            print *, "Difference (Expected - Computed): ", delta
        end if

        delta = 2.0d0 * abs(f(ind2)) - x2
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "FFT_TRANSFORM FAILED."
            print *, "Expected: ", x2
            print *, "Computed: ", 2.0d0 * abs(f(ind2))
            print *, "Difference (Expected - Computed): ", delta
        end if

        ! Compute the inverse transform
        finv = ifft(f)

        ! Compare the results
        do i = 1, npts
            delta = real(finv(i)) - x(i)
            if (abs(delta) > tol) then
                rst = .false.
                print '(A)', "FFT_TRANFORM FAILED."
                print *, "Expected: ", x(i)
                print *, "Computed: ", real(finv(i))
                print *, "Difference (Expected - Computed): ", delta
                print *, "Index: ", i
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    function finite_diff_test() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter:: npts = 1024
        real(real64), parameter :: fs = 1024.0d0
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        real(real64), parameter :: f1 = 5.0d0
        real(real64), parameter :: f2 = 50.0d0
        real(real64), parameter :: x1 = 3.0d0
        real(real64), parameter :: x2 = 0.5d0
        real(real64), parameter :: tol = 1.0d-1

        ! Local Variables
        integer(int32) :: i
        real(real64) :: w1, w2, dt, delta, maxX, scaledTol, &
            t(npts), x(npts), dxdt(npts), ans(npts), xi(npts)

        ! Initialization
        rst = .true.
        dt = 1.0d0 / fs
        t(1) = 0.0d0
        do i = 2, npts
            t(i) = t(i-1) + dt
        end do
        w1 = 2.0d0 * pi * f1
        w2 = 2.0d0 * pi * f2
        x = x1 * sin(w1 * t) + x2 * sin(w2 * t)
        ans = w1 * x1 * cos(w1 * t) + w2 * x2 * cos(w2 * t)
        
        ! Utilize finite_difference
        dxdt = finite_difference(t, x)

        ! Test
        maxX = maxval(ans)
        scaledTol = tol * maxX
        do i = 1, npts
            delta = ans(i) - dxdt(i)
            if (abs(delta) > scaledTol) then
                rst = .false.
                print '(A)', "FINITE_DIFF_TEST FAILED."
                print *, "Expected: ", ans(i)
                print *, "Computed: ", dxdt(i)
                print *, "Difference (Expected - Computed): ", delta
                print *, "Index: ", i
            end if
        end do

        ! Integrate the computed dydx to see if we arrive back at x(t)
        xi = integrate(t, ans, 0.0d0) ! initial value at t = 0 is 0.

        ! Test
        maxX = maxval(x)
        scaledTol = tol * maxX
        do i = 1, npts
            delta = x(i) - xi(i)
            if (abs(delta) > scaledTol) then
                rst = .false.
                print '(A)', "FINITE_DIFF_TEST FAILED - INTEGRATION CHECK."
                print *, "Expected: ", x(i)
                print *, "Computed: ", xi(i)
                print *, "Difference (Expected - Computed): ", delta
                print *, "Index: ", i
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    function remove_offset_test() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter:: npts = 1024
        real(real64), parameter :: fs = 1024.0d0
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        real(real64), parameter :: f1 = 5.0d0
        real(real64), parameter :: f2 = 50.0d0
        real(real64), parameter :: x1 = 3.0d0
        real(real64), parameter :: x2 = 0.5d0
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        integer(int32) :: i
        real(real64) :: w1, w2, dt, delta, &
            t(npts), x(npts), ans(npts), xi(npts)

        ! Initialization
        rst = .true.
        dt = 1.0d0 / fs
        t(1) = 0.0d0
        do i = 2, npts
            t(i) = t(i-1) + dt
        end do
        w1 = 2.0d0 * pi * f1
        w2 = 2.0d0 * pi * f2
        ans = x1 * sin(w1 * t) + x2 * sin(w2 * t)
        x = ans + 0.5d0 ! add in a DC offset

        ! Remove the offset
        xi = remove_dc_offset(x)

        ! Test
        do i = 1, npts
            delta = ans(i) - xi(i)
            if (abs(delta) > tol) then
                rst = .false.
                print '(A)', "REMOVE_OFFSET_TEST FAILED."
                print *, "Expected: ", ans(i)
                print *, "Computed: ", xi(i)
                print *, "Difference (Expected - Computed): ", delta
                print *, "Index: ", i
            end if
        end do
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
