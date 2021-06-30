! stats_tests.f90

module stats_tests
    use iso_fortran_env
    use measurements_core
    use linalg_core
    implicit none
contains
! ------------------------------------------------------------------------------
function mean_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 100000
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), ans, computed, delta

    ! Initialization
    rst = .true.
    call random_number(x)

    ! Compute the actual solution the old-school way
    ans = sum(x) / npts

    ! Utilize the actual method
    computed = mean(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "MEAN_TEST FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function median_test_even() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), x1, x2, ans, computed, delta

    ! Initialization
    rst = .true.

    ! Populate an array with random values, and then sort into ascending
    ! order
    call random_number(x)
    call sort(x, .true.)

    ! Compute the answer
    x1 = x(npts / 2)
    x2 = x(npts / 2 + 1)
    ans = 0.5d0 * (x1 + x2)

    ! Process
    computed = median(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "MEDIAN_TEST_EVEN FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function median_test_odd() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 10001
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), ans, computed, delta

    ! Initialization
    rst = .true.

    ! Populate an array with random values, and then sort into ascending
    ! order
    call random_number(x)
    call sort(x, .true.)

    ! Compute the answer
    ans = x(floor(npts / 2.0d0) + 1)

    ! Process
    computed = median(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "MEDIAN_TEST_ODD FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function variance_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), ans, computed, delta, avg

    ! Initialization
    rst = .true.

    ! Create a random data set
    call random_number(x)

    ! Compute the solution
    avg = mean(x)
    ans = sum((x - avg)**2) / (npts - 1.0d0)

    ! Process
    computed = variance(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "VARIANCE_TEST FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function std_dev_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), ans, computed, delta, avg

    ! Initialization
    rst = .true.

    ! Create a random data set
    call random_number(x)

    ! Compute the solution
    avg = mean(x)
    ans = sqrt(sum((x - avg)**2) / (npts - 1.0d0))

    ! Process
    computed = standard_deviation(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "STANDARD_DEVIATION_TEST FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function range_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: x(npts), ans, computed, delta

    ! Initialization
    rst = .true.

    ! Create a random data set
    call random_number(x)

    ! Compute the solution
    ans = maxval(x) - minval(x)

    ! Process
    computed = data_range(x)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "RANGE_TEST FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function z_score_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameter
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: tol = 1.0d-8
    real(real64), parameter :: z80 = 1.281551565545d0
    real(real64), parameter :: z90 = 1.644853626951d0
    real(real64), parameter :: z95 = 1.959963984540d0
    real(real64), parameter :: z99 = 2.575829303549d0
    real(real64), parameter :: z999 = 3.290526731492d0

    ! Local Variables
    real(real64) :: ans, computed, delta

    ! Initialization
    rst = .true.

    ! 80%
    ans = z80
    computed = z_score(0.8d0)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "Z_SCORE_TEST 80% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! 90%
    ans = z90
    computed = z_score(0.9d0)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "Z_SCORE_TEST 90% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! 95%
    ans = z95
    computed = z_score(0.95d0)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "Z_SCORE_TEST 95% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! 99%
    ans = z99
    computed = z_score(0.99d0)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "Z_SCORE_TEST 99% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! 99.9%
    ans = z999
    computed = z_score(0.999d0)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "Z_SCORE_TEST 99.9% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function t_score_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameter
    integer(int32), parameter :: npts = 21
    real(real64), parameter :: tol = 1.0d-3
    real(real64), parameter :: t80 = 0.860d0
    real(real64), parameter :: t85 = 1.064d0
    real(real64), parameter :: t90 = 1.325d0
    real(real64), parameter :: t95 = 1.725d0

    ! Local Variables 
    real(real64) :: ans, computed, delta

    ! Initialization
    rst = .true.

    ! 80 %
    ans = t80
    computed = t_score(0.8d0, npts)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "T_SCORE_TEST 80% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! 85%
    ans = t85
    computed = t_score(0.85d0, npts)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "T_SCORE_TEST 85% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! 90%

    ans = t90
    computed = t_score(0.9d0, npts)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "T_SCORE_TEST 90% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! 95%
    ans = t95
    computed = t_score(0.95d0, npts)
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "T_SCORE_TEST 95% FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function confidence_interval_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: alpha = 0.95d0
    integer(int32), parameter :: npts = 10
    real(real64), parameter :: ans = 0.196943590756784000d0
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: zval, delta, computed, x(npts)

    ! Initialization
    rst = .true.
    x = [0.1266904993777120d0, 0.3827711768054340d0, 0.1370953805850570d0, &
        0.5852213531153070d0, 0.2267533281658030d0, 0.0999861358308985d0, &
        0.5851003510284570d0, 0.8136628645855180d0, 0.7400357894369070d0, &
        0.9787774755208680d0]

    ! Process
    zval = z_score(alpha)
    computed = confidence_interval(x, zval)

    ! Test
    delta = ans - computed
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "CONFIDENCE_INTERVAL_TEST FAILED."
        print *, "Expected: ", ans
        print *, "Computed: ", computed
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function normal_distribution_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: avg = 5.264049719320700d-1
    real(real64), parameter :: sigma = 2.110340893383650d-1
    integer(int32), parameter :: npts = 21
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: delta, x(npts), ans(npts), computed(npts)
    integer(int32) :: i

    ! Initialization
    rst = .true.
    x = [-5.0d0, -4.5d0, -4.0d0, -3.5d0, -3.0d0, -2.5d0, -2.0d0, -1.5d0, &
        -1.0d0, -0.5d0, 0.0d0, 0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, &
        3.5d0, 4.0d0, 4.5d0, 5.0d0]
    ans = [2.306272166314290d-149, 1.229690115816890d-123, &
        2.392009083889280d-100, 1.697508606279160d-79, 4.394840953499330d-61, &
        4.151034651360500d-45, 1.430380472227950d-31, 1.798161951816050d-20, &
        8.246849202359420d-12, 1.379841874725850d-5, 8.422724695133530d-02, &
        1.875676375640870d0, 1.523860560727160d-1, 4.516630549640910d-5, &
        4.883890565618130d-11, 1.926634346379860d-19, 2.772775386838060d-30, &
        1.455834812310320d-43, 2.788634062219750d-59, 1.948735843727060d-77, &
        4.968169870317840d-98]

    ! Process
    computed = normal_distribution_pdf(avg, sigma, x)

    ! Test
    do i = 1, npts
        delta = ans(i) - computed(i)
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "NORMAL_DISTRIBUTION_TEST FAILED."
            print *, "Expected: ", ans(i)
            print *, "Computed: ", computed(i)
            print *, "Difference (Expected - Computed): ", delta
            print *, "Index: ", i
        end if
    end do
end function

! ------------------------------------------------------------------------------
function t_distribution_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: dof = 20.0d0
    integer(int32), parameter :: npts = 21
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: delta, x(npts), ans(npts), computed(npts)
    integer(int32) :: i

    ! Initialization
    rst = .true.
    x = [-5.0d0, -4.5d0, -4.0d0, -3.5d0, -3.0d0, -2.5d0, -2.0d0, -1.5d0, &
        -1.0d0, -0.5d0, 0.0d0, 0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, &
        3.5d0, 4.0d0, 4.5d0, 5.0d0]
    ans = [0.0000789891062440353d0, 0.0002548336678335860d0, &
        0.0008224743001331390d0, 0.0026105772275963500d0, &
        0.0079637866461806600d0, 0.0226694437191449000d0, &
        0.0580872152473570000d0, 0.1286273829721460000d0, &
        0.2360456491267010000d0, 0.3458086123837420000d0, &
        0.3939885857114330000d0, 0.3458086123837420000d0, &
        0.2360456491267010000d0, 0.1286273829721460000d0, &
        0.0580872152473570000d0, 0.0226694437191449000d0, &
        0.0079637866461806600d0, 0.0026105772275963500d0, &
        0.0008224743001331390d0, 0.0002548336678335860d0, &
        0.0000789891062440353d0]

    ! Process
    computed = t_distribution_pdf(dof, x)

    ! Test
    do i = 1, npts
        delta = ans(i) - computed(i)
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "T_DISTRIBUTION_TEST FAILED."
            print *, "Expected: ", ans(i)
            print *, "Computed: ", computed(i)
            print *, "Difference (Expected - Computed): ", delta
            print *, "Index: ", i
        end if
    end do
end function

! ------------------------------------------------------------------------------
function beta_distribution_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: a = 1.0d0
    real(real64), parameter :: b = 3.0d0
    integer(int32), parameter :: npts = 18
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: delta, x(npts), ans(npts), computed(npts)
    integer(int32) :: i

    ! Initialization
    rst = .true.
    x = [0.1d0, 0.15d0, 0.2d0, 0.25d0, 0.3d0, 0.35d0, 0.4d0, 0.45d0, 0.5d0, &
        0.55d0, 0.6d0, 0.65d0, 0.7d0, 0.75d0, 0.8d0, 0.85d0, 0.9d0, 0.95d0]
    ans = [2.43000d0, 2.16750d0, 1.92000d0, 1.68750d0, 1.47000d0, 1.26750d0, &
        1.08000d0, 0.90750d0, 0.75000d0, 0.60750d0, 0.48000d0, 0.36750d0, &
        0.27000d0, 0.18750d0, 0.12000d0, 0.06750d0, 0.03000d0, 0.00750d0]

    ! Process
    computed = beta_distribution_pdf(a, b, x)

    ! Test
    do i = 1, npts
        delta = ans(i) - computed(i)
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "BETA_DISTRIBUTION_TEST FAILED."
            print *, "Expected: ", ans(i)
            print *, "Computed: ", computed(i)
            print *, "Difference (Expected - Computed): ", delta
            print *, "Index: ", i
        end if
    end do
end function

! ------------------------------------------------------------------------------
function f_distribution_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: d1 = 10.0d0
    real(real64), parameter :: d2 = 1.0d0
    integer(int32), parameter :: npts = 18
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    real(real64) :: delta, x(npts), ans(npts), computed(npts)
    integer(int32) :: i

    ! Initialization
    rst = .true.
    x = [0.1d0, 0.15d0, 0.2d0, 0.25d0, 0.3d0, 0.35d0, 0.4d0, 0.45d0, 0.5d0, &
        0.55d0, 0.6d0, 0.65d0, 0.7d0, 0.75d0, 0.8d0, 0.85d0, 0.9d0, 0.95d0]
    ans = [0.27189774911348000d0, 0.40342757249598100d0, &
        0.46776063476011400d0, 0.48916613056974600d0, 0.48666000366211000d0, &
        0.47170875848615800d0, 0.45079130426395800d0, 0.42748989305298300d0, &
        0.40375575782592800d0, 0.38062550389189700d0, 0.35862153897298900d0, &
        0.33797693751338100d0, 0.31876293731516800d0, 0.30096192124722400d0, &
        0.28450947518162900d0, 0.26931867071908400d0, 0.25529401072011300d0, &
        0.24233930873721800d0]

    ! Process
    computed = f_distribution_pdf(d1, d2, x)

    ! Test
    do i = 1, npts
        delta = ans(i) - computed(i)
        if (abs(delta) > tol) then
            rst = .false.
            print '(A)', "F_DISTRIBUTION_TEST FAILED."
            print *, "Expected: ", ans(i)
            print *, "Computed: ", computed(i)
            print *, "Difference (Expected - Computed): ", delta
            print *, "Index: ", i
        end if
    end do
end function

! ------------------------------------------------------------------------------
function t_test_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n1 = 15
    integer(int32), parameter :: n2 = 12
    real(real64), parameter :: tol = 1.0d-8
    real(real64), parameter :: ans1 = 0.6574262172637550d0
    real(real64), parameter :: ans2 = 0.6621143019183570d0
    real(real64), parameter :: ans3 = 0.2286303186286610d0

    ! Local Variables
    real(real64) :: x1(n1), x2(n2), x3(n1), delta
    type(statistic) :: t

    ! Initialization
    rst = .true.
    x1 = [0.5875701638177520d0, &
        0.4925038626867790d0, &
        0.0997847738770978d0, &
        0.5540924574002610d0, &
        0.0833121626712929d0, &
        0.1738451549308330d0, &
        0.3521655274264620d0, &
        0.2239625528107020d0, &
        0.6871828620071030d0, &
        0.3248518075223050d0, &
        0.0551977898473518d0, &
        0.8648552295498370d0, &
        0.9239272586628300d0, &
        0.4917627939852090d0, &
        0.3508690262031490d0]
    x2 = [0.7557955531972870d0, &
        0.0482975843398515d0, &
        0.7609889442453010d0, &
        0.4898203045069780d0, &
        0.4382872343657070d0, &
        0.9676872466629530d0, &
        0.1167483190258670d0, &
        0.0399776180777329d0, &
        0.2528774837460510d0, &
        0.6824976673552180d0, &
        0.6602062072459940d0, &
        0.4015093296585650d0]
    x3 = [0.3201877837239090d0, &
        0.0980256595288177d0, &
        0.6897988918691660d0, &
        0.1785484851694640d0, &
        0.0991062800273234d0, &
        0.1195800744029930d0, &
        0.8476199670433790d0, &
        0.8536320559829150d0, &
        0.6394323340044970d0, &
        0.8848230532535040d0, &
        0.9300526849294520d0, &
        0.6703901525053320d0, &
        0.7168448453351630d0, &
        0.9870657922150660d0, &
        0.2874068518452400d0]

    ! Equal Variance Assumption
    t = t_test(x1, x2, EQUAL_VARIANCE_ASSUMPTION)
    delta = ans1 - t%probability
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "T_TEST_TEST - EQUAL VARIANCE ASSUMPTION FAILED."
        print *, "Expected: ", ans1
        print *, "Computed: ", t%probability
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! Unequal Variance Assumption
    t = t_test(x1, x2, UNEQUAL_VARIANCE_ASSUMPTION)
    delta = ans2 - t%probability
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "T_TEST_TEST - UNEQUAL VARIANCE ASSUMPTION FAILED."
        print *, "Expected: ", ans2
        print *, "Computed: ", t%probability
        print *, "Difference (Expected - Computed): ", delta
    end if

    ! Paired Data Assumption
    t = t_test(x1, x3, PAIRED_DATA_SET_ASSUMPTION)
    delta = ans3 - t%probability
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "T_TEST_TEST - PAIRED DATA ASSUMPTION FAILED."
        print *, "Expected: ", ans3
        print *, "Computed: ", t%probability
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
function f_test_test() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n1 = 15
    integer(int32), parameter :: n2 = 15
    real(real64), parameter :: tol = 1.0d-6
    real(real64), parameter :: fans = 0.95278924436d0
    real(real64), parameter :: ans = 0.46459183414d0

    ! Local Variables
    real(real64) :: x1(n1), x2(n2), delta
    type(statistic) :: f

    ! Initialization
    rst = .true.
    x1 = [0.58387698128673400d0, &
        0.94872643185602900d0, &
        0.17008749217405500d0, &
        0.05400035877767780d0, &
        0.49007414183526400d0, &
        0.91498304265046300d0, &
        0.27730662176870100d0, &
        0.14672540316540200d0, &
        0.75310066420122200d0, &
        0.50398721165895900d0, &
        0.59636223974370000d0, &
        0.47819755122118100d0, &
        0.96709705143154200d0, &
        0.86367902851430800d0, &
        0.61084279368466600d0]
    x2 = [0.2186373486915680d0, &
        0.7881189936359480d0, &
        0.8980582651872070d0, &
        0.2036242143271850d0, &
        0.2080903307925950d0, &
        0.0143712341343581d0, &
        0.7474501270172240d0, &
        0.4900316414528420d0, &
        0.6221226099553770d0, &
        0.4632448747935360d0, &
        0.6138713552592840d0, &
        0.9901017538082670d0, &
        0.2822073800781740d0, &
        0.3928656727796740d0, &
        0.0401029216162131d0]

    ! Process
    f = f_test(x1, x2)
    delta = ans - f%probability
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "F_TEST_TEST FAILED - PROBABILITY TERM."
        print *, "Expected: ", ans
        print *, "Computed: ", f%probability
        print *, "Difference (Expected - Computed): ", delta
    end if
    
    delta = fans - f%value
    if (abs(delta) > tol) then
        rst = .false.
        print '(A)', "F_TEST_TEST FAILED."
        print *, "Expected: ", fans
        print *, "Computed: ", f%value
        print *, "Difference (Expected - Computed): ", delta
    end if
end function

! ------------------------------------------------------------------------------
end module
