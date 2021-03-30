! measurements_anova.f90

! Info & References:
! - https://www.spcforexcel.com/knowledge/measurement-systems-analysis/anova-gage-rr-part-1
! - https://simongrund1.github.io/posts/anova-with-multiply-imputed-data-sets/

submodule (measurements_core) measurements_anova
    use ieee_arithmetic

contains
! ------------------------------------------------------------------------------
    ! REF: https://www.spcforexcel.com/knowledge/measurement-systems-analysis/anova-gage-rr-part-1
    module function gage_anova(x, err) result(rst)
        ! Arguments
        real(real64), intent(in), target, contiguous, dimension(:,:,:) :: x
        class(errors), intent(inout), optional, target :: err
        type(gage_anova_table) :: rst

        ! Local Variables
        integer(int32) :: i, j, k, nops, nresults, ntrials, n, nr, kr, flag
        real(real64) :: avg, nan
        real(real64), pointer, dimension(:) :: xptr
        real(real64), allocatable, dimension(:) :: xsub
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg
        
        ! Set up error handling
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Initialization
        nresults = size(x, 1)
        ntrials = size(x, 2)
        nops = size(x, 3)
        n = nresults * ntrials * nops
        nan = ieee_value(nan, ieee_quiet_nan)

        ! Input Checking
        if (nresults < 2) then
            write(errmsg, '(AI0A)') "Only ", nresults, &
                " part results provided.  A minimum of 2 results " // &
                "sets are required."
            call errmgr%report_error("gage_anova", trim(errmsg), &
                M_INSUFFICIENT_DATA_ERROR)
            return
        end if

        if (ntrials < 2) then
            write(errmsg, '(AI0A)') "Only ", ntrials, &
                " trials for each part were provided.  A minimum of 2 " // &
                "results sets are required."
            call errmgr%report_error("gage_anova", trim(errmsg), &
                M_INSUFFICIENT_DATA_ERROR)
            return
        end if

        ! Memory Allocation
        allocate(rst%operator(nops), stat = flag)
        if (flag == 0) allocate(rst%part(nresults), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("gage_anova", &
                "Insufficient memory available.", M_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Compute the overall mean of the input data
        xptr(1:n) => x
        rst%total%mean = mean(xptr)

        ! Compute the total sum of squares (variance of xptr * (n - 1))
        rst%total%mean_of_squares = variance(xptr)
        rst%total%sum_of_squares = rst%total%mean_of_squares * (n - 1.0d0)
        rst%total%dof = n - 1

        ! Compute operator information
        nr = nresults * ntrials
        rst%operators%sum_of_squares = 0.0d0
        rst%operators%dof = nops - 1
        do i = 1, nops
            ! Compute the results mean for operator i
            xptr(1:nr) => x(:,:,i)
            rst%operator(i)%mean = mean(xptr)
            rst%operator(i)%dof = nr - 1

            ! Compute the square of the difference between the operator mean 
            ! and the overall mean
            rst%operator(i)%sum_of_squares = nr * &
                (rst%operator(i)%mean - rst%total%mean)**2
            rst%operator(i)%mean_of_squares = &
                rst%operator(i)%sum_of_squares / rst%operator(i)%dof

            ! Add to the existing sum of squares
            rst%operators%sum_of_squares = &
                rst%operators%sum_of_squares + rst%operator(i)%sum_of_squares
        end do
        if (nops == 1) then
            rst%operators%mean_of_squares = 0.0d0
        else
            rst%operators%mean_of_squares = &
                rst%operators%sum_of_squares / rst%operators%dof
        end if
        rst%operators%mean = rst%total%mean

        ! Compute the part information
        kr = nops * ntrials
        rst%parts%sum_of_squares = 0.0d0
        rst%parts%dof = nresults - 1
        do i = 1, nresults
            ! Compute the part average across all operators
            xsub = reshape(x(i,:,:), [kr]) ! Organize for calculation purposes

            ! Compute the mean for part 1
            rst%part(i)%mean = mean(xsub)
            rst%part(i)%dof = kr - 1

            ! Compute the square of the difference between the part mean
            ! and the overall mean
            rst%part(i)%sum_of_squares = kr * &
                (rst%part(i)%mean - rst%total%mean)**2
            rst%part(i)%mean_of_squares = &
                rst%part(i)%sum_of_squares / rst%part(i)%dof

            ! Add to the existing sum of squares
            rst%parts%sum_of_squares = &
                rst%parts%sum_of_squares + rst%part(i)%sum_of_squares
        end do
        rst%parts%mean_of_squares = rst%parts%sum_of_squares / rst%parts%dof
        rst%parts%mean = rst%total%mean

        ! Compute the measurement equipment information
        rst%equipment%sum_of_squares = 0.0d0
        rst%equipment%dof = nresults * nops * (ntrials - 1)
        do i = 1, nops
            do j = 1, nresults
                ! Compute the mean of the results from a given part for a
                ! given operator
                avg = mean(x(j,:,i))

                ! Compute the sum of the squares of the difference between
                ! the measured value and the above calculated mean
                do k = 1, ntrials
                    rst%equipment%sum_of_squares = &
                        rst%equipment%sum_of_squares + (x(j,k,i) - avg)**2
                end do
            end do
        end do
        rst%equipment%mean_of_squares = &
            rst%equipment%sum_of_squares / rst%equipment%dof
        rst%equipment%mean = rst%total%mean

        ! Compute the operator-by-part interaction
        rst%operator_by_part%dof = (nops - 1) * (nresults - 1)
        rst%operator_by_part%sum_of_squares = rst%total%sum_of_squares - &
            (rst%operators%sum_of_squares + rst%parts%sum_of_squares + &
            rst%equipment%sum_of_squares)
        rst%operator_by_part%mean_of_squares = &
            rst%operator_by_part%sum_of_squares / rst%operator_by_part%dof

        ! Compute the F statistics
        rst%total%f_stat = nan
        rst%equipment%f_stat = nan
        rst%operator_by_part%f_stat = &
            rst%operator_by_part%mean_of_squares / rst%equipment%mean_of_squares
        rst%parts%f_stat = &
            rst%parts%mean_of_squares / rst%operator_by_part%mean_of_squares
        rst%operators%f_stat = &
            rst%operators%mean_of_squares / rst%operator_by_part%mean_of_squares

        ! Compute the probability terms
        rst%total%probability = nan
        rst%equipment%probability = nan
        
        rst%operator_by_part%probability = 1.0d0 - ftest_probability( &
            rst%operator_by_part%f_stat, &
            rst%operator_by_part%dof, &
            rst%equipment%dof)
        
        rst%parts%probability = 1.0d0 - ftest_probability(rst%parts%f_stat, &
            rst%parts%dof, rst%operator_by_part%dof)

        rst%operators%probability = 1.0d0 - ftest_probability( &
            rst%operators%f_stat, rst%operators%dof, rst%operator_by_part%dof)
    end function

! ------------------------------------------------------------------------------
    ! https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_HypothesisTesting-ANOVA/BS704_HypothesisTesting-Anova3.html#:~:text=The%20ANOVA%20table%20breaks%20down,Source%20of%20Variation
    module function anova(x, err) result(rst)
        ! Arguments
        real(real64), intent(in), target, contiguous, dimension(:,:) :: x
        class(errors), intent(inout), optional, target :: err
        type(anova_table) :: rst

        ! Local Variables
        integer(int32) :: i, j, npoints, nsets, n, flag
        real(real64), pointer, dimension(:) :: xptr
        real(real64), allocatable, dimension(:) :: avgs
        real(real64) :: nan
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg
        
        ! Set up error handling
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Initialization
        npoints = size(x, 1)
        nsets = size(x, 2)
        n = npoints * nsets
        nan = ieee_value(nan, ieee_quiet_nan)

        ! Input Checking
        if (nsets < 2) then
            write(errmsg, '(AI0A)') "Only ", nsets, " data sets " // &
                "were provided.  A minimum of 2 data sets are required."
            call errmgr%report_error("anova", trim(errmsg), &
                M_INSUFFICIENT_DATA_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(avgs(nsets), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("anova", &
                "Insufficient memory available.", &
                M_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Compute the overall mean
        xptr(1:n) => x
        rst%total%mean = mean(xptr)

        ! Compute the overall variance
        rst%total%mean_of_squares = variance(xptr)
        rst%total%sum_of_squares = rst%total%mean_of_squares * (n - 1.0d0)
        rst%total%dof = n - 1

        ! Look at each data set
        rst%between%sum_of_squares = 0.0d0
        rst%between%dof = nsets - 1
        do i = 1, nsets
            ! Compute the mean for the data set
            avgs(i) = mean(x(:,i))

            ! Compute thet sum of squares contribution
            rst%between%sum_of_squares = rst%between%sum_of_squares + &
                npoints * (avgs(i) - rst%total%mean)**2
        end do
        rst%between%mean_of_squares = &
            rst%between%sum_of_squares / rst%between%dof

        ! Compute the residuals
        rst%residual%sum_of_squares = 0.0d0
        rst%residual%dof = n - nsets
        do j = 1, nsets
            do i = 1, npoints
                rst%residual%sum_of_squares = rst%residual%sum_of_squares + &
                    (x(i,j) - avgs(j))**2
            end do
        end do
        rst%residual%mean_of_squares = rst%residual%sum_of_squares / &
            rst%residual%dof

        ! Compute the F statistic
        rst%between%f_stat = rst%between%mean_of_squares / &
            rst%residual%mean_of_squares
        rst%residual%f_stat = 0.0d0
        rst%total%f_stat = 0.0d0

        ! Compute the probability term
        rst%between%probability = 1.0d0 - ftest_probability( &
            rst%between%f_stat, rst%between%dof, rst%residual%dof)
        rst%residual%probability = nan
        rst%total%probability = nan
    end function

! ------------------------------------------------------------------------------
end submodule
