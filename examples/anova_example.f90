! anova_example.f90

program main
    use iso_fortran_env
    use measurements_core
    implicit none

    ! Parameters
    integer(int32), parameter :: nops = 3
    integer(int32), parameter :: nparts = 5
    integer(int32), parameter :: ntrials = 3

    ! Local Variables
    real(real64) :: x(nparts, ntrials, nops)
    type(gage_anova_table) :: rst

    ! Operator 1 Results
    x(:,:,1) = reshape([&
        3.29d0, 2.44d0, 4.34d0, 3.47d0, 2.2d0, &
        3.41d0, 2.32d0, 4.17d0, 3.5d0, 2.08d0, &
        3.64d0, 2.42d0, 4.27d0, 3.64d0, 2.16d0], &
        [nparts, ntrials])

    ! Operator 2 Results
    x(:,:,2) = reshape([ &
        3.08d0, 2.53d0, 4.19d0, 3.01d0, 2.44d0, &
        3.25d0, 1.78d0, 3.94d0, 4.03d0, 1.8d0, &
        3.07d0, 2.32d0, 4.34d0, 3.2d0, 1.72d0], &
        [nparts, ntrials])

    ! Operator 3 Results
    x(:,:,3) = reshape([ &
        3.04d0, 1.62d0, 3.88d0, 3.14d0, 1.54d0, &
        2.89d0, 1.87d0, 4.09d0, 3.2d0, 1.93d0, &
        2.85d0, 2.04d0, 3.67d0, 3.11d0, 1.55d0], &
        [nparts, ntrials])

    ! Perform the ANOVA
    rst = gage_anova(x)

    ! Display the results
    print '(A)', "Operator Results:"
    print '(AI0)', achar(9) // "DOF: ", rst%operators%dof
    print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%operators%sum_of_squares
    print '(AF0.3)', achar(9) // "Mean of Squares (Variance): ", rst%operators%mean_of_squares
    print '(AF0.3)', achar(9) // "F Statistic: ", rst%operators%f_stat
    
    print '(A)', new_line('a') // "Part Results:"
    print '(AI0)', achar(9) // "DOF: ", rst%parts%dof
    print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%parts%sum_of_squares
    print '(AF0.3)', achar(9) // "Mean of Squares (Variance): ", rst%parts%mean_of_squares
    print '(AF0.3)', achar(9) // "F Statistic: ", rst%parts%f_stat

    print '(A)', new_line('a') // "Equipment Results:"
    print '(AI0)', achar(9) // "DOF: ", rst%equipment%dof
    print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%equipment%sum_of_squares
    print '(AF0.3)', achar(9) // "Mean of Squares (Variance): ", rst%equipment%mean_of_squares

    print '(A)', new_line('a') // "Operator-Part Interaction Results:"
    print '(AI0)', achar(9) // "DOF: ", rst%operator_by_part%dof
    print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%operator_by_part%sum_of_squares
    print '(AF0.3)', achar(9) // "Mean of Squares (Variance): ", rst%operator_by_part%mean_of_squares
    print '(AF0.3)', achar(9) // "F Statistic: ", rst%operator_by_part%f_stat

    print '(A)', new_line('a') // "Total Results:"
    print '(AI0)', achar(9) // "DOF: ", rst%total%dof
    print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%total%sum_of_squares
    print '(AF0.3)', achar(9) // "Mean of Squares (Variance): ", rst%total%mean_of_squares
    print '(AF0.3)', achar(9) // "Overall Mean: ", rst%total%mean
end program