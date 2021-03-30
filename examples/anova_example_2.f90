! anova_example_2.f90
! REF: https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_HypothesisTesting-ANOVA/BS704_HypothesisTesting-Anova3.html#:~:text=The%20ANOVA%20table%20breaks%20down,Source%20of%20Variation

program main
    use iso_fortran_env
    use measurements_core
    implicit none

    ! Parameters
    integer(int32), parameter :: nsamples = 5
    integer(int32), parameter :: ncategories = 4

    ! Local Variables
    real(real64) :: x(nsamples, ncategories)
    type(anova_table) :: rst

    ! Data by category
    x = reshape([ &
        8.0d0, 9.0d0, 6.0d0, 7.0d0, 3.0d0, &
        2.0d0, 4.0d0, 3.0d0, 5.0d0, 1.0d0, &
        3.0d0, 5.0d0, 4.0d0, 2.0d0, 3.0d0, &
        2.0d0, 2.0d0, -1.0d0, 0.0d0, 3.0d0], [nsamples, ncategories])

    ! Perform the ANOVA
    rst = anova(x)

    ! Display the results
    print '(A)', "Between Category Results:"
    print '(AI0)', achar(9) // "DOF: ", rst%between%dof
    print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%between%sum_of_squares
    print '(AF0.3)', achar(9) // "Mean of Squares: ", rst%between%mean_of_squares
    print '(AF0.3)', achar(9) // "F Statistic: ", rst%between%f_stat
    print '(AF0.5)', achar(9) // "Probability: ", rst%between%probability

    print '(A)', new_line('a') // "Residual Results:"
    print '(AI0)', achar(9) // "DOF: ", rst%residual%dof
    print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%residual%sum_of_squares
    print '(AF0.3)', achar(9) // "Mean of Squares: ", rst%residual%mean_of_squares

    print '(A)', new_line('a') // "Total Results:"
    print '(AI0)', achar(9) // "DOF: ", rst%total%dof
    print '(AF0.3)', achar(9) // "Sum of Squares: ", rst%total%sum_of_squares
    print '(AF0.3)', achar(9) // "Mean of Squares: ", rst%total%mean_of_squares
end program
