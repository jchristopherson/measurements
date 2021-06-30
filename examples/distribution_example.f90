! distribution_example.f90

program main
    use iso_fortran_env
    use measurements_core
    use fplot_core
    implicit none

    ! Local Variables
    integer(int32), parameter :: npts = 100
    real(real64) :: x(npts), pdf(npts), cdf(npts)
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(legend), pointer :: lgnd
    type(normal_distribution) :: nd
    type(log_normal_distribution) :: lnd
    type(beta_distribution) :: bd
    type(f_distribution) :: fd

    ! Initialize the plot & set up font properties to improve readability
    call plt%initialize()
    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_horizontal_position(LEGEND_LEFT)
    call plt%set_use_y2_axis(.true.)

    ! Look at a Normal Distribution
    x = linspace(-5.0d0, 5.0d0, npts)
    call nd%set_mu(0.0d0)
    call nd%set_sigma(0.5d0)
    pdf = nd%pdf(x)
    cdf = nd%cdf(x)

    ! Plot the data
    call plt%set_title("Normal Distribution")

    call pd1%define_data(x, pdf)
    call pd1%set_name("PDF")
    call pd1%set_line_width(2.0)
    call plt%push(pd1)

    call pd2%define_data(x, cdf)
    call pd2%set_name("CDF (Y2)")
    call pd2%set_line_width(2.0)
    call pd2%set_line_style(LINE_DASHED)
    call pd2%set_draw_against_y2(.true.)
    call plt%push(pd2)

    call plt%draw()

    ! Look at Student's T Distribution
    pdf = t_distribution_pdf(10.0d0, x)
    cdf = t_distribution_cdf(10.0d0, x)

    ! Plot the data
    call plt%set_title("Student's T Distribution")
    call plt%clear_all()

    call pd1%define_data(x, pdf)
    call plt%push(pd1)

    call pd2%define_data(x, cdf)
    call plt%push(pd2)

    call plt%draw()

    ! Look at an F-Distribution
    x = linspace(0.0d0, 5.0d0, npts)
    call fd%set_d1(5.0d0)
    call fd%set_d2(4.5d0)
    pdf = fd%pdf(x)
    cdf = fd%cdf(x)

    ! Display Model Parameters
    print '(A)', new_line('a') // "F-Distribution:"
    print '(AF0.4)', achar(9) // "Mean = ", fd%mean()
    print '(AF0.4)', achar(9) // "Mode = ", fd%mode()
    print '(AF0.4)', achar(9) // "Variance = ", fd%variance()

    ! Plot the data
    call plt%set_title("F-Distribution")
    call lgnd%set_horizontal_position(LEGEND_RIGHT)
    call lgnd%set_vertical_position(LEGEND_CENTER)
    call plt%clear_all()

    call pd1%define_data(x, pdf)
    call plt%push(pd1)

    call pd2%define_data(x, cdf)
    call plt%push(pd2)

    call plt%draw()

    ! Look at a Beta Distribution
    x = linspace(0.0d0, 1.0d0, npts)
    call bd%set_alpha(2.0d0)
    call bd%set_beta(5.0d0)
    pdf = bd%pdf(x)
    cdf = bd%cdf(x)

    ! Display Model Parameters
    print '(A)', new_line('a') // "Beta Distribution:"
    print '(AF0.4)', achar(9) // "Mean = ", bd%mean()
    print '(AF0.4)', achar(9) // "Geometric Mean = ", bd%geometric_mean()
    print '(AF0.4)', achar(9) // "Median = ", bd%median()
    print '(AF0.4)', achar(9) // "Mode = ", bd%mode()
    print '(AF0.4)', achar(9) // "Variance = ", bd%variance()

    ! Plot the data
    call plt%set_title("Beta Distribution")
    call lgnd%set_horizontal_position(LEGEND_RIGHT)
    call lgnd%set_vertical_position(LEGEND_CENTER)
    call plt%clear_all()

    call pd1%define_data(x, pdf)
    call plt%push(pd1)

    call pd2%define_data(x, cdf)
    call plt%push(pd2)

    call plt%draw()

    ! Look at the log normal distribution
    x = linspace(0.0d0, 3.0d0, npts)
    call lnd%set_mu(0.0d0)
    call lnd%set_sigma(0.5d0)
    pdf = lnd%pdf(x)
    cdf = lnd%cdf(x)

    ! Plot the data
    call plt%set_title("Log Normal Distribution")
    call plt%clear_all()

    call pd1%define_data(x, pdf)
    call plt%push(pd1)

    call pd2%define_data(x, cdf)
    call plt%push(pd2)

    call plt%draw()
end program