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

    ! Initialize the plot & set up font properties to improve readability
    call plt%initialize()
    call plt%set_font_size(14)
    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_horizontal_position(LEGEND_LEFT)
    call plt%set_use_y2_axis(.true.)

    ! Look at a Normal Distribution
    x = linspace(-5.0d0, 5.0d0, npts)
    call nd%set_mean(0.0d0)
    call nd%set_standard_deviation(0.5d0)
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

    ! Look at an F Distribution
    x = linspace(0.0d0, 5.0d0, npts)
    pdf = f_distribution_pdf(5.0d0, 2.0d0, x)
    cdf = f_distribution_cdf(5.0d0, 2.0d0, x)

    ! Plot the data
    call plt%set_title("F Distribution")
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
    pdf = beta_distribution_pdf(2.0d0, 5.0d0, x)
    cdf = beta_distribution_cdf(2.0d0, 5.0d0, x)

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
    pdf = log_normal_distribution_pdf(0.0d0, 0.5d0, x)
    cdf = log_normal_distribution_cdf(0.0d0, 0.5d0, x)

    ! Plot the data
    call plt%set_title("Log Normal Distribution")
    call plt%clear_all()

    call pd1%define_data(x, pdf)
    call plt%push(pd1)

    call pd2%define_data(x, cdf)
    call plt%push(pd2)

    call plt%draw()
end program