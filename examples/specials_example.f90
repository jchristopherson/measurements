! specials_example.f90

program main
    use iso_fortran_env
    use measurements_core
    use fplot_core
    implicit none

    ! Variables
    integer(int32), parameter :: npts = 100
    real(real64) :: x(npts), f(npts), fa(npts)
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Evaluate the digamma function
    x = linspace(0.1d0, 5.0d0, npts)
    f = digamma(x)
    
    ! Evaluate an approximation to the function
    fa = log(x) - 1.0d0 / (2.0d0 * x)

    ! Plot the functions
    call plt%initialize()
    call plt%set_font_size(14)
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    lgnd => plt%get_legend()

    call xAxis%set_title("x")
    call yAxis%set_title("\psi(x)")

    call lgnd%set_is_visible(.true.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)

    call d1%define_data(x, f)
    call d1%set_name("Actual")
    call d1%set_line_width(2.0)
    call plt%push(d1)

    call d2%define_data(x, fa)
    call d2%set_name("Approximation")
    call d2%set_line_width(2.0)
    call d2%set_line_style(LINE_DASHED)
    call d2%set_line_color(CLR_RED)
    call plt%push(d2)

    call plt%draw()
end program