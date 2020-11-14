set style line 1 \
    linecolor rgb "royalblue" \
    linetype 1 \
    linewidth 2

set style line 2 \
    linecolor rgb '#5e9c36' \
    linetype 1 \
    linewidth 2

set style line 3 \
    linecolor rgb "gray30" \
    linetype 1 \
    linewidth 2

set term pngcairo

# RK2

set output "rk2_x(t).png"
set title "Metoda RK2"
set xlabel "t"
set ylabel "x(t)"
plot 'rk2_tol1.dat' u 1:3 with lines linestyle 1 t "TOL = 10e-2", \
     'rk2_tol2.dat' u 1:3 with lines linestyle 2 t "TOL = 10e-5"

set output "rk2_v(t).png"
set title "Metoda RK2"
set xlabel "t"
set ylabel "v(t)"
plot 'rk2_tol1.dat' u 1:4 with lines linestyle 1 t "TOL = 10e-2", \
     'rk2_tol2.dat' u 1:4 with lines linestyle 2 t "TOL = 10e-5"

set output "rk2_v(x).png"
set title "Metoda RK2"
set xlabel "x"
set ylabel "v"
plot 'rk2_tol1.dat' u 3:4 with lines linestyle 1 t "TOL = 10e-2", \
     'rk2_tol2.dat' u 3:4 with lines linestyle 2 t "TOL = 10e-5"

set output "rk2_dt(t).png"
set title "Metoda RK2"
set xlabel "t"
set ylabel "dt(t)"
plot 'rk2_tol1.dat' u 1:2 with lines linestyle 1 t "TOL = 10e-2", \
     'rk2_tol2.dat' u 1:2 with lines linestyle 2 t "TOL = 10e-5"

#  Trapezy

set output "trap_x(t).png"
set title "Metoda trapezow"
set xlabel "t"
set ylabel "x(t)"
plot 'trap_tol1.dat' u 1:3 with lines linestyle 1 t "TOL = 10e-2", \
     'trap_tol2.dat' u 1:3 with lines linestyle 2 t "TOL = 10e-5"

set output "trap_v(t).png"
set title "Metoda trapezow"
set xlabel "t"
set ylabel "v(t)"
plot 'trap_tol1.dat' u 1:4 with lines linestyle 1 t "TOL = 10e-2", \
     'trap_tol2.dat' u 1:4 with lines linestyle 2 t "TOL = 10e-5"

set output "trap_v(x).png"
set title "Metoda trapezow"
set xlabel "x"
set ylabel "v"
plot 'trap_tol1.dat' u 3:4 with lines linestyle 1 t "TOL = 10e-2", \
     'trap_tol2.dat' u 3:4 with lines linestyle 2 t "TOL = 10e-5"

set output "trap_dt(t).png"
set title "Metoda trapezow"
set xlabel "t"
set ylabel "dt(t)"
plot 'trap_tol1.dat' u 1:2 with lines linestyle 1 t "TOL = 10e-2", \
     'trap_tol2.dat' u 1:2 with lines linestyle 2 t "TOL = 10e-5"
