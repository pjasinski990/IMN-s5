set term pngcairo

set output "global.png"
set title "Globalna relaksacja"
set xlabel "nr iteracji"
set ylabel "S"
set logscale x
plot [:100000] 'global_sum_omega1.dat' u 1:2 with lines linestyle 1 t "omega = 0.6", \
     'global_sum_omega2.dat' u 1:2 with lines linestyle 2 t "omega = 1.0"

#set output "rk2_v(t).png"
#set title "Metoda RK2"
#set xlabel "t"
#set ylabel "v(t)"
#plot 'rk2_tol1.dat' u 1:4 with lines linestyle 1 t "TOL = 10e-2", \
#     'rk2_tol2.dat' u 1:4 with lines linestyle 2 t "TOL = 10e-5"
#
#set output "rk2_v(x).png"
#set title "Metoda RK2"
#set xlabel "x"
#set ylabel "v"
#plot 'rk2_tol1.dat' u 3:4 with lines linestyle 1 t "TOL = 10e-2", \
#     'rk2_tol2.dat' u 3:4 with lines linestyle 2 t "TOL = 10e-5"
#
#set output "rk2_dt(t).png"
#set title "Metoda RK2"
#set xlabel "t"
#set ylabel "dt(t)"
#plot 'rk2_tol1.dat' u 1:2 with lines linestyle 1 t "TOL = 10e-2", \
#     'rk2_tol2.dat' u 1:2 with lines linestyle 2 t "TOL = 10e-5"
#
##  Trapezy
#
#set output "trap_x(t).png"
#set title "Metoda trapezow"
#set xlabel "t"
#set ylabel "x(t)"
#plot 'trap_tol1.dat' u 1:3 with lines linestyle 1 t "TOL = 10e-2", \
#     'trap_tol2.dat' u 1:3 with lines linestyle 2 t "TOL = 10e-5"
#
#set output "trap_v(t).png"
#set title "Metoda trapezow"
#set xlabel "t"
#set ylabel "v(t)"
#plot 'trap_tol1.dat' u 1:4 with lines linestyle 1 t "TOL = 10e-2", \
#     'trap_tol2.dat' u 1:4 with lines linestyle 2 t "TOL = 10e-5"
#
#set output "trap_v(x).png"
#set title "Metoda trapezow"
#set xlabel "x"
#set ylabel "v"
#plot 'trap_tol1.dat' u 3:4 with lines linestyle 1 t "TOL = 10e-2", \
#     'trap_tol2.dat' u 3:4 with lines linestyle 2 t "TOL = 10e-5"
#
#set output "trap_dt(t).png"
#set title "Metoda trapezow"
#set xlabel "t"
#set ylabel "dt(t)"
#plot 'trap_tol1.dat' u 1:2 with lines linestyle 1 t "TOL = 10e-2", \
#     'trap_tol2.dat' u 1:2 with lines linestyle 2 t "TOL = 10e-5"
