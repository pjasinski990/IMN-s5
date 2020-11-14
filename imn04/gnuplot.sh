set term pngcairo

set xlabel "nr iteracji"
set ylabel "S"
set logscale x

set output "global_sum.png"
set title "Globalna relaksacja"
plot [:100000] 'global_sum_omega1.dat' u 1:2 with lines linestyle 1 t "omega = 0.6", \
     'global_sum_omega2.dat' u 1:2 with lines linestyle 2 t "omega = 1.0"

set output "local_sum.png"
set title "Lokalna relaksacja"
plot [:100000] 'local_sum_omega1.dat' u 1:2 with lines linestyle 1 t "omega = 1.0", \
     'local_sum_omega2.dat' u 1:2 with lines linestyle 2 t "omega = 1.4", \
     'local_sum_omega3.dat' u 1:2 with lines linestyle 3 t "omega = 1.8", \
     'local_sum_omega4.dat' u 1:2 with lines linestyle 4 t "omega = 1.9"


unset logscale x
set xlabel "x"
set ylabel "y"

set output "solution1.png"
set title "Potencjal V (omega = 0.6)"
plot 'global_solution_omega1.dat' u 1:2:3 with points pt 5 lc palette

set output "sigma1.png"
set title "Blad rozwiazania relaksacja globalna, omega = 0.6"
plot 'global_sigma_omega1.dat' u 1:2:3 with points pt 5 lc palette


set output "solution2.png"
set title "Potencjal V (omega = 1.0)"
plot 'global_solution_omega2.dat' u 1:2:3 with points pt 5 lc palette

set output "sigma2.png"
set title "Blad rozwiazania relaksacja globalna, omega = 1.0"
plot 'global_sigma_omega2.dat' u 1:2:3 with points pt 5 lc palette
