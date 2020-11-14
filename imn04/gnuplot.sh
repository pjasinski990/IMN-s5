set term pngcairo

set xlabel "nr iteracji"
set ylabel "S"
set logscale x

set output "global.png"
set title "Globalna relaksacja"
plot [:100000] 'global_sum_omega1.dat' u 1:2 with lines linestyle 1 t "omega = 0.6", \
     'global_sum_omega2.dat' u 1:2 with lines linestyle 2 t "omega = 1.0"

set output "local.png"
set title "Lokalna relaksacja"
plot [:100000] 'local_sum_omega1.dat' u 1:2 with lines linestyle 1 t "omega = 1.0", \
     'local_sum_omega2.dat' u 1:2 with lines linestyle 2 t "omega = 1.4", \
     'local_sum_omega3.dat' u 1:2 with lines linestyle 3 t "omega = 1.8", \
     'local_sum_omega4.dat' u 1:2 with lines linestyle 4 t "omega = 1.9"
