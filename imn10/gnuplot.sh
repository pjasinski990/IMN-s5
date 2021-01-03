set term pngcairo
set tics out
set tics nomirror

set output "E_0.0_0.0.png"
plot "E_0.0_0.0.dat" w l, \
    "E_0.0_0.1.dat" w l, \
    "E_0.0_1.0.dat" w l

set size ratio -1
set xlabel "t"
set ylabel "x"

set view map

set output "v_0.0_0.0.png"
set title "v, alpha = 0.0, beta = 0.0"
splot 'v_0.0_0.0.dat' u 1:2:3 with pm3d

set output "v_0.0_0.1.png"
set title "v, alpha = 0.0, beta = 0.1"
splot 'v_0.0_0.1.dat' u 1:2:3 with pm3d

set output "v_0.0_1.0.png"
set title "v, alpha = 0.0, beta = 1.0"
splot 'v_0.0_1.0.dat' u 1:2:3 with pm3d

set output "v_1.0_1.0.png"
set title "v, alpha = 1.0, beta = 1.0"
splot 'v_1.0_1.0.dat' u 1:2:3 with pm3d
