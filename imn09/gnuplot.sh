set term pngcairo
set size ratio -1
set tics out
set tics nomirror

set xlabel "x"
set ylabel "y"

set view map

set output "T_it100.png"
set title "T(x, y), it=100"
splot 'T_it100.dat' u 1:2:3 with pm3d

set output "T_diff_it100.png"
set title "diff T(x, y), it=100"
splot 'T_diff_it100.dat' u 1:2:3 with pm3d

set output "T_it200.png"
set title "T(x, y), it=200"
splot 'T_it200.dat' u 1:2:3 with pm3d

set output "T_diff_it200.png"
set title "diff T(x, y), it=200"
splot 'T_diff_it200.dat' u 1:2:3 with pm3d

set output "T_it500.png"
set title "T(x, y), it=500"
splot 'T_it500.dat' u 1:2:3 with pm3d

set output "T_diff_it500.png"
set title "diff T(x, y), it=500"
splot 'T_diff_it500.dat' u 1:2:3 with pm3d

set output "T_it1000.png"
set title "T(x, y), it=1000"
splot 'T_it1000.dat' u 1:2:3 with pm3d

set output "T_diff_it1000.png"
set title "diff T(x, y), it=1000"
splot 'T_diff_it1000.dat' u 1:2:3 with pm3d

set output "T_it2000.png"
set title "T(x, y), it=2000"
splot 'T_it2000.dat' u 1:2:3 with pm3d

set output "T_diff_it2000.png"
set title "diff T(x, y), it=2000"
splot 'T_diff_it2000.dat' u 1:2:3 with pm3d
