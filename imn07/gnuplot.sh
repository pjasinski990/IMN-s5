set term pngcairo

set xlabel "x"
set ylabel "y"
set tics out
set tics nomirror

#set palette defined (-5 0 0 1, 0 1 1 1, 5 1 0 0)

set xrange[0:200]
set yrange[0:90]
set output "psi.png"
set title "psi"
plot 'psi.dat' u 1:2:3 with image



set xrange[0:200]
set yrange[0:90]
set output "zeta.png"
set title "zeta"
plot 'zeta.dat' u 1:2:3 with image
