set term pngcairo

set xlabel "x"
set ylabel "y"
set xrange[0:25]
set yrange[0:25]
set tics out
set tics nomirror

set output "k16.png"
set title "k = 16"
plot 'k16.dat' u 1:2:3 with image

set output "k8.png"
set title "k = 8"
plot 'k8.dat' u 1:2:3 with image

set output "k4.png"
set title "k = 4"
plot 'k4.dat' u 1:2:3 with image

set output "k2.png"
set title "k = 2"
plot 'k2.dat' u 1:2:3 with image

set output "k1.png"
set title "k = 1"
plot 'k1.dat' u 1:2:3 with image
