set term pngcairo

set xlabel "x"
set ylabel "y"
set tics out
set tics nomirror

set palette defined (-5 0 0 1, 0 1 1 1, 5 1 0 0)

set xrange[0:5]
set yrange[0:5]
set output "n50.png"
set title "n = 50"
plot 'n50.dat' u 1:2:3 with image

set xrange[0:10]
set yrange[0:10]
set output "n100.png"
set title "n = 100"
plot 'n100.dat' u 1:2:3 with image

set xrange[0:20]
set yrange[0:20]
set output "n200.png"
set title "n = 200"
plot 'n200.dat' u 1:2:3 with image

set xrange[0:10]
set yrange[0:10]
set cbrange[-0.8:0.8]
set output "eps1.png"
set title "eps1 = 1, eps2 = 1"
plot 'eps1.dat' u 1:2:3 with image

set output "eps2.png"
set title "eps1 = 1, eps2 = 2"
plot 'eps2.dat' u 1:2:3 with image

set output "eps10.png"
set title "eps1 = 1, eps2 = 10"
plot 'eps10.dat' u 1:2:3 with image
