set term pngcairo

set output "xsr_ct_d0.png"
set xlabel "t"
set ylabel "xsr, c(t)"
set yrange[0:4]

set title "D=0.0"
plot 'no_diff.dat' u 1:2 with lines t "c (D=0)", \
   'no_diff.dat' u 1:3 with lines t "xsr (D=0)"

set title "D=0.1"
set output "xsr_ct_c1.png"
plot 'with_diff.dat' u 1:2 with lines t "c (D=0.1)", \
   'with_diff.dat' u 1:3 with lines t "xsr (D=0.1)"

unset yrange

set title "vx"
set output "vx.png"
set view map
splot "vx.dat" u 1:2:3 w pm3d

set title "vy"
set output "vy.png"
set view map
splot "vy.dat" u 1:2:3 w pm3d

set title "D=0.0"
set output "dist_no_diff_1.png"
set view map
splot "dist_no_diff_1.dat" u 1:2:3 w pm3d

set output "dist_no_diff_2.png"
set view map
splot "dist_no_diff_2.dat" u 1:2:3 w pm3d

set output "dist_no_diff_3.png"
set view map
splot "dist_no_diff_3.dat" u 1:2:3 w pm3d

set output "dist_no_diff_4.png"
set view map
splot "dist_no_diff_4.dat" u 1:2:3 w pm3d

set output "dist_no_diff_5.png"
set view map
splot "dist_no_diff_5.dat" u 1:2:3 w pm3d


set title "D=0.1"
set output "dist_with_diff_1.png"
set view map
splot "dist_with_diff_1.dat" u 1:2:3 w pm3d

set output "dist_with_diff_2.png"
set view map
splot "dist_with_diff_2.dat" u 1:2:3 w pm3d

set output "dist_with_diff_3.png"
set view map
splot "dist_with_diff_3.dat" u 1:2:3 w pm3d

set output "dist_with_diff_4.png"
set view map
splot "dist_with_diff_4.dat" u 1:2:3 w pm3d

set output "dist_with_diff_5.png"
set view map
splot "dist_with_diff_5.dat" u 1:2:3 w pm3d
