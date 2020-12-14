set term pngcairo

set output "xsr_ct.png"
set xlabel "t"
set ylabel "xsr, c(t)"
set yrange[0:4]

plot 'no_diff.dat' u 1:2 with lines t "c (D=0)", \
   'no_diff.dat' u 1:3 with lines t "xsr (D=0)", \
   'with_diff.dat' u 1:2 with lines t "c (D=0.1)", \
   'with_diff.dat' u 1:3 with lines t "xsr (D=0.1)"

unset yrange
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
