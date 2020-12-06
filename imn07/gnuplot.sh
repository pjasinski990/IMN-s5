set term pngcairo

set size ratio -1
set xlabel 'x'
set ylabel 'y'
set tics out
set tics nomirror

#contours
set contour 
unset surface
set view map
unset key

# psi contours 
set output "psi_contour_neg1000.png"
set title "psi, Q = -1000"
set cbr[-55:-50]
set cntrparam levels increment -55, 0.2, -50
splot 'Q=neg1000.dat' u 1:2:3:3 with lines palette

set output "psi_contour_neg4000.png"
set title "psi, Q = -4000"
set cbr[-218:-202]
set cntrparam levels increment -218, 1, -202
splot 'Q=neg4000.dat' u 1:2:3:3 with lines palette

set output "psi_contour_pos4000.png"
set title "psi, Q = 4000"
set cbr[202:218]
set cntrparam levels increment 202, 1, 218
splot 'Q=pos4000.dat' u 1:2:3:3 with lines palette

# zeta contours 
set output "zeta_contour_neg1000.png"
set title "zeta, Q = -1000"
set cbr[-200:300]
set cntrparam levels increment -200, 10, 300
splot 'Q=neg1000.dat' u 1:2:4:4 with lines palette

set output "zeta_contour_neg4000.png"
set title "zeta, Q = -4000"
set cbr[-500:1000]
set cntrparam levels increment -500, 50, 1000
splot 'Q=neg4000.dat' u 1:2:4:4 with lines palette


# maps
unset contour 
set surface
unset cbr

# neg1000
set output "u_map_neg1000.png" 
set title "u(x, y), Q = -1000"
set view map
splot "Q=neg1000.dat" u 1:2:5 w pm3d

set output "v_map_neg1000.png" 
set title "v(x, y), Q = -1000"
set view map
splot "Q=neg1000.dat" u 1:2:6 w pm3d

# neg 4000
set output "u_map_neg4000.png" 
set title "u(x, y), Q = -4000"
set view map
splot "Q=neg4000.dat" u 1:2:5 w pm3d

set output "v_map_neg4000.png" 
set title "v(x, y), Q = -4000"
set view map
splot "Q=neg4000.dat" u 1:2:6 w pm3d

