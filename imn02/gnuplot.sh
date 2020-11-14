set style line 1 \
    linecolor rgb "royalblue" \
    linetype 1 linewidth 2
set style line 2 \
    linecolor rgb "green" \
    linetype 1 linewidth 2

set term png

set output "picard.png"
set title "Metoda Picarda"
set xlabel "t"
set ylabel "u(t), z(t)"
plot 'picard.dat' u 1:2 with lines linestyle 1 t "u(t)", \
     'picard.dat' u 1:3 with lines linestyle 2 t "z(t)"

set output "newton.png"
set title "Iteracja Newtona"
set xlabel "t"
set ylabel "u(t), z(t)"
plot 'newton.dat' u 1:2 with lines linestyle 1 t "u(t)", \
     'newton.dat' u 1:3 with lines linestyle 2 t "z(t)"

set output "rk2.png"
set title "Niejawna metoda RK2"
set xlabel "t"
set ylabel "u(t), z(t)"
plot 'rk2.dat' u 1:2 with lines linestyle 1 t "u(t)", \
     'rk2.dat' u 1:3 with lines linestyle 2 t "z(t)"
