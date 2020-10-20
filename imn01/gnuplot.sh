set style line 1 \
    linecolor rgb "red" \
    linetype 1 linewidth 2
set style line 2 \
    linecolor rgb "green" \
    linetype 1 linewidth 2
set style line 3 \
    linecolor rgb "royalblue" \
    linetype 1 linewidth 2
set style line 3 \
    linecolor rgb "violet" \
    linetype 1 linewidth 2

set term png

set output "1a.png"
plot '1a.dat' using 1:2 index 0 with linespoints linestyle 1 t "dt = 0.01", \
    ''      using 1:2 index 1 with linespoints linestyle 2 t "dt = 0.1", \
    ''      using 1:2 index 2 with linespoints linestyle 3 t "dt = 1" 

set output "1a_err.png"
plot '1a.dat' using 1:3 index 0 with linespoints linestyle 1 t "dt = 0.01", \
    ''      using 1:3 index 1 with linespoints linestyle 2 t "dt = 0.1", \
    ''      using 1:3 index 2 with linespoints linestyle 3 t "dt = 1" 

set output "1b.png"
plot '1b.dat' using 1:2 index 0 with linespoints linestyle 1 t "dt = 0.01", \
    ''      using 1:2 index 1 with linespoints linestyle 2 t "dt = 0.1", \
    ''      using 1:2 index 2 with linespoints linestyle 3 t "dt = 1" 

set output "1b_err.png"
plot '1b.dat' using 1:3 index 0 with linespoints linestyle 1 t "dt = 0.01", \
    ''      using 1:3 index 1 with linespoints linestyle 2 t "dt = 0.1", \
    ''      using 1:3 index 2 with linespoints linestyle 3 t "dt = 1" 

set output "1c.png"
plot '1c.dat' using 1:2 index 0 with linespoints linestyle 1 t "dt = 0.01", \
    ''      using 1:2 index 1 with linespoints linestyle 2 t "dt = 0.1", \
    ''      using 1:2 index 2 with linespoints linestyle 3 t "dt = 1" 

set output "1c_err.png"
plot '1c.dat' using 1:3 index 0 with linespoints linestyle 1 t "dt = 0.01", \
    ''      using 1:3 index 1 with linespoints linestyle 2 t "dt = 0.1", \
    ''      using 1:3 index 2 with linespoints linestyle 3 t "dt = 1" 

set output "2_Q.png"
plot '2.dat' using 1:2 index 0 with linespoints linestyle 1 t "Q = f(t), 0.5*omega", \
    ''      using 1:2 index 1 with linespoints linestyle 2 t "Q = f(t), 0.8*omega", \
    ''      using 1:2 index 2 with linespoints linestyle 3 t "Q = f(t), 1.0*omega", \
    ''      using 1:2 index 3 with linespoints linestyle 4 t "Q = f(t), 1.2*omega" 

set output "2_I.png"
plot '2.dat' using 1:3 index 0 with linespoints linestyle 1 t "I = f(t), 0.5*omega", \
    ''      using 1:3 index 1 with linespoints linestyle 2 t "I = f(t), 0.8*omega", \
    ''      using 1:3 index 2 with linespoints linestyle 3 t "I = f(t), 1.0*omega", \
    ''      using 1:3 index 3 with linespoints linestyle 4 t "I = f(t), 1.2*omega" 