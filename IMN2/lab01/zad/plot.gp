#!/usr/bin/env gnuplot

set term png size 800,600

set grid
set xlabel 'E'
set ylabel 'R_n(E)'

set xrange[0:150]

set out 'zad1.png'
plot 0*x title "0", 'zad1.txt' u 1:2 w l lw 2 title "R"
unset out

unset xrange
set out 'zad3.png'
plot 0*x title "0", 'zad3_E0.txt' u 1:2 w l lw 2 title "E0", 'zad3_E1.txt' u 1:2 w l lw 2 title "E1", 'zad3_E2.txt' u 1:2 w l lw 2 title "E2", 'zad3_E3.txt' u 1:2 w l lw 2 title "E3"
unset out


set out 'zad4.png'
plot 0*x title "0", 'zad4_E0.txt' u 1:2 w l lw 2 title "E0", 'zad4_E1.txt' u 1:2 w l lw 2 title "E1", 'zad4_E2.txt' u 1:2 w l lw 2 title "E2", 'zad4_E3.txt' u 1:2 w l lw 2 title "E3"
unset out
