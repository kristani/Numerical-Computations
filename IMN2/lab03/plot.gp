#!/usr/bin/env gnuplot

set term png size 1900,1080

set xlabel 'x'
set ylabel 'u(x)'

set out 'results/wynik.png'
set multiplot layout 2,3

set title "M = 5"
plot 'results/M_5/u_0.txt' u 1:2 w l lw 2 title "mu = 0", 'results/M_5/u_1.txt' u 1:2 w l lw 2 title "mu = 1", 'results/M_5/u_2.txt' u 1:2 w l lw 2 title "mu = 2", 'results/M_5/u_3.txt' u 1:2 w l lw 2 title "mu = 3", 'results/M_5/u_4.txt' u 1:2 w l lw 2 title "mu = 4"

set title "M = 10"
plot 'results/M_10/u_0.txt' u 1:2 w l lw 2 title "mu = 0", 'results/M_10/u_1.txt' u 1:2 w l lw 2 title "mu = 1", 'results/M_10/u_2.txt' u 1:2 w l lw 2 title "mu = 2", 'results/M_10/u_3.txt' u 1:2 w l lw 2 title "mu = 3", 'results/M_10/u_4.txt' u 1:2 w l lw 2 title "mu = 4"

set title "M = 30"
plot 'results/M_30/u_0.txt' u 1:2 w l lw 2 title "mu = 0", 'results/M_30/u_1.txt' u 1:2 w l lw 2 title "mu = 1", 'results/M_30/u_2.txt' u 1:2 w l lw 2 title "mu = 2", 'results/M_30/u_3.txt' u 1:2 w l lw 2 title "mu = 3", 'results/M_30/u_4.txt' u 1:2 w l lw 2 title "mu = 4"

set xlabel 'alpha'
set ylabel 'E'

set title "M = 5"
plot 'results/M_5/E.txt' u 1:2 w l lw 2 title "mu = 0", 'results/M_5/E.txt' u 1:3 w l lw 2 title "mu = 1", 'results/M_5/E.txt' u 1:4 w l lw 2 title "mu = 2", 'results/M_5/E.txt' u 1:5 w l lw 2 title "mu = 3", 'results/M_5/E.txt' u 1:6 w l lw 2 title "mu = 4"

set title "M = 10"
plot 'results/M_10/E.txt' u 1:2 w l lw 2 title "mu = 0", 'results/M_10/E.txt' u 1:3 w l lw 2 title "mu = 1", 'results/M_10/E.txt' u 1:4 w l lw 2 title "mu = 2", 'results/M_10/E.txt' u 1:5 w l lw 2 title "mu = 3", 'results/M_10/E.txt' u 1:6 w l lw 2 title "mu = 4"

set title "M = 30"
plot 'results/M_30/E.txt' u 1:2 w l lw 2 title "mu = 0", 'results/M_30/E.txt' u 1:3 w l lw 2 title "mu = 1", 'results/M_30/E.txt' u 1:4 w l lw 2 title "mu = 2", 'results/M_30/E.txt' u 1:5 w l lw 2 title "mu = 3", 'results/M_30/E.txt' u 1:6 w l lw 2 title "mu = 4"

unset out
