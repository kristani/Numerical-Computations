#!/usr/bin/env gnuplot

set term png size 1200,800

set grid
set xlabel 'x'
set ylabel 'y'

unset key
set out "siatka.png"
set title 'n_x = n_y = 3'
plot "siatka.txt" u 1:2 w p pt 7 ps 3 
unset out

set term png size 800,800
set style increment userstyles
set contour
set cntrparam levels 40
set view map
unset surface
unset key
unset grid
set xrange [0:3.14]
set yrange [0:3.14]

set palette rgb 33,13,10

set out "u_nx_3.png"
set title 'n_x = n_y = 3'
splot "u_nx_3.txt" u 1:2:3 w l lt -1 lw 3 palette t ''
unset out

set out "u_nx_10.png"
set title 'n_x = n_y = 10'
splot "u_nx_10.txt" u 1:2:3 w l lt -1 lw 3 palette t ''
unset out