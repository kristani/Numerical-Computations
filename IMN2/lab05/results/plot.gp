#!/usr/bin/env gnuplot

set term png size 800,800

set grid
set xlabel 'x'
set ylabel 'y'

unset key
set out "mesh.png"
set title 'n_x = n_y = 10'
plot "mesh.txt" u 1:2 w l lt -1 lw 2  t '' , "" u ($1+0.1):($2+0.1):3 w labels t ''  
unset out

set term png size 800,800
set style increment userstyles
set contour
set cntrparam levels 30
set view map
unset surface
unset key
unset grid
set xrange [-5:5]
set yrange [-5:5]

set palette rgb 33,13,10

set out "u_10.png"
set title 'n_x = n_y = 10'
splot "u_10.txt" u 1:2:3 w l lt -1 lw 3 palette t ''
unset out