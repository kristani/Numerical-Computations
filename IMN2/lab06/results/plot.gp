#!/usr/bin/env gnuplot
unset multiplot
set term png size 1200,800

set grid
set xlabel 't'
set ylabel ''
unset key

set out "zad2/wartosci.png"
set multiplot layout 2,2

set title 'Wartosc 1'
plot "zad2/val_1.txt" u 1:2 w l lw 2 

set title 'Wartosc 2'
plot "zad2/val_2.txt" u 1:2 w l lw 2 

set title 'Wartosc 3'
plot "zad2/val_3.txt" u 1:2 w l lw 2 

set title 'Wartosc 4'
plot "zad2/val_4.txt" u 1:2 w l lw 2 
unset multiplot 

unset out


set term png size 800,800
set xlabel 'x'
set ylabel 'y'

unset key
set out "mesh.png"
set title 'n_x = n_y = 10'
plot "mesh.txt" u 1:2 w l lt -1 lw 2  t '' , "" u ($1+0.1):($2+0.1):3 w labels t ''  
unset out

set term png size 800,800
set view map
set surface
unset key
unset grid
set xrange [-5:5]
set yrange [-5:5]

set palette rgb 33,13,10

set out "zad1/u_0.png"
set title 'n_x = n_y = 10'
splot "zad1/u_0.txt" u 1:2:3 palette w pm3d t ''
unset out

set out "zad1/u_1.png"
set title 'n_x = n_y = 10'
splot "zad1/u_1.txt" u 1:2:3 palette w pm3d t ''
unset out

set out "zad1/u_2.png"
set title 'n_x = n_y = 10'
splot "zad1/u_2.txt" u 1:2:3 palette w pm3d t ''
unset out

set out "zad1/u_3.png"
set title 'n_x = n_y = 10'
splot "zad1/u_3.txt" u 1:2:3 palette w pm3d t ''
unset out

set out "zad1/u_4.png"
set title 'n_x = n_y = 10'
splot "zad1/u_4.txt" u 1:2:3 palette w pm3d t ''
unset out

set out "zad1/u_5.png"
set title 'n_x = n_y = 10'
splot "zad1/u_5.txt" u 1:2:3 palette w pm3d t ''
unset out

set out "zad1/u_6.png"
set title 'n_x = n_y = 10'
splot "zad1/u_6.txt" u 1:2:3 palette w pm3d t ''
unset out

set out "zad1/u_7.png"
set title 'n_x = n_y = 10'
splot "zad1/u_7.txt" u 1:2:3 palette w pm3d t ''
unset out

set out "zad1/u_8.png"
set title 'n_x = n_y = 10'
splot "zad1/u_8.txt" u 1:2:3 palette w pm3d t ''
unset out

set out "zad1/u_9.png"
set title 'n_x = n_y = 10'
splot "zad1/u_9.txt" u 1:2:3 palette w pm3d t ''
unset out

set terminal gif small animate size 800, 800
set pm3d
set view map
set size ratio -1

set output "zad3/u_animation.gif"
do for [i=0:50] {
splot "zad3/u.txt" in i u 1:2:3 w pm3d not
}

unset pm3d
unset view map