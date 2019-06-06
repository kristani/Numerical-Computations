#!/usr/bin/env gnuplot

unset grid
unset key

reset

set palette rgb 33,13,10

set term gif animate size 800,800
set output "u_animate.gif"
n=100    #n frames
set pm3d
set view map
set size ratio -1
set xrange[-5:5]
set yrange[-5:5]
set zr [-0.25:0.25]
set cbr [-0.25:0.25]

i=1
load "plot_prepare.gp"
set output