#!/usr/bin/env gnuplot

set term png size 800,600

set grid
set xlabel 'x'
set ylabel 'Blad'

set out 'zad1/kolokacji.png'
set title 'Metoda kolokacji'
plot 'zad1/kolokacji_6.txt' u 1:2 w l lw 1 title "N=6", 'zad1/kolokacji_7.txt' u 1:2 w l lw 1 title "N=7", 'zad1/kolokacji_8.txt' u 1:2 w l lw 1 title "N=8", 'zad1/kolokacji_9.txt' u 1:2 w l lw 1 title "N=9", 'zad1/kolokacji_10.txt' u 1:2 w l lw 1 title "N=10"
unset out


set out 'zad2/kwadratow.png'
set title 'Metoda najmniejszych kwadratow'
plot 'zad2/kwadratow_6.txt' u 1:2 w l lw 1 title "N=6", 'zad2/kwadratow_7.txt' u 1:2 w l lw 1 title "N=7", 'zad2/kwadratow_8.txt' u 1:2 w l lw 1 title "N=8", 'zad2/kwadratow_9.txt' u 1:2 w l lw 1 title "N=9", 'zad2/kwadratow_10.txt' u 1:2 w l lw 1 title "N=10"
unset out


set out 'zad3/galerkina.png'
set title 'Metoda Galerkina'
plot 'zad3/galerkin_6.txt' u 1:2 w l lw 1 title "N=6", 'zad3/galerkin_7.txt' u 1:2 w l lw 1 title "N=7", 'zad3/galerkin_8.txt' u 1:2 w l lw 1 title "N=8", 'zad3/galerkin_9.txt' u 1:2 w l lw 1 title "N=9", 'zad3/galerkin_10.txt' u 1:2 w l lw 1 title "N=10"
unset out


set out 'zad4/kolokacji.png'
set title 'Metoda kolokacji'
plot 'zad4/kolokacji_6.txt' u 1:2 w l lw 1 title "N=6", 'zad4/kolokacji_7.txt' u 1:2 w l lw 1 title "N=7", 'zad4/kolokacji_8.txt' u 1:2 w l lw 1 title "N=8", 'zad4/kolokacji_9.txt' u 1:2 w l lw 1 title "N=9", 'zad4/kolokacji_10.txt' u 1:2 w l lw 1 title "N=10"
unset out


set out 'zad4/kwadratow.png'
set title 'Metoda najmniejszych kwadratow'
plot 'zad4/kwadratow_6.txt' u 1:2 w l lw 1 title "N=6", 'zad4/kwadratow_7.txt' u 1:2 w l lw 1 title "N=7", 'zad4/kwadratow_8.txt' u 1:2 w l lw 1 title "N=8", 'zad4/kwadratow_9.txt' u 1:2 w l lw 1 title "N=9", 'zad4/kwadratow_10.txt' u 1:2 w l lw 1 title "N=10"
unset out


set out 'zad4/galerkina.png'
set title 'Metoda Galerkina'
plot 'zad4/galerkin_6.txt' u 1:2 w l lw 1 title "N=6", 'zad4/galerkin_7.txt' u 1:2 w l lw 1 title "N=7", 'zad4/galerkin_8.txt' u 1:2 w l lw 1 title "N=8", 'zad4/galerkin_9.txt' u 1:2 w l lw 1 title "N=9", 'zad4/galerkin_10.txt' u 1:2 w l lw 1 title "N=10"
unset out