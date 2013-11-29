#!/usr/bin/gnuplot

set terminal postscript eps enhanced color 18

set pm3d map
#set palette gray
set palette
#set palette gamma 0.5
set xrange [0:3]
set yrange [0:3]
set cbrange [-4:5]
set ytics border nomirror

set output "ndlinear.eps"
set xlabel "r"
set ylabel "{/Symbol \161}"

set colorbox horizontal user origin 0.1,0.05 size 1.8,0.05

set size 2.0, 1.1
set multiplot

set size 1.0, 1.0
set origin 0, 0.1
set title "Exact solution"
unset key
splot 'ndlinear.dat' us 1:2:3

unset colorbox

set size 1.0, 1.0
set origin 1, 0.1
set title "Model solution"
splot 'ndlinear.dat' us 1:2:4

unset multiplot
