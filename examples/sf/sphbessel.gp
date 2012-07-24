set multiplot
set size 0.5, 0.5
set xrange [0:12]
set grid
set xlabel 'x'

set origin 0, 0.5
set yrange [-0.5:1.2]
set ylabel 'j(x)'
plot 'sphbessel.dat' u 1:2 title 'j0' w l lw 2, '' u 1:3 title 'j1' w l lw 2, '' u 1:4 title 'j2' w l lw 2, '' u 1:5 title 'j3' w l lw 2

set origin 0.5, 0.5
set ylabel 'y(x)'
set yrange [-1:0.6]
plot 'sphbessel.dat' u 1:6 title 'y0' w l lw 2, '' u 1:7 title 'y1' w l lw 2, '' u 1:8 title 'y2' w l lw 2

set xrange [0:4]
set yrange [0:1]
set origin 0, 0
set ylabel 'exp(-|x|) i(x)'
plot 'sphbessel.dat' u 1:9 title 'i0' w l lw 2, '' u 1:10 title 'i1' w l lw 2, '' u 1:11 title 'i2' w l lw 2

set origin 0.5, 0
set yrange [0:10]
set ylabel 'exp(x) k(x)'
plot 'sphbessel.dat' u 1:12 title 'k0' w l lw 2, '' u 1:13 title 'k1' w l lw 2, '' u 1:14 title 'k2' w l lw 2

