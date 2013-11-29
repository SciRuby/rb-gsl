set pointsize 2
set multiplot
set xlabel 'x'
set ylabel 'y'
set yrange [0:14]
set size 0.5, 0.5
set grid
set origin 0, 0.5
set key box
plot 'points'

set origin 0.5, 0.5
plot 'points', 'linear.dat' w l, 'polynomial.dat' w l

set origin 0, 0
plot 'points', 'cspline.dat' w l, 'cspline_periodic.dat' w l

set origin 0.5, 0
plot 'points', 'akima.dat' w l, 'akima_periodic.dat' w l

