set multiplot
set size 0.5, 0.5
set grid

set origin 0, 0.5
plot "alf.dat" u 1:2, sin(x)

set origin 0.5, 0.5
plot "alf.dat" u 1:3, 3*cos(x)*sin(x)

set origin 0, 0
plot "alf.dat" u 1:4, 3*sin(x)*sin(x)

set origin 0.5, 0
plot "alf.dat" u 1:5, 1.5*sin(x)*(5*cos(x)*cos(x)-1)
