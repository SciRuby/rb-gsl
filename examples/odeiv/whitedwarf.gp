#set terminal postscript eps enhanced color
#set output 'whitedwarf.eps'
set multiplot
set size 0.5, 0.5

set grid
set origin 0, 0.5
#set xlabel 'r/R'
#set ylabel '{/Symbol r}/{/Symbol r}_c'
set xlabel 'Radius (normalized)'
set ylabel 'Density (normalized)'
set xrange [0:1]
plot 'polytrope.dat' u ($1/6.8969):4 title 'Density Dist' w l lt 1 lw 5

set origin 0.5, 0.5
set logscale x
set xrange [1e6:1e11]
#set xlabel '{/Symbol r}_c [g/cm^3]'
set xlabel 'Central Density [g/cm^3]'
set ylabel 'Radius [km]'
#plot 'whitedwarf.dat' u 1:3 title '{/Symbol r}_c - Radius' w l lw 5
plot 'whitedwarf.dat' u 1:3 title 'Central Density - Radius' w l lt 2 lw 5

set origin 0, 0
set logscale x
set xrange [1e6:1e11]
set yrange [0.5:1.6]
#set xlabel '{/Symbol r}_c [g/cm^3]'
set xlabel 'Central Density [g/cm^3]'
set ylabel 'WD Mass [Msolar]'
plot 'whitedwarf.dat' u 1:4 title 'Central Density - Mass' w l lt 3 lw 5

set origin 0.5, 0
unset logscale x
set xrange [0.5:1.5]
set yrange [100:25000]
set xlabel 'Mass [Msolar]'
set ylabel 'Radius [km]'
plot 'whitedwarf.dat' u 4:3 title 'Mass - Radius' w l lt 4 lw 5

