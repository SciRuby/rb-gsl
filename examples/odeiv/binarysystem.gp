set title 'Orbital decay of the binary system J0730-3039'
set xlabel 'Time since birth [Myr]'
set ylabel 'Orbital Period [Hr]'
set y2label 'Eccentricity'

set grid
set yrange [1e-2:10]
set xrange [90:200]
set y2range [1e-4:1e-1]
#set xrange [90:2500]
#set y2range [1e-4:2e-1]

set ytics nomirror
set y2tics
set logscale y
set logscale y2

set pointsize 1
set label 1 'Present' at first 96, 0.04
set label 2 'Plunge!' at first 150, 0.04
set arrow from 100,0.03 to 100,0.013 lw 2 lt 3
set arrow from 160, 0.035 to 182, 0.011 lw 2 lt 4
plot 'binarysystem.dat' u 1:2 title 'Orbital period' w lp lw 1, '' u 1:3 axes x1y2 title 'Eccentricity' w lp lw 1
