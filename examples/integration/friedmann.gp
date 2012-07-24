#set terminal postscript eps enhanced color
#set output 'tmp.eps'
set xlabel 'Redshift z'
set ylabel 'Luminosity Distance [1/H0]'
set ytics nomirror
set grid
set y2tics
set yrange [0:13]
set y2range [0:54.89]
#set logscale xy
#set logscale y2
#set yrange [0.01:11]
#set y2range [0.04:60]

set y2label 'Luminosity Distance [Gpc] (h = 0.71)'
plot 'friedmann.dat' u 1:3 title 'Matter = 1, Lambda = 0 (Einstein - de Sitter)' w l lw 3, '' u 1:5 title 'Matter = 0.27, Lambda = 0.73 (WMAP, Lambda-CDM)' w l lw 3, '' u 1:7 title 'Matter = 0, Lambda = 1 (de Sitter)' w l lw 3
