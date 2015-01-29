set logscale y
set xtics nomirror
set yrange [1e-5:1e1]
set xrange [2100:900]
set x2range [2100*2.725:900*2.725]
set x2tics
set xlabel '1 + z (Redshift)'
set ylabel 'Fractional ionization'
set x2label 'Temperature of the Universe [K]'
set grid
plot 'recombination.dat' u 2:3 title '' w l lw 3, 1 title ''
