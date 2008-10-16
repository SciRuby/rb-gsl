set multiplot
set size 0.33, 0.25

set origin 0, 0.75
plot 'rk2.dat' title 'rk2', exp(-2*x)
set origin 0.33, 0.75
plot 'rk4.dat' title 'rk4', exp(-2*x)
set origin 0.66, 0.75
plot 'rkf45.dat' title 'rkf45', exp(-2*x)
set origin 0, 0.5
plot 'rkck.dat' title 'rkck', exp(-2*x)
set origin 0.33, 0.5
plot 'rk8pd.dat' title 'rk8pd', exp(-2*x)
set origin 0.66, 0.5
plot 'rk2imp.dat' title 'rk2imp', exp(-2*x)
set origin 0, 0.25
plot 'rk4imp.dat' title 'rk4imp', exp(-2*x)
set origin 0.33, 0.25
plot 'bsimp.dat' title 'bsimp', exp(-2*x)
set origin 0.66, 0.25
plot 'gear1.dat' title 'gear1', exp(-2*x)
set origin 0, 0
plot 'gear2.dat' title 'gear2', exp(-2*x)
set origin 0.33, 0
plot 'rk2simp.dat' title 'rk2simp', exp(-2*x)
