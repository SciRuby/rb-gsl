set multiplot
set size 0.33, 0.33

set origin 0, 0.66
plot 'rk2.dat' title 'rk2', exp(-2*x)
set origin 0.33, 0.66
plot 'rk4.dat' title 'rk4', exp(-2*x)
set origin 0.66, 0.66
plot 'rkf45.dat' title 'rkf45', exp(-2*x)
set origin 0, 0.33
plot 'rkck.dat' title 'rkck', exp(-2*x)
set origin 0.33, 0.33
plot 'rk8pd.dat' title 'rk8pd', exp(-2*x)
set origin 0.66, 0.33
plot 'rk2imp.dat' title 'rk2imp', exp(-2*x)
set origin 0, 0
plot 'rk4imp.dat' title 'rk4imp', exp(-2*x)
set origin 0.33, 0
plot 'bsimp.dat' title 'bsimp', exp(-2*x)
set origin 0.66, 0
plot 'gear1.dat' title 'gear1', exp(-2*x)
#set origin 0.5, 0
#plot 'gear2.dat' title 'gear2', exp(-2*x)
