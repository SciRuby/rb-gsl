set multiplot
set size 0.5, 0.5
set grid
set key box
set xrange [-5.5:5.5]

# Hermite Polynomials
H0(x) = 1
H1(x) = 2*x
H2(x) = 4*x*x-2
H10(x) = -30240 + x*(0 + x*(302400 + x*(0 + x*(-403200 + x*(0 + x*(161280 + x*(0 + x*(-23040 + x*(0 + x*1024)))))))))

# Normalization constant
coef(n) = sqrt(1.0/2**n/gamma(n+1)/sqrt(pi))

psi0(x) = coef(0)*exp(-x*x/2)*H0(x)
psi1(x) = coef(1)*exp(-x*x/2)*H1(x)
psi2(x) = coef(2)*exp(-x*x/2)*H2(x)
psi10(x) = coef(10)*exp(-x*x/2)*H10(x)

set ylabel 'psi(x)'

set pointsize 1

set origin 0, 0.5
plot psi0(x) title 'Exact: n = 0', "qhoscillator.dat" u 1:2 pt 7 title 'Numerical'

set origin 0.5, 0.5
plot psi1(x) title 'Exact: n = 1', "qhoscillator.dat" u 1:3 pt 7 title 'Numerical'

set origin 0, 0
plot psi2(x) title 'Exact: n = 2', "qhoscillator.dat" u 1:4 pt 7title 'Numerical'

set origin 0.5, 0
plot psi10(x) title 'Exact: n = 10', "qhoscillator.dat" u 1:5 pt 7 title 'Numerical'
