#!/usr/bin/env ruby
require("gsl")
require("ool/conmin.rb")

def iteration_echo(m)
  f = m.f
  x = m.x
  printf( "f( " );
  for i in 0...3 do
    printf( "%+6.3e, ", x[i])
  end
  printf( "... ) = %+6.3e\n", f )

end

NN = 100

fun = Proc.new { |x, params|
  n = x.size
  f = 0
  for i in 0...n do
    xi = x[i] - i.to_f/10.0
    f += (i + 1.0)*xi*xi
  end
  f
}

fun_df = Proc.new { |x, params, g|
  n = x.size
  for i in 0...n do
    xi = x[i] - i.to_f/10.0
    g[i] = 2.0*(i+1.0)*xi
  end
}

fun_Hv = Proc.new { |x, params, v, hv|
  n = x.size
  for i in 0...n do
    hv[i] = 2.0*(i+1.0)*v[i]
  end
}

a = GSL::Vector.alloc(NN)
for i in 0...NN do
  a[i] = (i.to_f + 1.0)/10.0
end
f = OOL::Conmin::Function.alloc()
f.set(NN, fun, fun_df, nil, fun_Hv, a)

L = GSL::Vector.alloc(NN)
U = GSL::Vector.alloc(NN)
L.set_all(-3.0)
U.set_all(3.0)
c = OOL::Conmin::Constraint.alloc()
c.set(NN, L, U)

#m = OOL::Conmin::Minimizer.alloc(OOL::Conmin::Minimizer::Spg, NN)
m = OOL::Conmin::Minimizer.alloc(OOL::Conmin::Minimizer::Pgrad, NN)
#m = OOL::Conmin::Minimizer.alloc("pgrad", NN)
#m = OOL::Conmin::Minimizer.alloc("spg", NN)

#params = OOL::Conmin::Minimizer::Spg.parameters_default
params = m.parameters_default
p params.class
p params

x = GSL::Vector.alloc(NN)
for i in 0...NN do
  x[i] = 1.0 + i
end
m.set(f, c, x, params)

ii = 0
NMAX = 1000
status = OOL::CONTINUE;

printf( "%4d : ", ii )
iteration_echo ( m )
while ii < NMAX && status == OOL::CONTINUE
  ii+=1
  m.iterate
  status = m.is_optimal

  printf( "%4d : ", ii )
  iteration_echo( m )
end

printf("%s method\n", m.name)
if status == OOL::SUCCESS
  printf("\nConvergence in %i iterations\n\n", ii);
else
  printf("Stopped with %i iterations\n", ii);
end
printf("variables................: %6d\n", NN)
printf("function evaluations.....: %6d\n", m.fcount)
printf("gradient evaluations.....: %6d\n", m.gcount)
printf("function value...........: %.6e\n", m.minimum)
printf("projected gradient norm..: %.6e\n", m.size)


