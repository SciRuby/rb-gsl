#!/usr/bin/env ruby
# Solve
#   dydt = -2y     ---> y(t) = exp(-2t)

require("gsl")
include GSL::Odeiv

dim = 1

PARAM = 2.0

func = Proc.new { |t, y, dydt, param|
  dydt[0] = -param*y[0]
}

jac = Proc.new { |t, y, dfdy, dfdt, param|
  dfdy = -param
  dfdt[0] = 0.0
}

def odeiv_compare_algorithm(solver, steptype, t0, tend, h0, y0, outfile)
  t = t0*1.0
  h = h0*1.0
  y = GSL::Vector.alloc([y0])
  s = Step.alloc(steptype, solver.dim)
  solver.set_step(s)
  solver.evolve.reset
#  p solver.params

  i = 0
  file = File.open(outfile, "w")
  while t < tend
    t, h, status = solver.apply(t, tend, h, y)
    break if status != GSL::SUCCESS
    file.printf("%.5e %.5e %.5e\n", t, y[0], h)
    i += 1
  end
  file.close
  return i
end

if GSL::VERSION >= "1.5.90"
  ALGORITHMS = ["rk2", "rk4", "rkf45", "rkck", "rk8pd", "rk2imp", "rk4imp",
                "bsimp", "gear1", "gear2", "rk2simp"]
  gpfile = "demo2.gp"
else
  ALGORITHMS = ["rk2", "rk4", "rkf45", "rkck", "rk8pd", "rk2imp", "rk4imp",
                "bsimp", "gear1", "gear2"]
  gpfile = "demo.gp"
end

solver = Solver.alloc(Step::RKF45, [1e-6, 0.0], func, jac, dim)
solver.set_params(PARAM)

puts("Solve dydt = -2y, y[0] = 1")
puts("y(t) = exp(-2t)\n\n")
ALGORITHMS.each do |alg|
  outfile = alg + ".dat"
  t0, tend, h0, y0 = 0.0, 3.0, 1e-6, 1.0
  i = odeiv_compare_algorithm(solver, alg, t0, tend, h0, y0, outfile)
  printf("%7s: Iteration %5d\n", alg, i)
end

system("gnuplot -persist #{gpfile}")

ALGORITHMS.each do |alg|
  outfile = alg + ".dat"
  File.delete(outfile)
end
