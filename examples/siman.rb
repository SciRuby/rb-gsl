#!/usr/bin/env ruby
require("gsl")
include GSL
include GSL::Siman

N_TRIES = 200
ITERS_FIXED_T = 10
STEP_SIZE = 10
K = 1.0
T_INITIAL = 0.002
MU_T = 1.005
T_MIN = 2.0e-6

params = Siman::Params.alloc(N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL,
                           MU_T, T_MIN)

efunc = Efunc.alloc { |vx|
  x = vx[0]
  s = (x - 1.0)
  Math::exp(-s*s)*Math::sin(8*x)
}

metric = Metric.alloc { |vx, vy|
  (x[0] - y[0]).abs
}

step = Step.alloc { |rng, vx, step_size|
  r = rng.uniform
  old_x = vx[0]
  a =  r * 2 * step_size - step_size + old_x
  vx[0] = a
}

simanprint = Print.alloc { |vx|
  printf("%12g", vx[0])
}

x = GSL::Vector.alloc([15.5])
GSL::Rng.env_setup()
r = GSL::Rng.alloc()

Siman::solve(r, x, efunc, step, metric, simanprint, params)
#Siman::solve(r, x, efunc, step, metric, nil, params)
p x
