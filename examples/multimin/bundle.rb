#!/usr/bin/env ruby
require("gsl")

demyanov_f = Proc.new { |x, params|
  GSL::MAX_DBL(5.0*x[0]+x[1], GSL::MAX_DBL(x[0]*x[0]+x[1]*x[1]+4.0*x[1], -5.0*x[0] + x[1]))
}

demyanov_sdf = Proc.new { |x, sdf|
  f = 5.0*x[0] + x[1]
  ff = x[0]*x[0] + x[1]*x[1] + 4.0*x[1]
  i_max = 1
  if f < ff
    f = ff
    i_max = 2
  end
  if f < -5.0*x[0] + x[1]
    i_max = 3
  end
  case i_max
  when 1
    sdf[0] = 5.0
    sdf[1] = 1.0
  when 2
    sdf[0] = 2.0*x[0]
    sdf[1] = 2.0*x[1] + 4.0
  when 3
    sdf[0] = -5.0
    sdf[1] = 1.0
  end
}

max_iter = 1000
n = 2
function = GSL::MultiMin::Function_fsdf.alloc(demyanov_f, demyanov_sdf, n)

start_point = GSL::Vector.alloc(n)
start_point.set_all(1.0)
bundle_size_max = n + 3

s = GSL::MultiMin::FsdfMinimizer.alloc("bundle_method", n)
s.set(function, start_point, bundle_size_max)

printf("********************  %s  ********************\n\n","Demyanov function")

printf("== k ===== f(x) ===== ||sgr_f(x)|| ======= eps ======  \n");

subgradient = s.subgradient

iter = 0;
printf("%4d  %14.7f  %13.8e  %13.8e\n", iter, s.f, subgradient.dnrm2, s.eps)

begin
  iter += 1
  status = s.iterate
  status = s.test_convergence(1e-5)
  printf("%4d  %14.7f  %13.8e  %13.8e\n", iter, s.f, subgradient.dnrm2, s.eps)
  if status == GSL::SUCCESS
    printf("\nMinimum is found at\n")
    x = s.x
    for j in 0...x.size do
      printf("%9.6f ", x[j])
    end
    printf("\n\n")
  end
end while status == GSL::CONTINUE and iter <= max_iter

