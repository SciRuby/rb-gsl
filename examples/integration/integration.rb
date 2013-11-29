#!/usr/bin/env ruby
require("gsl")
include Math

proc = Proc.new{ |x, alpha|
  log(alpha*x)/sqrt(x)
}

f = GSL::Function.alloc(proc)
f.set_params(1.0)

expected = -4.0

#result, error, neval = f.integration_qags([0.0, 1.0], 0.0, 1.0e-7, 1000)
#result, error, neval = f.integration_qags([0.0, 1.0])
result, error, neval = f.qags([0.0, 1.0])

printf("result          = %.18f\n", result);
printf("exact result    = %.18f\n", expected);
printf("estimated error = %.18f\n", error);
printf("actual error    = %.18f\n", result - expected);
printf("intervals =  %d\n", neval);
