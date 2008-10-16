#!/usr/bin/env ruby
require("gsl")
include GSL
include Monte
include Math

proc_f = Proc.new { |k, dim, params|
  pi = Math::PI
  a = 1.0/(pi*pi*pi)
  a/(1.0 - cos(k[0])*cos(k[1])*cos(k[2]))
}

def display_results(title, result, error)
  exact = 1.3932039296856768591842462603255

  diff = result - exact
  printf("%s ==================\n", title);
  printf("result = % .6f\n", result);
  printf("sigma  = % .6f\n", error);
  printf("exact  = % .6f\n", exact);
  printf("error  = % .6f = %.1g sigma\n", diff, diff.abs/error)
end

dim = 3
xl = Vector.alloc([0, 0, 0])
xu = Vector.alloc([PI, PI, PI])
G = Monte::Function.alloc(proc_f, dim)
calls = 500000
r = GSL::Rng.alloc(Rng::DEFAULT)

plain = Monte::Plain.alloc(dim)
result, error = G.integrate(xl, xu, dim, calls, r, plain)
display_results("plain", result, error)

miser = Monte::Miser.alloc(dim)
result, error = G.integrate(xl, xu, dim, calls, r, miser)
display_results("miser", result, error)

vegas = Monte::Vegas.alloc(dim)
result, error = G.integrate(xl, xu, dim, 10000, r, vegas)
display_results("vegas warm-up", result, error)
puts("converging...");
begin
  result, error = G.integrate(xl, xu, dim, calls/5, r, vegas)
  printf("result = % .6f sigma = % .6f chisq/dof = %.1f\n", result, error, vegas.chisq)
end while (vegas.chisq-1.0).abs > 0.5
display_results("vegas final", result, error)
