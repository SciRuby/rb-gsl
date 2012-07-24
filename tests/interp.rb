#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "interpolation/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
include Math

def test_bsearch()
  status = 0
  x_array = GSL::Vector.alloc(0.0, 1.0, 2.0, 3.0, 4.0)
  index_result = Interp.bsearch(x_array, 1.5, 0, 4)
  s = (index_result != 1) ? 1 : 0
  status += 1
  test(s, "simple bsearch")

  index_result = x_array.bsearch(4.0, 0, 4)
  s = (index_result != 3) ? 1 : 0
  status += s;
  test(s, "upper endpoint bsearch");

  index_result = Interp.bsearch(x_array, 0.0, 0, 4);
  s = (index_result != 0) ? 1 : 0
  status += s;
  test(s, "lower endpoint bsearch");

  index_result = Interp.bsearch(x_array, 2.0, 0, 4);
  s = (index_result != 2) ? 1 : 0
  status += s;
  test(s, "degenerate bsearch");

  index_result = Interp.bsearch(x_array, 10.0, 0, 4);
  s = (index_result != 3) ? 1 : 0
  status += s;
  test(s, "out of bounds bsearch +");

  index_result = Interp.bsearch(x_array, -10.0, 0, 4);
  s = (index_result != 0) ? 1 : 0
  status += s;
  test(s, "out of bounds bsearch -");

end



test_bsearch()
