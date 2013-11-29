#!/usr/bin/env ruby
require("gsl")

exit unless GSL.const_defined?("ALF")

lmax = 3
#w = GSL::ALF::Workspace.alloc(lmax)
#w = GSL::ALF.alloc(lmax)
w = GSL::alf_alloc(lmax)
w.params(0, 0, GSL::ALF::NORM_NONE)


result = GSL::Vector.alloc(GSL::ALF::array_size(lmax))
deriv = GSL::Vector.alloc(GSL::ALF::array_size(lmax))

ind11 = GSL::ALF::array_index(1, 1)
ind21 = GSL::ALF::array_index(2, 1)
ind22 = GSL::ALF::array_index(2, 2)
ind31 = GSL::ALF::array_index(3, 1)

File.open("alf.dat", "w") do |f|
  theta = 0.01
  while theta < 1.99*Math::PI do
    x = Math.cos(theta)
#    w.Plm_deriv_array(lmax, x, result, deriv)
    w.Plm_array(lmax, x, result, deriv)
    f.printf("%e %e %e %e %e\n", theta, result[ind11], result[ind21], result[ind22], result[ind31])
    theta += 0.01
  end
end

system("gnuplot -persist alf.gp")
