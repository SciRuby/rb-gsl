#!/usr/bin/env ruby
require("gsl")

unless GSL::MultiFit.const_defined?("Ndlinear")
  puts("The extension library NDLINEAR is not installed.")
  exit()
end

N_DIM = 3
N_SUM_R = 10
N_SUM_THETA = 11
N_SUM_PHI = 9
R_MAX = 3.0

def psi_real_exact(k, l, m, r, theta, phi)
  rr = GSL::pow(r, l)*Math::exp(-r*r)*GSL::Sf::laguerre_n(k, l + 0.5, 2 * r * r)   

  tt = GSL::Sf::legendre_sphPlm(l, m, Math::cos(theta))

  pp = Math::cos(m*phi)

  rr*tt*pp
end

basis_r = Proc.new { |r, y, params|
  params.eval(r, y)
}

basis_theta = Proc.new { |theta, y, params|
  for i in 0...N_SUM_THETA do
    y[i] = GSL::Sf::legendre_Pl(i, Math::cos(theta));
  end
}

basis_phi = Proc.new { |phi, y, params|
  for i in 0...N_SUM_PHI do
    if i%2 == 0
      y[i] = Math::cos(i*0.5*phi)
    else
      y[i] = Math::sin((i+1.0)*0.5*phi)
    end
  end
}


GSL::Rng::env_setup()

k = 5
l = 4
m = 2

NDATA = 3000

N = [N_SUM_R, N_SUM_THETA, N_SUM_PHI]
u = [basis_r, basis_theta, basis_phi]

rng = GSL::Rng.alloc()

bspline = GSL::BSpline.alloc(4, N_SUM_R - 2)
bspline.knots_uniform(0.0, R_MAX)

ndlinear = GSL::MultiFit::Ndlinear.alloc(N_DIM, N, u, bspline)
multifit = GSL::MultiFit.alloc(NDATA, ndlinear.n_coeffs)
vars = GSL::Matrix.alloc(NDATA, N_DIM)
data = GSL::Vector.alloc(NDATA)


for i in 0...NDATA do

  r = rng.uniform()*R_MAX
  theta = rng.uniform()*Math::PI
  phi = rng.uniform()*2*Math::PI

  psi = psi_real_exact(k, l, m, r, theta, phi)

  dpsi = rng.gaussian(0.05*psi)

  vars[i,0] = r
  vars[i,1] = theta
  vars[i,2] = phi    

  data[i] = psi + dpsi
end

#GSL::MultiFit::Ndlinear::design(vars, X, ndlinear)
X = GSL::MultiFit::Ndlinear::design(vars, ndlinear)

coeffs, cov, chisq, = GSL::MultiFit::linear(X, data, multifit)

rsq = GSL::MultiFit::linear_Rsq(data, chisq)
STDERR.printf("chisq = %e, Rsq = %f\n", chisq, rsq)

eps_rms = 0.0
volume = 0.0
dr = 0.05;
dtheta = 5.0 * Math::PI / 180.0
dphi = 5.0 * Math::PI / 180.0
x = GSL::Vector.alloc(N_DIM)

r = 0.01
while r < R_MAX do
  theta = 0.0
  while theta < Math::PI do
    phi = 0.0
    while phi < 2*Math::PI do
      dV = r*r*Math::sin(theta)*dr*dtheta*dphi
      x[0] = r
      x[1] = theta
      x[2] = phi

      psi_model, err = GSL::MultiFit::Ndlinear.calc(x, coeffs, cov, ndlinear)
      psi = psi_real_exact(k, l, m, r, theta, phi)
      err = psi_model - psi
      eps_rms += err * err * dV;
      volume += dV;

      if phi == 0.0
         printf("%e %e %e %e\n", r, theta, psi, psi_model)
      end

      phi += dphi
    end
    theta += dtheta
  end
  printf("\n");

  r += dr
end

eps_rms /= volume
eps_rms = Math::sqrt(eps_rms)
STDERR.printf("rms error over all parameter space = %e\n", eps_rms)

