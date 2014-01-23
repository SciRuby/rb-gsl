require 'test_helper'

class MultifitTest < GSL::TestCase

  def _test_lmder(fdf, x, xx, f, cov)
    s = GSL::MultiFit::FdfSolver.alloc('lmsder', fdf.n, fdf.p)
    s.set(fdf, x)

    20.times { |i|
      s.iterate

      fdf.p.times { |j|
        assert_rel s.x[j], xx[fdf.p * i + j], 1e-5, "lmsder,  iter=#{i}, x#{j}"
      }

      assert_rel GSL::Blas.dnrm2(s.f), f[i], 1e-5, "lmsder, iter=#{i}, f"
    }

    covar = s.covar(0.0)

    fdf.p.times { |i|
      fdf.p.times { |j|
        assert_rel covar[i, j], cov[i * fdf.p + j], 1e-7, "gsl_multifit_covar cov(#{i},#{j})"
      }
    }
  end

  def _test_fdf(name, fdf, x, x_final, f_sumsq, sigma)
    s = GSL::MultiFit::FdfSolver.alloc('lmsder', fdf.n, fdf.p)
    s.set(fdf, x)

    1000.times {
      s.iterate

      status = s.test_delta(0.0, 1e-7)
      break if status != GSL::CONTINUE
    }

    covar = s.covar(0.0)

    fdf.p.times { |i|
      assert_rel s.x[i], x_final[i], 1e-5, "#{name}, lmsder, x#{i}"
    }

    s2 = GSL.pow(GSL::Blas.dnrm2(s.f), 2.0)
    assert_rel s2, f_sumsq, 1e-5, "#{name}, lmsder, |f|^2"

    fdf.p.times { |i|
      ei = Math.sqrt(s2 / (fdf.n - fdf.p)) * Math.sqrt(covar[i, i])
      assert_rel ei, sigma[i], 1e-4, "#{name}, sigma(#{i})"
    }
  end

  def test_2dgauss
    maxiter = 10
    n = 33

    point = Struct.new(:x, :y)

    # model: a * exp(-((x - x0) ** 2 + (y - y0) ** 2) / 2 / sigma ** 2)
    gauss_f = lambda { |x, t, y, s, f|
      # x: parameters as a Vecor
      # t: observed points as an Array
      # y: observed data as a GSL::Vector
      # s: errorbar
      # f: result
      a = x[0]
      x0 = x[1]
      y0 = x[2]
      sigma2 = x[3] ** 2

      y.size.times { |i|
        f.set(i, (a * Math.exp(-((t[i].x - x0) ** 2 + (t[i].y - y0) ** 2) / 2 / sigma2) - y[i]) / s[i])
      }

      GSL::SUCCESS
    }

    gauss_df = lambda { |x, t, y, s, df|
      a = x[0]
      x0 = x[1]
      y0 = x[2]
      sigma = x[3]
      sigma2 = sigma ** 2

      y.size.times { |i|
        dx = t[i].x - x0; dx2 = dx ** 2
        dy = t[i].y - y0; dy2 = dy ** 2

        f = a * Math.exp(-(dx2 + dy2) / 2 / sigma2)

        df.set(i, 0, f / a / s[i])
        df.set(i, 1, f * dx / sigma2 / s[i])
        df.set(i, 2, f * dy / sigma2 / s[i])
        df.set(i, 3, f * (dx2 + dy2) / sigma2 / sigma / s[i])
      }

      GSL::SUCCESS
    }

    # goal
    xgoal = GSL::Vector.alloc([1, 0, 0, 1])
    parname = %w[a x0 y0 si]

    # data
    t = []
    tmin = -10.0
    tmax =  10.0

    n.times { |j|
      n.times { |i|
        t << point.new(tmin + (tmax - tmin) * i / (n - 1), tmin + (tmax - tmin) * j / (n - 1))
      }
    }

    stdev = xgoal[0] * 0.1

    s = GSL::Vector.alloc(Array.new(t.size, stdev))  # error bar of each datum
    r = GSL::Rng.alloc
    e = GSL::Vector.alloc(t.size)

    t.size.times { |i|
      e[i] = -r.gaussian(stdev)  # perturbation to data
    }

    y = GSL::Vector.alloc(t.size)
    n = GSL::Vector.alloc(Array.new(t.size, 1.0))
    gauss_f.call(xgoal, t, e, n, y)  # data: y = model - e

    # fitting
    x = GSL::Vector.alloc([0.5, 0.1, -0.1, 2.0])  # initial guess
    fdf = GSL::MultiFit::Function_fdf.alloc(gauss_f, gauss_df, x.size)
    fdf.set_data(t, y, s)

    solver = GSL::MultiFit::FdfSolver.alloc(GSL::MultiFit::FdfSolver::LMSDER, t.size, x.size)
    solver.set(fdf, x)

    #solver.print_state(0)

    maxiter.times { |i|
      solver.iterate

      status = solver.test_delta(1e-6, 1e-6)
      #solver.print_state(i + 1)

      break if status != GSL::CONTINUE
    }

    # results
    covar = solver.covar(0.0)
    xresult = solver.position
    dof = t.size - xresult.size
    chi2 = GSL.pow_2(solver.f.dnrm2)
    xsigma = GSL::Vector.alloc(xresult.size)

    xresult.size.times { |i|
      xsigma[i] = Math.sqrt(chi2 / dof * covar[i, i]) * 2.0
      # allow resulting parameters to differ two times than standard error
    }

    desc = "a*exp(-((x-x0)**2+(y-y0)**2)/2/si**2), chi2/N:%.3g" % (chi2 / t.size)

    xresult.size.times { |i|
      assert_rel xresult[i], xgoal[i], xsigma[i], '%s: %-2.2s' % [desc, parname[i]]

      refute((xresult[i] - xgoal[i]).abs > xsigma[i],
        '%s: %-2.2s is %s +- %s' % [desc, parname[i], xresult[i], xsigma[i]])
    }
  end

  def test_brown
    brown_N = 20
    brown_P = 4

    brown_X = GSL::Matrix.alloc(
      [24.3485677, 4.71448798, -2.19486633, 2.69405755],
      [22.4116222, 3.93075538, -1.42344852, 2.5233557],
      [17.88886, 2.9290853, 0.125174936, -3.96823353],
      [17.3237176, 2.99606803, 2.03285653, 2.28992327],
      [17.0906508, 3.02485425, 0.296995153, 0.0876226126],
      [16.578006, 3.1036312, -0.18617941, 0.103262914],
      [15.692993, 3.33088442, 0.0706406887, 1.05923955],
      [14.3232177, 3.85604218, -2.3762839, -3.09486813],
      [14.1279266, 3.97896121, 0.446109351, 1.40023753],
      [13.6081961, 4.16435075, -1.51250057, -1.52510626],
      [13.4295245, 4.22697223, -0.196985195, 0.532009293],
      [13.0176117, 4.3579261, -0.353131208, 0.301377627],
      [12.2713535, 4.62398535, -0.00183585584, 0.894170703],
      [11.0316144, 5.13967727, -2.38978772, -2.89510064],
      [10.8807981, 5.24558004, 0.230495952, 1.27315117],
      [10.4029264, 5.41141257, -1.5116632, -1.47615921],
      [10.2574435, 5.46211045, -0.299855732, 0.451893162],
      [9.87863876, 5.57914292, -0.368885288, 0.358086545],
      [9.1894983, 5.82082741, -0.230157969, 0.621476534],
      [8.00589008, 6.27788753, -1.46022815, -1.33468082]
    )

    brown_F = GSL::Vector.alloc(
      2474.05541, 1924.69004, 1280.63194, 1244.81867,
      1190.53739, 1159.34935, 1108.44426, 1090.11073,
      1015.92942, 1002.43533, 971.221084, 949.589435,
      911.359899, 906.522994, 840.525729, 833.950164,
      807.557511, 791.00924, 761.09598, 726.787783
    )

    brown_cov = GSL::Matrix.alloc(
      [ 1.8893186910e-01, -4.7099989571e-02,  5.2154168404e-01,  1.6608168209e-02],
      [-4.7099989571e-02,  1.1761534388e-02, -1.2987843074e-01, -4.1615942391e-03],
      [ 5.2154168404e-01, -1.2987843074e-01,  1.4653936514e+00,  1.5738321686e-02],
      [ 1.6608168209e-02, -4.1615942391e-03,  1.5738321686e-02,  4.2348042340e-02]
    )

    brown_x0 = GSL::Vector.alloc(25, 5, -5, -1)

    brown_f = lambda { |x, t, y, f|
      brown_N.times { |i|
        ti = 0.2 * (i + 1)
        ui = x[0] + x[1] * ti - Math.exp(ti)
        vi = x[2] + x[3] * Math.sin(ti) - Math.cos(ti)
        f[i] = ui * ui + vi * vi
      }

      GSL::SUCCESS
    }

    brown_df = lambda { |x, t, y, df|
      brown_N.times { |i|
        ti = 0.2 * (i + 1)
        ui = x[0] + x[1] * ti - Math.exp(ti)
        vi = x[2] + x[3] * Math.sin(ti) - Math.cos(ti)

        df.set(i, 0, 2.0 * ui)
        df.set(i, 1, 2.0 * ui * ti)
        df.set(i, 2, 2.0 * vi)
        df.set(i, 3, 2.0 * vi * Math.sin(ti))
      }

      GSL::SUCCESS
    }

    fdf = GSL::MultiFit::Function_fdf.alloc(brown_f, brown_df, brown_P)
    fdf.set_data(GSL::Vector.alloc(brown_N), GSL::Vector.alloc(brown_N))

    _test_lmder(fdf, brown_x0, brown_X.vector_view, brown_F, brown_cov.vector_view)
  end

  def test_enso
    enso_N = 168
    enso_P = 9

    enso_x0 = GSL::Vector.alloc(10.0, 3.0, 0.5, 44.0, -1.5, 0.5, 26.0, 0.1, 1.5)

    enso_x = GSL::Vector.alloc(
      1.0510749193E+01,  3.0762128085E+00, 5.3280138227E-01,
      4.4311088700E+01, -1.6231428586E+00, 5.2554493756E-01,
      2.6887614440E+01,  2.1232288488E-01, 1.4966870418E+00
    )

    enso_sumsq = 7.8853978668E+02

    enso_sigma = GSL::Vector.alloc(
      1.7488832467E-01, 2.4310052139E-01, 2.4354686618E-01,
      9.4408025976E-01, 2.8078369611E-01, 4.8073701119E-01,
      4.1612939130E-01, 5.1460022911E-01, 2.5434468893E-01
    )

    enso_F = GSL::Vector.alloc(
      12.90000, 11.30000, 10.60000, 11.20000, 10.90000,  7.50000,  7.70000,
      11.70000, 12.90000, 14.30000, 10.90000, 13.70000, 17.10000, 14.00000,
      15.30000,  8.50000,  5.70000,  5.50000,  7.60000,  8.60000,  7.30000,
       7.60000, 12.70000, 11.00000, 12.70000, 12.90000, 13.00000, 10.90000,
      10.40000, 10.20000,  8.00000, 10.90000, 13.60000, 10.50000,  9.20000,
      12.40000, 12.70000, 13.30000, 10.10000,  7.80000,  4.80000,  3.00000,
       2.50000,  6.30000,  9.70000, 11.60000,  8.60000, 12.40000, 10.50000,
      13.30000, 10.40000,  8.10000,  3.70000, 10.70000,  5.10000, 10.40000,
      10.90000, 11.70000, 11.40000, 13.70000, 14.10000, 14.00000, 12.50000,
       6.30000,  9.60000, 11.70000,  5.00000, 10.80000, 12.70000, 10.80000,
      11.80000, 12.60000, 15.70000, 12.60000, 14.80000,  7.80000,  7.10000,
      11.20000,  8.10000,  6.40000,  5.20000, 12.00000, 10.20000, 12.70000,
      10.20000, 14.70000, 12.20000,  7.10000,  5.70000,  6.70000,  3.90000,
       8.50000,  8.30000, 10.80000, 16.70000, 12.60000, 12.50000, 12.50000,
       9.80000,  7.20000,  4.10000, 10.60000, 10.10000, 10.10000, 11.90000,
      13.60000, 16.30000, 17.60000, 15.50000, 16.00000, 15.20000, 11.20000,
      14.30000, 14.50000,  8.50000, 12.00000, 12.70000, 11.30000, 14.50000,
      15.10000, 10.40000, 11.50000, 13.40000,  7.50000,  0.60000,  0.30000,
       5.50000,  5.00000,  4.60000,  8.20000,  9.90000,  9.20000, 12.50000,
      10.90000,  9.90000,  8.90000,  7.60000,  9.50000,  8.40000, 10.70000,
      13.60000, 13.70000, 13.70000, 16.50000, 16.80000, 17.10000, 15.40000,
       9.50000,  6.10000, 10.10000,  9.30000,  5.30000, 11.20000, 16.60000,
      15.60000, 12.00000, 11.50000,  8.60000, 13.80000,  8.70000,  8.60000,
       8.60000,  8.70000, 12.80000, 13.20000, 14.00000, 13.40000, 14.80000
    )

    enso_f = lambda { |x, t, y, f|
      b = x

      enso_N.times { |i|
        ti, pi = t[i], GSL::M_PI

        yy  = b[0]
        yy += b[1] * Math.cos(2.0 * pi * ti / 12)
        yy += b[2] * Math.sin(2.0 * pi * ti / 12)
        yy += b[4] * Math.cos(2.0 * pi * ti / b[3])
        yy += b[5] * Math.sin(2.0 * pi * ti / b[3])
        yy += b[7] * Math.cos(2.0 * pi * ti / b[6])
        yy += b[8] * Math.sin(2.0 * pi * ti / b[6])

        f[i] = y[i] - yy
      }

      GSL::SUCCESS
    }

    enso_df = lambda { |x, t, y, df|
      b = x

      enso_N.times { |i|
        ti, pi = t[i], GSL::M_PI

        df.set(i, 0, -1.0)
        df.set(i, 1, -Math.cos(2.0 * pi * ti / 12))
        df.set(i, 2, -Math.sin(2.0 * pi * ti / 12))
        df.set(i, 3, -b[4] * (2.0 * pi * ti / (b[3] * b[3])) * Math.sin(2 * pi * ti / b[3]) + b[5] * (2 * pi * ti / (b[3] * b[3])) * Math.cos(2 * pi * ti / b[3]))
        df.set(i, 4, -Math.cos(2 * pi * ti / b[3]))
        df.set(i, 5, -Math.sin(2 * pi * ti / b[3]))
        df.set(i, 6, -b[7] * (2 * pi * ti / (b[6] * b[6])) * Math.sin(2 * pi * ti / b[6]) + b[8] * (2 * pi * ti / (b[6] * b[6])) * Math.cos(2 * pi * ti / b[6]))
        df.set(i, 7, -Math.cos(2 * pi * ti / b[6]))
        df.set(i, 8, -Math.sin(2 * pi * ti / b[6]))
      }

      GSL::SUCCESS
    }

    fdf = GSL::MultiFit::Function_fdf.alloc(enso_f, enso_df, enso_P)

    #fdf.set_data(GSL::Vector.alloc(1..168), enso_F)
    fdf.set_data(GSL::Vector.indgen(168, 1), enso_F)

    _test_fdf('nist-ENSO', fdf, enso_x0, enso_x, enso_sumsq, enso_sigma)
  end

  def test_filip
    filip_n = 82
    filip_p = 11

    filip_x = GSL::Vector.alloc(
      -6.860120914, -4.324130045, -4.358625055, -4.358426747, -6.955852379,
      -6.661145254, -6.355462942, -6.118102026, -7.115148017, -6.815308569,
      -6.519993057, -6.204119983, -5.853871964, -6.109523091, -5.79832982,
      -5.482672118, -5.171791386, -4.851705903, -4.517126416, -4.143573228,
      -3.709075441, -3.499489089, -6.300769497, -5.953504836, -5.642065153,
      -5.031376979, -4.680685696, -4.329846955, -3.928486195, -8.56735134,
      -8.363211311, -8.107682739, -7.823908741, -7.522878745, -7.218819279,
      -6.920818754, -6.628932138, -6.323946875, -5.991399828, -8.781464495,
      -8.663140179, -8.473531488, -8.247337057, -7.971428747, -7.676129393,
      -7.352812702, -7.072065318, -6.774174009, -6.478861916, -6.159517513,
      -6.835647144, -6.53165267,  -6.224098421, -5.910094889, -5.598599459,
      -5.290645224, -4.974284616, -4.64454848,  -4.290560426, -3.885055584,
      -3.408378962, -3.13200249,  -8.726767166, -8.66695597,  -8.511026475,
      -8.165388579, -7.886056648, -7.588043762, -7.283412422, -6.995678626,
      -6.691862621, -6.392544977, -6.067374056, -6.684029655, -6.378719832,
      -6.065855188, -5.752272167, -5.132414673, -4.811352704, -4.098269308,
      -3.66174277,  -3.2644011
    )

    filip_y = GSL::Vector.alloc(
      0.8116, 0.9072, 0.9052, 0.9039, 0.8053, 0.8377, 0.8667, 0.8809, 0.7975,
      0.8162, 0.8515, 0.8766, 0.8885, 0.8859, 0.8959, 0.8913, 0.8959, 0.8971,
      0.9021, 0.909,  0.9139, 0.9199, 0.8692, 0.8872, 0.89,   0.891,  0.8977,
      0.9035, 0.9078, 0.7675, 0.7705, 0.7713, 0.7736, 0.7775, 0.7841, 0.7971,
      0.8329, 0.8641, 0.8804, 0.7668, 0.7633, 0.7678, 0.7697, 0.77,   0.7749,
      0.7796, 0.7897, 0.8131, 0.8498, 0.8741, 0.8061, 0.846,  0.8751, 0.8856,
      0.8919, 0.8934, 0.894,  0.8957, 0.9047, 0.9129, 0.9209, 0.9219, 0.7739,
      0.7681, 0.7665, 0.7703, 0.7702, 0.7761, 0.7809, 0.7961, 0.8253, 0.8602,
      0.8809, 0.8301, 0.8664, 0.8834, 0.8898, 0.8964, 0.8963, 0.9074, 0.9119,
      0.9228
    )

    work = GSL::MultiFit::Workspace.alloc(filip_n, filip_p)

    expected_c = GSL::Vector.alloc(
      -1467.48961422980,      -2772.17959193342, -2316.37108160893,
      -1127.97394098372,       -354.478233703349,  -75.1242017393757,
        -10.8753180355343,       -1.06221498588947, -0.670191154593408e-01,
         -0.246781078275479e-02, -0.402962525080404e-04
    )

    expected_sd = GSL::Vector.alloc(
      298.084530995537,      559.779865474950,  466.477572127796,
      227.204274477751,       71.6478660875927,  15.2897178747400,
        2.23691159816033,      0.221624321934227, 0.142363763154724e-01,
        0.535617408889821e-03, 0.896632837373868e-05
    )

    expected_chisq = 0.795851382172941e-03

    xx = GSL::Matrix.alloc(filip_n, filip_p)

    filip_n.times { |i|
      filip_p.times { |j|
        xx.set(i, j, GSL.pow(filip_x[i], j))
      }
    }

    c, cov, chisq, _ = GSL::MultiFit.linear(xx, filip_y, work)

    assert_rel c[0],  expected_c[0],  1e-7, 'filip gsl_fit_multilinear c0'
    assert_rel c[1],  expected_c[1],  1e-7, 'filip gsl_fit_multilinear c1'
    assert_rel c[2],  expected_c[2],  1e-7, 'filip gsl_fit_multilinear c2'
    assert_rel c[3],  expected_c[3],  1e-7, 'filip gsl_fit_multilinear c3'
    assert_rel c[4],  expected_c[4],  1e-7, 'filip gsl_fit_multilinear c4'
    assert_rel c[5],  expected_c[5],  1e-7, 'filip gsl_fit_multilinear c5'
    assert_rel c[6],  expected_c[6],  1e-7, 'filip gsl_fit_multilinear c6'
    assert_rel c[7],  expected_c[7],  1e-7, 'filip gsl_fit_multilinear c7'
    assert_rel c[8],  expected_c[8],  1e-7, 'filip gsl_fit_multilinear c8'
    assert_rel c[9],  expected_c[9],  1e-7, 'filip gsl_fit_multilinear c9'
    assert_rel c[10], expected_c[10], 1e-7, 'filip gsl_fit_multilinear c10'

    diag = cov.diagonal

    assert_rel diag[0],  GSL.pow(expected_sd[0],2.0),  1e-6, 'filip gsl_fit_multilinear cov00'
    assert_rel diag[1],  GSL.pow(expected_sd[1],2.0),  1e-6, 'filip gsl_fit_multilinear cov11'
    assert_rel diag[2],  GSL.pow(expected_sd[2],2.0),  1e-6, 'filip gsl_fit_multilinear cov22'
    assert_rel diag[3],  GSL.pow(expected_sd[3],2.0),  1e-6, 'filip gsl_fit_multilinear cov33'
    assert_rel diag[4],  GSL.pow(expected_sd[4],2.0),  1e-6, 'filip gsl_fit_multilinear cov44'
    assert_rel diag[5],  GSL.pow(expected_sd[5],2.0),  1e-6, 'filip gsl_fit_multilinear cov55'
    assert_rel diag[6],  GSL.pow(expected_sd[6],2.0),  1e-6, 'filip gsl_fit_multilinear cov66'
    assert_rel diag[7],  GSL.pow(expected_sd[7],2.0),  1e-6, 'filip gsl_fit_multilinear cov77'
    assert_rel diag[8],  GSL.pow(expected_sd[8],2.0),  1e-6, 'filip gsl_fit_multilinear cov88'
    assert_rel diag[9],  GSL.pow(expected_sd[9],2.0),  1e-6, 'filip gsl_fit_multilinear cov99'
    assert_rel diag[10], GSL.pow(expected_sd[10],2.0), 1e-6, 'filip gsl_fit_multilinear cov1010'

    assert_rel chisq, expected_chisq, 1e-7, 'filip gsl_fit_multilinear chisq'

    expected_c = GSL::Vector.alloc(
      -1467.48961422980,      -2772.17959193342,      -2316.37108160893, -1127.97394098372,
       -354.478233703349,       -75.1242017393757,      -10.8753180355343,  -1.06221498588947,
         -0.670191154593408e-01, -0.246781078275479e-02, -0.402962525080404e-04
    )

    expected_cov = GSL::Matrix.alloc(
      [ 7.9269341767252183262588583867942e9,  1.4880416622254098343441063389706e10,
        1.2385811858111487905481427591107e10, 6.0210784406215266653697715794241e9,
        1.8936652526181982747116667336389e9,  4.0274900618493109653998118587093e8,
        5.8685468011819735806180092394606e7,  5.7873451475721689084330083708901e6,
        3.6982719848703747920663262917032e5,  1.3834818802741350637527054170891e4,
        2.301758578713219280719633494302e2 ],
      [ 1.4880416622254098334697515488559e10, 2.7955091668548290835529555438088e10,
        2.3286604504243362691678565997033e10, 1.132895006796272983689297219686e10,
        3.5657281653312473123348357644683e9,  7.5893300392314445528176646366087e8,
        1.1066654886143524811964131660002e8,  1.0921285448484575110763947787775e7,
        6.9838139975394769253353547606971e5,  2.6143091775349597218939272614126e4,
        4.3523386330348588614289505633539e2 ],
      [ 1.2385811858111487890788272968677e10, 2.3286604504243362677757802422747e10,
        1.9412787917766676553608636489674e10, 9.4516246492862131849077729250098e9,
        2.9771226694709917550143152097252e9,  6.3413035086730038062129508949859e8,
        9.2536164488309401636559552742339e7,  9.1386304643423333815338760248027e6,
        5.8479478338916429826337004060941e5,  2.1905933113294737443808429764554e4,
        3.6493161325305557266196635180155e2 ],
      [ 6.0210784406215266545770691532365e9,  1.1328950067962729823273441573365e10,
        9.4516246492862131792040001429636e9,  4.6053152992000107509329772255094e9,
        1.4517147860312147098138030287038e9,  3.0944988323328589376402579060072e8,
        4.5190223822292688669369522708712e7,  4.4660958693678497534529855690752e6,
        2.8599340736122198213681258676423e5,  1.0720394998549386596165641244705e4,
        1.7870937745661967319298031044424e2 ],
      [ 1.8936652526181982701620450132636e9,  3.5657281653312473058825073094524e9,
        2.9771226694709917514149924058297e9,  1.451714786031214708936087401632e9,
        4.5796563896564815123074920050827e8,  9.7693972414561515534525103622773e7,
        1.427717861635658545863942948444e7,   1.4120161287735817621354292900338e6,
        9.0484361228623960006818614875557e4,  3.394106783764852373199087455398e3,
        5.6617406468519495376287407526295e1 ],
      [ 4.0274900618493109532650887473599e8,  7.589330039231444534478894935778e8,
        6.3413035086730037947153564986653e8,  3.09449883233285893390542947998e8,
        9.7693972414561515475770399055121e7,  2.0855726248311948992114244257719e7,
        3.0501263034740400533872858749566e6,  3.0187475839310308153394428784224e5,
        1.9358204633534233524477930175632e4,  7.2662989867560017077361942813911e2,
        1.2129002231061036467607394277965e1 ],
      [ 5.868546801181973559370854830868e7,   1.1066654886143524778548044386795e8,
        9.2536164488309401413296494869777e7,  4.5190223822292688587853853162072e7,
        1.4277178616356585441556046753562e7,  3.050126303474040051574715592746e6,
        4.4639982579046340884744460329946e5,  4.4212093985989836047285007760238e4,
        2.8371395028774486687625333589972e3,  1.0656694507620102300567296504381e2,
        1.7799982046359973175080475654123e0 ],
      [ 5.7873451475721688839974153925406e6,  1.0921285448484575071271480643397e7,
        9.1386304643423333540728480344578e6,  4.4660958693678497427674903565664e6,
        1.4120161287735817596182229182587e6,  3.0187475839310308117812257613082e5,
        4.4212093985989836021482392757677e4,  4.3818874017028389517560906916315e3,
        2.813828775753142855163154605027e2,   1.0576188138416671883232607188969e1,
        1.7676976288918295012452853715408e-1 ],
      [ 3.6982719848703747742568351456818e5,  6.9838139975394768959780068745979e5,
        5.8479478338916429616547638954781e5,  2.8599340736122198128717796825489e5,
        9.0484361228623959793493985226792e4,  1.9358204633534233490579641064343e4,
        2.8371395028774486654873647731797e3,  2.8138287757531428535592907878017e2,
        1.8081118503579798222896804627964e1,  6.8005074291434681866415478598732e-1,
        1.1373581557749643543869665860719e-2 ],
      [ 1.3834818802741350562839757244708e4,  2.614309177534959709397445440919e4,
        2.1905933113294737352721470167247e4,  1.0720394998549386558251721913182e4,
        3.3941067837648523632905604575131e3,  7.2662989867560016909534954790835e2,
        1.0656694507620102282337905013451e2,  1.0576188138416671871337685672492e1,
        6.8005074291434681828743281967838e-1, 2.5593857187900736057022477529078e-2,
        4.2831487599116264442963102045936e-4 ],
      [ 2.3017585787132192669801658674163e2,  4.3523386330348588381716460685124e2,
        3.6493161325305557094116270974735e2,  1.7870937745661967246233792737255e2,
        5.6617406468519495180024059284629e1,  1.2129002231061036433003571679329e1,
        1.7799982046359973135014027410646e0,  1.7676976288918294983059118597214e-1,
        1.137358155774964353146460100337e-2,  4.283148759911626442000316269063e-4,
        7.172253875245080423800933453952e-6 ]
    )

    expected_chisq = 0.795851382172941E-03

    filip_n.times { |i|
      filip_p.times { |j|
        xx.set(i, j, GSL.pow(filip_x[i], j))
      }
    }

    w = GSL::Vector.alloc(filip_n)
    w.set_all(1.0)

    c, cov, _, _ = GSL::MultiFit.wlinear(xx, w, filip_y, work)

    filip_p.times { |i|
      assert_rel c[i], expected_c[i], 1e-7, "filip gsl_fit_multilinear c#{i}"
    }

    filip_p.times { |i|
      filip_p.times { |j|
        assert_rel cov[i, j], expected_cov[i, j], 1e-6, "filip gsl_fit_wmultilinear cov(#{i},#{j})"
      }
    }
  end

  def test_gauss
    maxiter = 10
    n = 1000

    # model: a * exp(-(x - x0) ** 2 / 2 / sigma ** 2)
    gauss_p = 3
    gauss_f = lambda { |x, t, y, s, f|
      # x: parameters as a Vecor
      # t: observed points as a GSL::Vector
      # y: observed data as a GSL::Vector
      # s: errorbar
      # f: result
      a = x[0]
      x0 = x[1]
      sigma2 = x[2] ** 2

      y.size.times { |i|
        f.set(i, (a * Math.exp(-(t[i] - x0) ** 2 / 2 / sigma2) - y[i]) / s[i])
      }

      GSL::SUCCESS
    }

    gauss_df = lambda { |x, t, y, s, df|
      a = x[0]
      x0 = x[1]
      sigma = x[2]
      sigma2 = sigma ** 2

      y.size.times { |i|
        dx = t[i] - x0
        dx2 = dx ** 2
        f = a * Math.exp(-dx2 / 2 / sigma2)

        df.set(i, 0, f / a / s[i])
        df.set(i, 1, f * dx / sigma2 / s[i])
        df.set(i, 2, f * dx2 / sigma2 / sigma / s[i])
      }

      GSL::SUCCESS
    }

    # goal
    xgoal = GSL::Vector.alloc([1, 0, 1])
    parname = %w[a x0 si]

    # data
    t = GSL::Vector.alloc(n)  # positions of data
    tmin = -10.0
    tmax =  10.0

    t.size.times { |i|
      t[i] = tmin + (tmax - tmin) * i / (n - 1)
    }

    stdev = xgoal[0] * 0.1

    s = GSL::Vector.alloc(Array.new(t.size, stdev))  # error bar of each datum
    r = GSL::Rng.alloc
    e = GSL::Vector.alloc(t.size)

    t.size.times { |i|
      e[i] = -r.gaussian(stdev)  # perturbation to data
    }

    y = GSL::Vector.alloc(t.size)
    n = GSL::Vector.alloc(Array.new(t.size, 1.0))
    gauss_f.call(xgoal, t, e, n, y)  # data: y = model - e

    # fitting
    x = GSL::Vector.alloc([0.5, 0.1, 2])  # initial guess

    fdf = GSL::MultiFit::Function_fdf.alloc(gauss_f, gauss_df, gauss_p)
    fdf.set_data(t, y, s)

    solver = GSL::MultiFit::FdfSolver.alloc(GSL::MultiFit::FdfSolver::LMSDER, t.size, gauss_p)
    solver.set(fdf, x)

    #solver.print_state(0)

    maxiter.times { |i|
      solver.iterate

      status = solver.test_delta(1e-6, 1e-6)
      #solver.print_state(i + 1)

      break if status != GSL::CONTINUE
    }

    # results
    covar = solver.covar(0.0)
    xresult = solver.position
    dof = t.size - gauss_p
    chi2 = GSL.pow_2(solver.f.dnrm2)
    xsigma = GSL::Vector.alloc(xresult.size)

    xresult.size.times { |i|
      xsigma[i] = Math.sqrt(chi2 / dof * covar[i, i]) * 2.0
      # resulting parameters to differ two times than standard error
    }

    desc = 'a*exp(-(x-x0)**2/2/si**2), chi2/N:%.3g' % (chi2 / t.size)

    xresult.size.times { |i|
      assert_rel xresult[i], xgoal[i], xsigma[i], '%s: %-2.2s' % [desc, parname[i]]
      refute((xresult[i] - xgoal[i]).abs > xsigma[i], desc)
    }
  end

  def test_longley
    longley_n = 16
    longley_p = 7

    longley_x = GSL::Vector.alloc(
      1,  83.0, 234289, 2356, 1590, 107608, 1947,
      1,  88.5, 259426, 2325, 1456, 108632, 1948,
      1,  88.2, 258054, 3682, 1616, 109773, 1949,
      1,  89.5, 284599, 3351, 1650, 110929, 1950,
      1,  96.2, 328975, 2099, 3099, 112075, 1951,
      1,  98.1, 346999, 1932, 3594, 113270, 1952,
      1,  99.0, 365385, 1870, 3547, 115094, 1953,
      1, 100.0, 363112, 3578, 3350, 116219, 1954,
      1, 101.2, 397469, 2904, 3048, 117388, 1955,
      1, 104.6, 419180, 2822, 2857, 118734, 1956,
      1, 108.4, 442769, 2936, 2798, 120445, 1957,
      1, 110.8, 444546, 4681, 2637, 121950, 1958,
      1, 112.6, 482704, 3813, 2552, 123366, 1959,
      1, 114.2, 502601, 3931, 2514, 125368, 1960,
      1, 115.7, 518173, 4806, 2572, 127852, 1961,
      1, 116.9, 554894, 4007, 2827, 130081, 1962
    )

    longley_y = GSL::Vector.alloc(
      60323, 61122, 60171, 61187, 63221, 63639, 64989, 63761,
      66019, 67857, 68169, 66513, 68655, 69564, 69331, 70551
    )

    work = GSL::MultiFit::Workspace.alloc(longley_n, longley_p)

    x = GSL::Matrix.alloc(longley_x, longley_n, longley_p).view
    y = longley_y.view

    expected_c = GSL::Vector.alloc(
      -3482258.63459582,
            15.0618722713733,
            -0.358191792925910e-01,
            -2.02022980381683,
            -1.03322686717359,
            -0.511041056535807e-01,
          1829.15146461355
    )

    expected_sd = GSL::Vector.alloc(
      890420.383607373,
          84.9149257747669,
           0.334910077722432e-01,
           0.488399681651699,
           0.214274163161675,
           0.226073200069370,
         455.478499142212
    )

    expected_chisq = 836424.055505915

    c, cov, chisq, _ = GSL::MultiFit.linear(x, y, work)

    7.times { |i|
      assert_rel c[i], expected_c[i], 1e-10, "longley gsl_fit_multilinear c#{i}"
    }

    diag = cov.diagonal

    assert_rel diag[0], GSL.pow(expected_sd[0],2.0), 1e-10, 'longley gsl_fit_multilinear cov00'
    assert_rel diag[1], GSL.pow(expected_sd[1],2.0), 1e-10, 'longley gsl_fit_multilinear cov11'
    assert_rel diag[2], GSL.pow(expected_sd[2],2.0), 1e-10, 'longley gsl_fit_multilinear cov22'
    assert_rel diag[3], GSL.pow(expected_sd[3],2.0), 1e-10, 'longley gsl_fit_multilinear cov33'
    assert_rel diag[4], GSL.pow(expected_sd[4],2.0), 1e-10, 'longley gsl_fit_multilinear cov44'
    assert_rel diag[5], GSL.pow(expected_sd[5],2.0), 1e-10, 'longley gsl_fit_multilinear cov55'
    assert_rel diag[6], GSL.pow(expected_sd[6],2.0), 1e-10, 'longley gsl_fit_multilinear cov66'

    assert_rel chisq, expected_chisq, 1e-10, 'longley gsl_fit_multilinear chisq'

    expected_cov = GSL::Matrix.alloc(
      [ 8531122.56783558, -166.727799925578, 0.261873708176346, 3.91188317230983,
        1.1285582054705, -0.889550869422687, -4362.58709870581 ],
      [ -166.727799925578, 0.0775861253030891, -1.98725210399982e-05, -0.000247667096727256,
        -6.82911920718824e-05, 0.000136160797527761, 0.0775255245956248 ],
      [ 0.261873708176346, -1.98725210399982e-05, 1.20690316701888e-08, 1.66429546772984e-07,
        3.61843600487847e-08, -6.78805814483582e-08, -0.00013158719037715 ],
      [ 3.91188317230983, -0.000247667096727256, 1.66429546772984e-07, 2.56665052544717e-06,
        6.96541409215597e-07, -9.00858307771567e-07, -0.00197260370663974 ],
      [ 1.1285582054705, -6.82911920718824e-05, 3.61843600487847e-08, 6.96541409215597e-07,
        4.94032602583969e-07, -9.8469143760973e-08, -0.000576921112208274 ],
      [ -0.889550869422687, 0.000136160797527761, -6.78805814483582e-08, -9.00858307771567e-07,
        -9.8469143760973e-08, 5.49938542664952e-07, 0.000430074434198215 ],
      [ -4362.58709870581, 0.0775255245956248, -0.00013158719037715, -0.00197260370663974,
        -0.000576921112208274, 0.000430074434198215, 2.23229587481535 ]
    )

    expected_chisq = 836424.055505915

    w = GSL::Vector.alloc(longley_n)
    w.set_all(1.0)

    c, cov, chisq, _ = GSL::MultiFit.wlinear(x, w, y, work)

    7.times { |i|
      assert_rel c[i], expected_c[i], 1e-10, "longley gsl_fit_wmultilinear c#{i}"
    }

    longley_p.times { |i|
      longley_p.times { |j|
        assert_rel cov[i, j], expected_cov[i, j], 1e-7, "longley gsl_fit_wmultilinear cov(#{i},#{j})"
      }
    }

    assert_rel chisq, expected_chisq, 1e-10, 'longley gsl_fit_wmultilinear chisq'
  end

end
