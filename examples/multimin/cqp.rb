#!/usr/bin/env ruby
require("gsl")

module GSL::CQP
  class Test_Problem
    def Test_Problem.gould()
      cqp_data = GSL::CQP::Data.alloc()
      cqp_data.Q = GSL::Matrix.eye(2, 2)
      cqp_data.q = GSL::Vector[-2.0, -1.0]
      cqp_data.A = GSL::Matrix[[3, 1], 1, 2]
      cqp_data.b = GSL::Vector.alloc(1); cqp_data.b[0] = 1.5
      cqp_data.C = GSL::Matrix[[1, 0, 0, 1], 2, 2]
      cqp_data.d = GSL::Vector.calloc(2)

      Test_Problem.new("Goulds's problem", 2, 1, 2, cqp_data, 0.0)
    end
    def Test_Problem.betts()
      cqp_data = GSL::CQP::Data.alloc()
      cqp_data.Q = GSL::Matrix[[0.02, 0, 0, 2], 2, 2]
      cqp_data.q = GSL::Vector.calloc(2)
      cqp_data.A = GSL::Matrix[[10, -1.0], 1, 2]
      cqp_data.b = GSL::Vector.alloc(1); cqp_data.b[0] = 20
      cqp_data.C = GSL::Matrix.calloc(4, 2)
      cqp_data.C[0,0] = 1.0; cqp_data.C[1,1] = 1.0
      cqp_data.C[2,0] = -1.0; cqp_data.C[3,1] = -11.0
      cqp_data.d = GSL::Vector[2.0, -50, -50, -50]


      Test_Problem.new("Betts's problem", 2, 1, 4, cqp_data, 0.04)
    end
    def Test_Problem.beale()
      cqp_data = GSL::CQP::Data.alloc()
      cqp_data.Q = GSL::Matrix[[4, 2, 2, 2, 4, 0, 2, 0, 2], 3, 3]
      cqp_data.q = GSL::Vector[-8, -6, -4]
      cqp_data.A = GSL::Matrix[[-1, -1, -2], 1, 3]
      cqp_data.b = GSL::Vector.alloc(1); cqp_data.b[0] = -3.0
      cqp_data.C = GSL::Matrix.eye(3, 3)
      cqp_data.d = GSL::Vector.calloc(3)

      cqp_data.Q[0,0] = 4.0; cqp_data.Q[0,1] = 2.0
      cqp_data.Q[1,0] = 2.0; cqp_data.Q[1,1] = 4.0
      cqp_data.Q[2,0] = 2.0; cqp_data.Q[2,2] = 2.0
      Test_Problem.new("Beale's problem", 3, 1, 3, cqp_data, 9.0+1.0/9.0)
    end

    def initialize(name, n, me, mi, cqp, opt_value)
      @name = name
      @n = n
      @me = me
      @mi = mi
      @cqp = cqp
      @opt_value = opt_value
    end

    def solve()
      max_iter = 1000
      iter = 1
      status = GSL::CONTINUE
      s = GSL::CQP::Minimizer.alloc("mg_pdip", @n, @me, @mi)
      s.set(@cqp)
      printf("********************  %s  ********************\n\n", @name)

      printf("== Itn ======= f ======== ||gap|| ==== ||residual||\n\n")

      begin
        status = s.iterate
        status = s.test_convergence(1e-10, 1e-10)
        printf("%4d   %14.8f  %13.6e  %13.6e\n", iter, s.f, s.gap, s.residuals_norm)
        if status == GSL::SUCCESS
          printf("\nMinimum is found at\n");
          x = s.x
          lm_eq = s.lm_eq
          lm_ineq = s.lm_ineq
          for j in 0...x.size do
            printf("%9.6f ", x[j])
          end
          printf("\n\n")
          printf("\nLagrange-multipliers for Ax=b\n")
          for j in 0...lm_eq.size do
            printf("%9.6f ", lm_eq[j])
          end
          printf("\n\n")
          printf("\nLagrange-multipliers for Cx>=d\n");
          for j in 0...lm_ineq.size do
            printf("%9.6f ", lm_ineq[j])
          end
          printf("\n\n")
        else
          iter += 1
        end
     end while status == GSL::CONTINUE and iter <= max_iter
      GSL::SUCCESS
    end

    attr_accessor :name, :n, :me, :mi
    attr_accessor :cqp, :opt_value
  end

end

tp = Array.new(3)
tp[0] = GSL::CQP::Test_Problem.gould()
tp[1] = GSL::CQP::Test_Problem.betts()
tp[2] = GSL::CQP::Test_Problem.beale()

tp.each do |t|
  t.solve()
end

