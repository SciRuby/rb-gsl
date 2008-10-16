require("gsl")

# This script is a copy of the one found in the NArray package, 
# written by M. Tanaka.
#   http://www.ir.isas.ac.jp/~masa/ruby/index-e.html
module GSL::Bench

  def bench_time(n=REPEAT)
#    t1 = T.times.utime
    t1 = Time.new
    n.times do
      yield
    end
#    t2 = T.times.utime 
    t2 = Time.new
    puts " Time: %.2f sec\n\n" % [t2-t1]
  end
end
