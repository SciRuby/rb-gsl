require "./lib/nmatrix.rb"

puts n = NMatrix.new(:yale, [2,3], :float64)

puts n[0,0] = 0.01
puts n[1,1] = 0.1
puts n[0,1] = 0.2
puts n[1,0] = 0.3
puts n[1,2] = 0.4

puts n[0,0] == 0.01
puts n[1,1] == 0.1
puts n[0,1] == 0.2
puts n[1,0] == 0.3
puts n[1,2] == 0.4

#puts m = n.dup
#m.pretty_print

#puts m[0,0] = -0.1
#puts m[0,1] = -0.2
#puts m[1,0] = -1.0
#m.pretty_print


#puts m.shape == n.shape
#puts m.rank == n.rank
#puts m.object_id == n.object_id
#puts m.stype
#puts m.dtype
#m.pretty_print
