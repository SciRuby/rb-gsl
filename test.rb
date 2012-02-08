require "./lib/nmatrix.rb"

a = NMatrix.new(:yale, 4, :float32)
a[0,1] = 4.0
a[1,2] = 1.0
a[1,3] = 1.0
a[3,1] = 2.0
a.__yale_print__
a.pretty_print

b = a.dup
c = a.multiply(b)
