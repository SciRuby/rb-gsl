#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create test Vector with 15 elements initialized to 0.0
v = GSL::Vector::alloc(15)

# Set elements of v from file smpv.dat
File.open("smpv.dat") do |f|
  v.fscanf(f)
end

# Show v
v

# Create another test vector with three elements
# initialized to 10.0, 1.0, and 3.0.
v2 = GSL::Vector.alloc(10, 1, 3)
END

__END__

# TODO Remove examples below that are duplicated elsewhere
# TODO Include examples below that are not duplicated elsewhere
v = GSL::Vector::alloc([9, 1, 2, 3, 12, 6, 0.1, 0.56, 5, 7, 2])
v.print
p v.size
v2 = GSL::Vector::alloc(11)
v2.set(1, 123)
p v2.to_a

v = GSL::Vector::alloc([1, 2, 3, 4, 5])
v.print
a = v.to_a
p a
a[2] = 12.0
v2 = a.to_gv
v2.print
__END__

p = v.sort_index
p p.to_a

v2 = GSL::Vector::alloc([5, 6, 7])

v3 = v2.scale!(2)
v3.print
v2.print

p v3.minmax

p v3.max_index
p v3.min_index
p v3.minmax_index

a = [1, 2, 3]
v = a.to_gv
p v.to_a

__END__

v3 = v * v2
v3.print

v.print

v.mul!(v2)
v.print

__END__

v = GSL::Vector::alloc(10)
for i in 1...10 do
  v.set(i, i.to_f)
end
v.print

v.swap_elements(3, 5)
v.print

v.reverse.print

__END__

v2 = v.subvector(2, 3)
v2.print

v3 = v.subvector_with_stride(2, 2, 3)
v3.print

v3.set([2, 3, 9])
v3.print

v4 = GSL::Vector::alloc([1, 2, 3, 5])
v4.print

a = v4.to_a
p a

__END__
v.each do |x|
  p x
end
