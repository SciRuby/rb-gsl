#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create two test Vectors
a = GSL::Vector.alloc(1, 2, 3, 4, 5)
b = GSL::Vector.alloc(6, 7, 8, 9, 10)

# Show element-wise multiplication
a * b

# Show that *= is not an in-place operation
# (it sets receiver to new object).
a.object_id
a *= b
a.object_id

# Show that #mul! is a true in-place operation
# (object_id does not change).
a.object_id
a.mul!(b)
a.object_id

# Same as "a = a / b" (element-wise division)
a /= b

# Show multiply by constant
a * 2

# Show coersion
2 * a

# Show coersion
5 / a

# Show that Vector * Vector::Col produces dot product
a * b.trans
a.dot(b)

# Show that Vector::Col * Vector produces Matrix
a.trans * b
END
