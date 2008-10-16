#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "combination/test.c"
require("gsl")
require("./gsl_test.rb")
include GSL::Test

c63 = GSL::Matrix.alloc([0, 1, 2], [0, 1, 3], [0, 1, 4], [0, 1, 5], 
                 [0, 2, 3], [0, 2, 4], [0, 2, 5], 
                 [0, 3, 4], [0, 3, 5],
                 [0, 4, 5], 
                 [1, 2, 3], [1, 2, 4], [1, 2, 5],
                 [1, 3, 4], [1, 3, 5],
                 [1, 4, 5],
                 [2, 3, 4], [2, 3, 5], [2, 4, 5],
                 [3, 4, 5])

GSL::IEEE.env_setup()

status = false
c = GSL::Combination.alloc(6, 3)
c.init_first

i = 0
begin
  if i >= 20
    status = true
    break
  end
  for j in 0...3
    status |= (c.data[j] != c63[i][j])
  end

  s1 = c.valid?
  desc = sprintf("GSL::Combination#valid\(%u\)", i)
  GSL::Test.test(s1, desc)
  i += 1
end while c.next == GSL::SUCCESS

GSL::Test.test(status, "GSL::Combination#next, 6 choose 3 combination, 20 steps")

c.next
c.next
c.next

for j in 0...3 
  status |= (c.data[j] != c63[19][j])
end
GSL::Test.test(status, "GSL::Combination#next on the last combination")

s1 = c.valid?
GSL::Test.test(s1, "GSL::Combination#valid on the last combination")

d = GSL::Combination.alloc(6, 3)
GSL::Combination.memcpy(d, c)

status = false
for j in 0...3 
  status |= (d.data[j] != c.data[j])
end
GSL::Test.test(status, "GSL::Combination.memcpy, 6 choose 3 combination")

c.init_last
i = 20
begin
  if i == 0
    status = true
    break;
  end
  i -= 1
  for j in 0...3
    status |= (c.data[j] != c63[i][j])
  end
  s1 = c.valid?
  desc = sprintf("GSL::Combination#valid\(%u\)", i)
  GSL::Test.test(s1, desc)
end while c.prev == GSL::SUCCESS

GSL::Test.test(status, "GSL::Combination#prev, 6 choose 3 combination, 20 steps")

c.prev
c.prev
c.prev

for j in 0...3
  status |= (c.data[j] != c63[0][j])
end
GSL::Test.test(status, "GSL::Combination#prev on the first combination")

s1 = c.valid?
GSL::Test.test(s1, "GSL::Combination#valid on the first combination")
d = GSL::Combination.alloc(6, 3)
GSL::Combination.memcpy(d, c)

status = false
for j in 0...3 
  status |= (d.data[j] != c.data[j])
end
GSL::Test.test(status, "GSL::Combination.memcpy, 6 choose 3 combination")

c = GSL::Combination.calloc(7, 0)
status |= (c.next != GSL::FAILURE)
status |= (c.next != GSL::FAILURE)
status |= (c.prev != GSL::FAILURE)
status |= (c.prev != GSL::FAILURE)
GSL::Test.test(status, "GSL::Combination 7 choose 0")

c = GSL::Combination.calloc(7, 7)
for j in 0...7
  status |= (c.get(j) != j)
end
status |= (c.next != GSL::FAILURE)
for j in 0...7
  status |= (c.get(j) != j)
end
status |= (c.next != GSL::FAILURE)
for j in 0...7
  status |= (c.get(j) != j)
end
status |= (c.next != GSL::FAILURE)
for j in 0...7
  status |= (c.get(j) != j)
end
GSL::Test.test(status, "GSL::Combination 7 choose 7")
