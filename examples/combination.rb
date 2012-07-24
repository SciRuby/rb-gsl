#!/usr/bin/env ruby
require("gsl")

c = 0
printf("All subsets of {0,1,2,3} by size:\n") ;
for i in 0...4 do
  c = GSL::Combination.calloc(4, i);
  begin
    printf("{");
    c.fprintf(STDOUT, " %u");
    printf(" }\n");
  end while c.next == GSL::SUCCESS
end

p c
c2 = c.clone
p c2

a = c.data
p a.class
p a
p a.print

