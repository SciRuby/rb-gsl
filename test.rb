require "./lib/nmatrix.rb"

n = NMatrix.new(:yale, [10,300], :float64)
n.__yale_print__
STDERR.puts n.__yale_ija__
STDERR.puts n.__yale_a__
n[5,1]   = 1.0
puts n[5,1] == 1.0
n[5,0]   = 1.5
puts n[5,0] == 1.5
n[5,15] = 2.0
puts n[5,15] == 2.0
n[5,291] = 3.0
n[5,292] = 4.0
n[5,289] = 5.0
n[5,290] = 6.0
n[5,293] = 2.0
n[5,299] = 7.0
n[5,100] = 8.0

puts n[5,290] == 2.0
puts n[5,291] == 3.0
puts n[5,292] == 4.0
puts n[5,289] == 5.0
puts n[5,290] == 6.0
puts n[5,293] == 2.0
puts n[5,299] == 7.0
puts n[5,100] == 8.0
puts n[5,0] == 1.5
puts n[5,1] == 1.0
