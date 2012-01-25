require File.join(File.dirname(__FILE__), "nmatrix/nmatrix.so")

class NMatrix
  VERSION = '0.0.1'

  #def inspect
  #
  #end

  # TODO: Make this actually pretty.
  def pretty_print
    raise(NotImplementedError, "can only print rank 2 matrices") unless rank == 2
    (0...shape[0]).each do |i|
      arr = []
      (0...shape[1]).each do |j|
        arr << self[i,j]
      end
      puts arr.join("  ")
    end
    nil
  end


  def inspect
    original_inspect = super
    original_inspect = original_inspect[0...original_inspect.size-1]
    ary = [original_inspect]
    ary << "shape:#{shape}" << dtype << stype
    ary.join(" ") + ">"
  end


end
