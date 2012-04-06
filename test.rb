require "./lib/nmatrix.rb"

    original = NMatrix.new(:dense, [3,3], [0,0,1,0,2,0,3,4,5], :int64).
        scast(:yale, :int32).
        scast(:dense, :float64).
        scast(:list, :int32).
        scast(:dense, :int16).
        scast(:list, :int32).
        scast(:yale, :int64)

