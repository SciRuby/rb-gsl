# nmatrix_spec.rb
require "./lib/nmatrix"

describe NVector do
  it "correctly initializes" do
    v = NVector.new 5, :float64
    v.shape[0].should == 5
    v.shape[1].should == 1
  end

  it "permits setting and getting contents" do
    v = NVector.new 5, :float64
    v[0] = 1.555
    v[0].should == 1.555
  end

end