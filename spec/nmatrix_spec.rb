# nmatrix_spec.rb
require "./src/nmatrix.so"

describe NMatrix do
  [:dense, :list].each do |storage_type|
    context storage_type do
      it "enforces shape boundaries" do
        lambda { NMatrix.new(storage_type, [1,10], 0)[-1,0] }.should raise_error(ArgumentError, "out of range")
        lambda { NMatrix.new(storage_type, [1,10], 0)[0,0]  }.should_not raise_error
        lambda { NMatrix.new(storage_type, [1,10], 0)[1,0]  }.should raise_error(ArgumentError, "out of range")
        lambda { NMatrix.new(storage_type, [1,10], 0)[0,9]  }.should_not raise_error
        lambda { NMatrix.new(storage_type, [1,10], 0)[0,10] }.should raise_error(ArgumentError, "out of range")
      end

      it "correctly gets default value" do
        NMatrix.new(storage_type, 3, 0)[1,1].should   == 0
        NMatrix.new(storage_type, 3, 0.1)[1,1].should == 0.1
        NMatrix.new(storage_type, 3, 1)[1,1].should   == 1
      end

      it "correctly sets and gets" do
        n = NMatrix.new(storage_type, 2, 0)
        (n[0,1] = 1).should == 1
        n[0,0].should == 0
        n[1,0].should == 0
        n[0,1].should == 1
        n[1,1].should == 0
      end
    end
  end

  it "correctly handles dense construction" do
    NMatrix.new(3,0)[1,1].should == 0
    NMatrix.new(3,:int8)[1,1].should == 0
  end

end