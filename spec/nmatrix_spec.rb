# nmatrix_spec.rb
require "./src/nmatrix.so"

describe NMatrix do

  it "list correctly handles missing initialization value" do
    NMatrix.new(:list, 3, :int8)[0,0].should    == 0
    NMatrix.new(:list, 4, :float64)[0,0].should == 0.0
  end


  it "dense correctly handles missing initialization value" do
    NMatrix.new(3, :int8)[0,0]
    NMatrix.new(4, :float64)[0,0]
  end


  [:dense, :list].each do |storage_type|
    context "(storage: #{storage_type})" do
      it "can be duplicated" do
        n = NMatrix.new(storage_type, [2,3], 1.1)
        m = n.dup
        m.shape.should == n.shape
        m.rank.should == n.rank
        m.object_id.should_not == n.object_id
        m[0,0].should == n[0,0]
        m[0,0] = 3.0
        m[0,0].should_not == n[0,0]
      end

      it "enforces shape boundaries" do
        lambda { NMatrix.new(storage_type, [1,10], 0)[-1,0] }.should raise_error
        lambda { NMatrix.new(storage_type, [1,10], 0)[1,0]  }.should raise_error(ArgumentError, "out of range")
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

      it "returns correct shape and rank" do
        NMatrix.new(storage_type, [3,2,8], 0).shape.should == [3,2,8]
        NMatrix.new(storage_type, [3,2,8], 0).rank.should  == 3
      end
    end
  end

  it "correctly handles dense construction" do
    NMatrix.new(3,0)[1,1].should == 0
    lambda { NMatrix.new(3,:int8)[1,1] }.should_not raise_error
  end

end