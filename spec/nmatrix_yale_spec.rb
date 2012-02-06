require "./lib/nmatrix"

describe NMatrix do
  it "correctly sets diagonal values" do
    n = NMatrix.new(:yale, [2,3], :float64)
    n[1,1] = 0.1
    n[0,0] = 0.2
    n.__yale_d__.should == [0.2, 0.1]
  end

  it "does not resize until necessary" do
    n = NMatrix.new(:yale, [2,3], :float64)
    n.__yale_size__.should == 3
    n.capacity.should == 5
    n[0,0] = 0.1
    n[0,1] = 0.2
    n[1,0] = 0.3
    n.__yale_size__.should == 5
    n.capacity.should == 5
  end


  it "correctly sets when not resizing" do
    n = NMatrix.new(:yale, [2,3], :float64)
    n[0,0] = 0.1
    n[0,1] = 0.2
    n[1,0] = 0.3
    n.__yale_a__ == [0.1, 0.0, 0.0, 0.2, 0.3]
    n.__yale_ija__ == [3,4,5,1,0]
  end

  it "correctly sets when resizing" do
    n = NMatrix.new(:yale, [2,3], :float64)
    n[0,0] = 0.01
    n[1,1] = 0.1
    n[0,1] = 0.2
    n[1,0] = 0.3
    n[1,2] = 0.4
    n.__yale_d__.should == [0.01, 0.1]
    n.__yale_ia__.should == [3,4,6]
    n.__yale_ja__.should == [1,0,2,nil]
    n.__yale_lu__.should == [0.2, 0.3, 0.4, nil]
  end

  it "correctly sets values within large rows" do
    n = NMatrix.new(:yale, [10,300], :float64)
    n[5,1]   = 1.0
    n[5,1].should == 1.0
    n[5,0]   = 1.5
    n[5,0].should == 1.5
    n[5,15] = 2.0
    n[5,15].should == 2.0
    n[5,291] = 3.0
    n[5,292] = 4.0
    n[5,289] = 5.0
    n[5,290] = 6.0
    n[5,293] = 2.0
    n[5,299] = 7.0
    n[5,100] = 8.0

    n[5,290].should == 2.0
    n[5,291].should == 3.0
    n[5,292].should == 4.0
    n[5,289].should == 5.0
    n[5,290].should == 6.0
    n[5,293].should == 2.0
    n[5,299].should == 7.0
    n[5,100].should == 8.0
    n[5,0].should == 1.5
    n[5,1].should == 1.0
  end

end