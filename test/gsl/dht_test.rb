require 'test_helper'

class DhtTest < GSL::TestCase

  N = 128

  def test_dht
    test_data = [GSL::Vector.alloc(N)]
    test_data << NMatrix.new([N], dtype: :float64) if ENV['NMATRIX']
    test_data.each do |vin|
      dht2 vin
      dht3 vin
      dht4 vin
    end
  end
  
  def test_dht1
    test_data = [GSL::Vector.alloc(1, 2, 3)]
    test_data << NMatrix.new([3], [1,2,3], dtype: :float64) if ENV['NMATRIX']
    test_data.each do |vin|
      dht = GSL::Dht.alloc(3, 1.0, 1.0)

      vout = dht.apply(vin)

      assert_in_delta  0.3752546494075203,  vout[0], 0.001
      assert_in_delta(-0.13350787269556064, vout[1], 0.001)
      assert_in_delta  0.0446799251438404,  vout[2], 0.001

      vin2 = dht.apply(vout)
      vin2 = vin2 * (13.323691936314223 ** 2)

      assert_in_delta 1.0000119186762644, vin2[0], 0.001
      assert_in_delta 1.9999790476647084, vin2[1], 0.001
      assert_in_delta 3.000035803234503,  vin2[2], 0.001
    end
  end

  def dht2 vin
    dht = GSL::Dht.alloc(N, 0.0, 100.0)

    N.times { |i|
      x = dht.x_sample(i)
      vin[i] = 1.0 / (1.0 + x * x)
    }

    vout = dht.apply(vin)

    assert_in_delta 3.999613382195876,   vout[0], 0.001
    assert_in_delta 1.8387637474026606,  vout[5], 0.001
    assert_in_delta 1.2677885358829588,  vout[10], 0.001
    assert_in_delta 0.3521910403797792,  vout[35], 0.001
    assert_in_delta 0.02373661279695407, vout[100], 0.001
  end

  def dht3 vin
    dht = GSL::Dht.alloc(N, 1.0, 20.0)

    N.times { |i|
      x = dht.x_sample(i)
      vin[i] = Math.exp(-x)
    }

    vout = dht.apply(vin)

    assert_in_delta 0.18148296716239096,   vout[0], 0.001
    assert_in_delta 0.35680451269699853,   vout[5], 0.001
    assert_in_delta 0.21101009980456306,   vout[10], 0.001
    assert_in_delta 0.02892068100516861,   vout[35], 0.001
    assert_in_delta 0.0022121119664674426, vout[100], 0.001
  end

  def dht4 vin
    dht = GSL::Dht.alloc(N, 1.0, 1.0)

    N.times { |i|
      x = dht.x_sample(i)
      vin[i] = x * (1.0 - x * x)
    }

    vout = dht.apply(vin)

    assert_in_delta  0.05727421417071144,    vout[0], 0.001
    assert_in_delta -0.0001908501261051786,  vout[5], 0.001
    assert_in_delta  2.434180086051901e-05,  vout[10], 0.001
    assert_in_delta -4.0392713194195724e-07, vout[35], 0.001
    assert_in_delta  8.255662619348403e-09,  vout[100], 0.001
  end
end
