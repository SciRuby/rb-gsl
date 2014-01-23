require 'test_helper'

class DhtTest < GSL::TestCase

  N = 128

  def test_dht1
    vin = GSL::Vector.alloc(1, 2, 3)
    dht = GSL::Dht.alloc(3, 1.0, 1.0)

    vout = dht.apply(vin)

    assert_equal  0.3752546494075203,  vout[0]
    assert_equal(-0.13350787269556064, vout[1])
    assert_equal  0.0446799251438404,  vout[2]

    vin2 = dht.apply(vout)
    vin2.scale!(13.323691936314223 ** 2)

    assert_equal 1.0000119186762644, vin2[0]
    assert_equal 1.9999790476647084, vin2[1]
    assert_equal 3.000035803234503,  vin2[2]
  end

  def test_dht2
    vin = GSL::Vector.alloc(N)
    dht = GSL::Dht.alloc(N, 0.0, 100.0)

    N.times { |i|
      x = dht.x_sample(i)
      vin[i] = 1.0 / (1.0 + x * x)
    }

    vout = dht.apply(vin)

    assert_equal 3.999613382195876,   vout[0]
    assert_equal 1.8387637474026606,  vout[5]
    assert_equal 1.2677885358829588,  vout[10]
    assert_equal 0.3521910403797792,  vout[35]
    assert_equal 0.02373661279695407, vout[100]
  end

  def test_dht3
    vin = GSL::Vector.alloc(N)
    dht = GSL::Dht.alloc(N, 1.0, 20.0)

    N.times { |i|
      x = dht.x_sample(i)
      vin[i] = Math.exp(-x)
    }

    vout = dht.apply(vin)

    assert_equal 0.18148296716239096,   vout[0]
    assert_equal 0.35680451269699853,   vout[5]
    assert_equal 0.21101009980456306,   vout[10]
    assert_equal 0.02892068100516861,   vout[35]
    assert_equal 0.0022121119664674426, vout[100]
  end

  def test_dht4
    vin = GSL::Vector.alloc(N)
    dht = GSL::Dht.alloc(N, 1.0, 1.0)

    N.times { |i|
      x = dht.x_sample(i)
      vin[i] = x * (1.0 - x * x)
    }

    vout = dht.apply(vin)

    assert_equal  0.05727421417071144,    vout[0]
    assert_equal(-0.0001908501261051786,  vout[5])
    assert_equal  2.434180086051901e-05,  vout[10]
    assert_equal(-4.0392713194195724e-07, vout[35])
    assert_equal  8.255662619348403e-09,  vout[100]
  end

end
