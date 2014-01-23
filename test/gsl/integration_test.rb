require 'test_helper'

class IntegrationTest < GSL::TestCase

  def test_integration1
    f = GSL::Function.alloc { |x| Math.exp(x) * Math.cos(x) }

    xmin = 0.0
    xmax = 1.0
    limit = 1000

    # QNG
    assert f.integration_qng(xmin, xmax, 0.0, 1.0e-7)
    assert f.integration_qng(xmin, xmax)
    assert f.integration_qng([xmin, xmax])
    assert f.integration_qng([xmin, xmax], [0.0, 1.0e-7])
    assert f.integration_qng([xmin, xmax], 0.0, 1.0e-7)
    assert f.integration_qng(xmin, xmax, [0.0, 1.0e-7])

    # QAG
    assert f.integration_qag(xmin, xmax, 0.0, 1.0e-7, limit, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], 0.0, 1.0e-7, limit, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15)
    assert f.integration_qag(xmin, xmax, [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15)
    assert f.integration_qag(xmin, xmax, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], GSL::Integration::GAUSS15)

    w = GSL::Integration::Workspace.alloc
    assert f.integration_qag(xmin, xmax, [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15, w)
    assert f.integration_qag(xmin, xmax, 0.0, 1.0e-7, GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], 0.0, 1.0e-7, GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], w)

    # QAGS
    assert f.integration_qags(xmin, xmax)
    assert f.integration_qags([xmin, xmax])
    assert f.integration_qags(xmin, xmax, 0.0, 1e-7)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7])
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7])
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7)
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7, limit)
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7], limit)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7], limit)
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7, limit, w)
    assert f.integration_qags(xmin, xmax, 0.0, 1e-7, limit, w)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7], w)
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7], w)

    assert f.integration_qags([xmin, xmax], limit)
    assert f.integration_qags(xmin, xmax, limit)
    assert f.integration_qags([xmin, xmax], w)
    assert f.integration_qags(xmin, xmax, w)
    assert f.integration_qags(xmin, xmax, limit, w)
    assert f.integration_qags([xmin, xmax], limit, w)

    # QAGP
    assert f.integration_qagp([xmin, xmax])
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7])
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7)
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7, limit)
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7], limit)
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7, limit, w)
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7], w)

    assert f.integration_qagp([xmin, xmax], limit)
    assert f.integration_qagp([xmin, xmax], w)
    assert f.integration_qagp([xmin, xmax], limit, w)
  end

  def test_integration2
    f = GSL::Function.alloc { |x| Math.sin(x) / x }

    xmin = 0.0
    xmax = 2.0 * Math::PI
    limit = 1000

    # QNG
    assert f.integration_qng(xmin, xmax, 0.0, 1.0e-7)
    assert f.integration_qng(xmin, xmax)
    assert f.integration_qng([xmin, xmax])
    assert f.integration_qng([xmin, xmax], [0.0, 1.0e-7])
    assert f.integration_qng([xmin, xmax], 0.0, 1.0e-7)
    assert f.integration_qng(xmin, xmax, [0.0, 1.0e-7])

    # QAG
    assert f.integration_qag(xmin, xmax, 0.0, 1.0e-7, limit, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], 0.0, 1.0e-7, limit, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15)
    assert f.integration_qag(xmin, xmax, [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15)
    assert f.integration_qag(xmin, xmax, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], GSL::Integration::GAUSS15)

    w = GSL::Integration::Workspace.alloc(2000)
    assert f.integration_qag(xmin, xmax, [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15, w)
    assert f.integration_qag(xmin, xmax, 0.0, 1.0e-7, GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], 0.0, 1.0e-7, GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15, w)

    # QAGS
    assert f.integration_qags(xmin, xmax)
    assert f.integration_qags([xmin, xmax])
    assert f.integration_qags(xmin, xmax, 0.0, 1e-7)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7])
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7])
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7)
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7, limit)
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7], limit)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7], limit)
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7, limit, w)
    assert f.integration_qags(xmin, xmax, 0.0, 1e-7, limit, w)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7], w)
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7], w)

    assert f.integration_qags([xmin, xmax], limit)
    assert f.integration_qags(xmin, xmax, limit)
    assert f.integration_qags([xmin, xmax], w)
    assert f.integration_qags(xmin, xmax, w)
    assert f.integration_qags(xmin, xmax, limit, w)
    assert f.integration_qags([xmin, xmax], limit, w)

    # QAGP
    assert f.integration_qagp([xmin, xmax])
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7])
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7)
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7, limit)
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7], limit)
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7, limit, w)
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7], w)

    assert f.integration_qagp([xmin, xmax], limit)
    assert f.integration_qagp([xmin, xmax], w)
    assert f.integration_qagp([xmin, xmax], limit, w)
  end

  def test_integration3
    f = GSL::Function.alloc { |x| Math.exp(-x) / Math.sqrt(x) }

    xmin = 0.0
    xmax = 1.0
    limit = 1000

    # QNG
    # XXX GSL::ERROR::ETOL: Ruby/GSL error code 14, failed to reach tolerance with
    # highest-order rule (file qng.c, line 189), failed to reach the specified tolerance
    #assert f.integration_qng(xmin, xmax, 0.0, 1.0e-7)
    #assert f.integration_qng(xmin, xmax)
    #assert f.integration_qng([xmin, xmax])
    #assert f.integration_qng([xmin, xmax], [0.0, 1.0e-7])
    #assert f.integration_qng([xmin, xmax], 0.0, 1.0e-7)
    #assert f.integration_qng(xmin, xmax, [0.0, 1.0e-7])

    # QAG
    assert f.integration_qag(xmin, xmax, 0.0, 1.0e-7, limit, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], 0.0, 1.0e-7, limit, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15)
    assert f.integration_qag(xmin, xmax, [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15)
    assert f.integration_qag(xmin, xmax, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], GSL::Integration::GAUSS15)

    w = GSL::Integration::Workspace.alloc(2000)
    assert f.integration_qag(xmin, xmax, [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15, w)
    assert f.integration_qag(xmin, xmax, 0.0, 1.0e-7, GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], 0.0, 1.0e-7, GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15, w)

    # QAGS
    assert f.integration_qags(xmin, xmax)
    assert f.integration_qags([xmin, xmax])
    assert f.integration_qags(xmin, xmax, 0.0, 1e-7)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7])
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7])
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7)
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7, limit)
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7], limit)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7], limit)
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7, limit, w)
    assert f.integration_qags(xmin, xmax, 0.0, 1e-7, limit, w)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7], w)
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7], w)

    assert f.integration_qags([xmin, xmax], limit)
    assert f.integration_qags(xmin, xmax, limit)
    assert f.integration_qags([xmin, xmax], w)
    assert f.integration_qags(xmin, xmax, w)
    assert f.integration_qags(xmin, xmax, limit, w)
    assert f.integration_qags([xmin, xmax], limit, w)

    # QAGP
    assert f.integration_qagp([xmin, xmax])
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7])
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7)
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7, limit)
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7], limit)
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7, limit, w)
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7], w)

    assert f.integration_qagp([xmin, xmax], limit)
    assert f.integration_qagp([xmin, xmax], w)
    assert f.integration_qagp([xmin, xmax], limit, w)
  end

  def test_integration4
    f = GSL::Function.alloc { |x| (1.0 - Math.exp(-x)) / Math.sqrt(x) }

    xmin = 0.0
    xmax = 0.2
    limit = 1000

    # QNG
    # XXX GSL::ERROR::ETOL: Ruby/GSL error code 14, failed to reach tolerance with
    # highest-order rule (file qng.c, line 189), failed to reach the specified tolerance
    #assert f.integration_qng(xmin, xmax, 0.0, 1.0e-7)
    #assert f.integration_qng(xmin, xmax)
    #assert f.integration_qng([xmin, xmax])
    #assert f.integration_qng([xmin, xmax], [0.0, 1.0e-7])
    #assert f.integration_qng([xmin, xmax], 0.0, 1.0e-7)
    #assert f.integration_qng(xmin, xmax, [0.0, 1.0e-7])

    # QAG
    assert f.integration_qag(xmin, xmax, 0.0, 1.0e-7, limit, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], 0.0, 1.0e-7, limit, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15)
    assert f.integration_qag(xmin, xmax, [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15)
    assert f.integration_qag(xmin, xmax, GSL::Integration::GAUSS15)
    assert f.integration_qag([xmin, xmax], GSL::Integration::GAUSS15)

    w = GSL::Integration::Workspace.alloc(2000)
    assert f.integration_qag(xmin, xmax, [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15, w)
    assert f.integration_qag(xmin, xmax, 0.0, 1.0e-7, GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], 0.0, 1.0e-7, GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], GSL::Integration::GAUSS15, w)
    assert f.integration_qag([xmin, xmax], [0.0, 1.0e-7], limit, GSL::Integration::GAUSS15, w)

    # QAGS
    assert f.integration_qags(xmin, xmax)
    assert f.integration_qags([xmin, xmax])
    assert f.integration_qags(xmin, xmax, 0.0, 1e-7)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7])
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7])
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7)
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7, limit)
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7], limit)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7], limit)
    assert f.integration_qags([xmin, xmax], 0.0, 1e-7, limit, w)
    assert f.integration_qags(xmin, xmax, 0.0, 1e-7, limit, w)
    assert f.integration_qags(xmin, xmax, [0.0, 1e-7], w)
    assert f.integration_qags([xmin, xmax], [0.0, 1e-7], w)

    assert f.integration_qags([xmin, xmax], limit)
    assert f.integration_qags(xmin, xmax, limit)
    assert f.integration_qags([xmin, xmax], w)
    assert f.integration_qags(xmin, xmax, w)
    assert f.integration_qags(xmin, xmax, limit, w)
    assert f.integration_qags([xmin, xmax], limit, w)

    # QAGP
    assert f.integration_qagp([xmin, xmax])
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7])
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7)
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7, limit)
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7], limit)
    assert f.integration_qagp([xmin, xmax], 0.0, 1e-7, limit, w)
    assert f.integration_qagp([xmin, xmax], [0.0, 1e-7], w)

    assert f.integration_qagp([xmin, xmax], limit)
    assert f.integration_qagp([xmin, xmax], w)
    assert f.integration_qagp([xmin, xmax], limit, w)
  end

end
