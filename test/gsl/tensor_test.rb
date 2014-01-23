require 'test_helper'

class TensorTest < GSL::TestCase

  RANK = 3
  DIMENSION = 5

  def test_tensor
    return unless GSL.const_defined?(:Tensor)

    t = GSL::Tensor.alloc(RANK, DIMENSION)

    assert t.rank == RANK, "#{t.class}.alloc returns valid rank"
    assert t.dimension == DIMENSION, "#{t.class}_alloc returns valid dimension"

    counter = 0
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          counter += 1
          t.set(i, j, k, counter)
    #      t[i, j, k] = counter
        end
      end
    end

    status = 0
    counter = 0
    data = t.data    # GSL::Vector::View
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          counter += 1
          if data[DIMENSION * DIMENSION * i + DIMENSION * j + k] != counter
            status += 1
          end
        end
      end
    end

    assert status.zero?, "#{t.class}#set writes into array"

    status = 0
    counter = 0
    data = t.data
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          counter += 1
    #      if t.get(i, j, k) != counter
          if t[i, j, k] != counter
            status += 1
          end
        end
      end
    end

    assert status.zero?, "#{t.class}#get reads from array"

    t = GSL::Tensor.calloc(RANK, DIMENSION)

    counter = 0
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          counter += 1
          t.set(i, j, k, counter)
        end
      end
    end

    exp_max = t[0, 0, 0]
    exp_min = t[0, 0, 0]
    exp_imax = exp_jmax = exp_kmax = 0
    exp_imin = exp_jmin = exp_kmin = 0

    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          value = t[i, j, k]
          if value > exp_max
            exp_max = value
            exp_imax = i; exp_jmax = j; exp_kmax = k
          end
          if value < exp_min
            exp_min = t[i, j, k]
            exp_imin = i; exp_jmin = j; exp_kmin = k
          end
        end
      end
    end

    max = t.max
    assert max == exp_max, "#{t.class}#max returns correct maximum value"
    min = t.min
    assert min == exp_min, "#{t.class}#min returns correct minimum value"

    min, max = t.minmax
    assert max == exp_max, "#{t.class}#minmax returns correct maximum value"
    assert min == exp_min, "#{t.class}#minmax returns correct minimum value"

    imax = t.max_index
    status = 0
    if imax[0] != exp_imax; status += 1; end
    if imax[1] != exp_jmax; status += 1; end
    if imax[2] != exp_kmax; status += 1; end
    assert status.zero?, "#{t.class}#max_index returns correct maximum indices"

    imin = t.min_index
    status = 0
    if imin[0] != exp_imin; status += 1; end
    if imin[1] != exp_jmin; status += 1; end
    if imin[2] != exp_kmin; status += 1; end
    assert status.zero?, "#{t.class}#min_index returns correct minimum indices"


    imin, imax = t.minmax_index
    status = 0
    if imin[0] != exp_imin; status += 1; end
    if imin[1] != exp_jmin; status += 1; end
    if imin[2] != exp_kmin; status += 1; end
    if imax[0] != exp_imax; status += 1; end
    if imax[1] != exp_jmax; status += 1; end
    if imax[2] != exp_kmax; status += 1; end
    assert status.zero?, "#{t.class}#minmax_index returns correct indices"

    ##### Operations
    a = GSL::Tensor.new(RANK, DIMENSION)
    b = GSL::Tensor.new(RANK, DIMENSION)
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          a[i, j, k] = 3 + i + 5 * j + 2 * k
          b[i, j, k] = 3 + 2 * i + 4 * j + k
        end
      end
    end

    # Addition
    c = a + b
    #c = a.add(b)
    status = 0
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          r = c[i, j, k]
          x = a[i, j, k]
          y = b[i, j, k]
          z = x + y
          status += 1 if r != z
        end
      end
    end
    assert status.zero?, "#{t.class}#add tensor addition"

    # Subtraction
    c = a - b
    # c = a.sub(b)
    status = 0
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          r = c[i, j, k]
          x = a[i, j, k]
          y = b[i, j, k]
          z = x - y
          status += 1 if r != z
        end
      end
    end
    assert status.zero?, "#{t.class}#sub tensor subtraction"

    # Element multiplication
     c = a.mul_elements(b)
    status = 0
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          r = c[i, j, k]
          x = a[i, j, k]
          y = b[i, j, k]
          z = x * y
          status += 1 if r != z
        end
      end
    end
    assert status.zero?, "#{t.class}#mul_elements element multiplication"

    # Element division
    c = a.div_elements(b)
    status = 0
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          r = c[i, j, k]
          x = a[i, j, k]
          y = b[i, j, k]
          z = x / y
          if (r - z).abs > 2 * GSL::FLT_EPSILON * z.abs; status += 1; end
        end
      end
    end
    assert status.zero?, "#{t.class}#div_elements element division"

    ### Tensor product
    c = a * b
    status = 0
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          for l in 0...DIMENSION do
            for m in 0...DIMENSION do
              for n in 0...DIMENSION do
                r = c[i, j, k, l, m, n]
                x = a[i, j, k]
                y = b[l, m, n]
                z = x * y
                if r != z; status += 1; end
              end
            end
          end
        end
      end
    end
    assert status.zero?, "#{t.class}#product tensorial product"

    ### Index contraction
    tt = a.contract(0, 1)
    assert tt.rank == RANK - 2, "#{t.class}.contract returns valid rank"
    assert tt.dimension == DIMENSION, "#{t.class}_contract returns valid dimension"

    ### Swap indices
    a_102 = a.swap_indices(0, 1)
    a_210 = a.swap_indices(0, 2)
    a_021 = a.swap_indices(1, 2)
    status = 0
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          x = a[i, j, k]
          x_102 = a_102[j, i, k]
          x_210 = a_210[k, j, i]
          x_021 = a_021[i, k, j]
          if x != x_102 or x != x_210 or x != x_021; status += 1; end
        end
      end
    end
    assert status.zero?, "#{t.class}#swap_indices swap indices"

    ### Test text IO
    file = 'tensor_test.txt'

    t = GSL::Tensor.alloc(RANK, DIMENSION)
    counter = 0
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          counter += 1
          t[i, j, k] = counter
        end
      end
    end

    t.fprintf(file, '%g')
    tt = GSL::Tensor.alloc(RANK, DIMENSION)
    status = 0
    tt.fscanf(file)
    counter = 0
    data = tt.data
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          counter += 1
          if data[DIMENSION * DIMENSION * i + DIMENSION * j + k] != counter
            status += 1
          end
        end
      end
    end
    assert status.zero?, "#{t.class}#fprintf and fscanf"
    File.delete(file)

    ### Test binary IO
    file = 'tensor_test.dat'
    t.fwrite(file)
    tt = GSL::Tensor.alloc(RANK, DIMENSION)
    status = 0
    tt.fread(file)
    counter = 0
    data = tt.data
    for i in 0...DIMENSION do
      for j in 0...DIMENSION do
        for k in 0...DIMENSION do
          counter += 1
          if data[DIMENSION * DIMENSION * i + DIMENSION * j + k] != counter
            status += 1
          end
        end
      end
    end
    assert status.zero?, "#{t.class}#fwrite and fread"
    File.delete(file)

    ### Trap
    i = j = k = 0
    t = GSL::Tensor.calloc(RANK, DIMENSION)
    assert_nothing_raised("#{t.class}#set traps 1st index above upper bound") {
      t[DIMENSION + 1, 0, 0] = 1.2
    }

    assert_nothing_raised("#{t.class}#set traps 2nd index above upper bound") {
      t[0, DIMENSION + 1, 0] = 1.2
    }

    assert_nothing_raised("#{t.class}#set traps 3rd index above upper bound") {
      t[0, 0, DIMENSION + 1] = 1.2
    }

    assert_nothing_raised("#{t.class}#set traps 2nd index at upper bound") {
      t[0, DIMENSION, 0] = 1.2
    }

    assert_nothing_raised("#{t.class}#set traps 2nd index below lower bound") {
      t[0, i - 1, 0] = 1.2
    }

    assert_nothing_raised("#{t.class}#get traps 1st index above upper bound") {
      t[DIMENSION + 1, 0, 0]
    }

    assert_nothing_raised("#{t.class}#get traps 2nd index above upper bound") {
      t[0, DIMENSION + 1, 0]
    }

    assert_nothing_raised("#{t.class}#get traps 3rd index above upper bound") {
      t[0, 0, DIMENSION + 1]
    }

    assert_nothing_raised("#{t.class}#get traps 2nd index at upper bound") {
      t[0, DIMENSION, 0]
    }

    assert_nothing_raised("#{t.class}#get traps 2nd index below lower bound") {
      t[0, i - 1, 0]
    }

    #####
    # Vector and Tensor, subtensors
    #v = GSL::Vector.new(0...125)
    v = GSL::Vector.indgen(125)
    t = v.to_tensor(3, 5)
    assert t.rank == RANK, "#{v.class}.to_tensor(#{RANK}, #{DIMENSION}) returns valid rank"
    assert t.dimension == DIMENSION, "#{v.class}.to_tensor(#{RANK}, #{DIMENSION}) returns valid dimension"

    m0_exp = GSL::Matrix[0...25, 5, 5]
    m1_exp = GSL::Matrix[25...50, 5, 5]
    m2_exp = GSL::Matrix[50...75, 5, 5]
    m3_exp = GSL::Matrix[75...100, 5, 5]
    m4_exp = GSL::Matrix[100...125, 5, 5]

    # Create tensors of rank 2
    t0 = t.subtensor(0)
    t1 = t[1]
    t2 = t.subtensor(2)
    t3 = t[3]
    t4 = t.subtensor(4)

    # 2-tensors can be compared directly with matrices
    assert t0 == m0_exp, "#{t.class}#subtensor(0) returns valid tensor"
    assert t1 == m1_exp, "#{t.class}#subtensor(1) returns valid tensor"
    assert t2 == m2_exp, "#{t.class}#subtensor(2) returns valid tensor"
    assert t3 == m3_exp, "#{t.class}#subtensor(3) returns valid tensor"
    assert t4 == m4_exp, "#{t.class}#subtensor(4) returns valid tensor"

    v0_exp = GSL::Vector[100...105]
    v1_exp = GSL::Vector[105...110]
    v2_exp = GSL::Vector[110...115]
    v3_exp = GSL::Vector[115...120]
    v4_exp = GSL::Vector[120...125]

    # Create tensors of rank1
    v0 = t[4, 0]
    v1 = t[4][1]
    v2 = t.subtensor(4, 2)
    v3 = t4[3]
    v4 = t4.subtensor(4)

    # 1-tensors can be compared directly with vectors
    assert v0 == v0_exp, "#{t.class}#subtensor(4,0) returns valid tensor"
    assert v1 == v1_exp, "#{t.class}#subtensor(4,1) returns valid tensor"
    assert v2 == v2_exp, "#{t.class}#subtensor(4,2) returns valid tensor"
    assert v3 == v3_exp, "#{t.class}#subtensor(4,3) returns valid tensor"
    assert v4 == v4_exp, "#{t.class}#subtensor(4,4) returns valid tensor"
  end

end
