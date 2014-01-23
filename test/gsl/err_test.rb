require 'test_helper'

class ErrTest < GSL::TestCase

  MAX_ERRS = 64

  ERRORS = %w[
    SUCCESS FAILURE CONTINUE EDOM ERANGE EFAULT EINVAL EFAILED
    EFACTOR ESANITY ENOMEM EBADFUNC ERUNAWAY EMAXITER EZERODIV
    EBADTOL ETOL EUNDRFLW EOVRFLW ELOSS EROUND EBADLEN ENOTSQR
    ESING EDIVERGE EUNSUP EUNIMPL ECACHE ETABLE ENOPROG ENOPROGJ
    ETOLF ETOLX ETOLG EOF
  ].map { |name| GSL.const_get(name) }

  def test_number
    assert ERRORS.uniq == ERRORS
  end

  def test_message
    assert ERRORS.map { |e| GSL.strerror(e) }.uniq.size == ERRORS.size
  end

end
