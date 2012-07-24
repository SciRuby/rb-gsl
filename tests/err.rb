#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "sum/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test

MAX_ERRS = 64

class ErrCode
  attr_accessor :number
  attr_accessor :name
  def initialize(number, name)
    @number = number
    @name = name
  end
end

errors = ["GSL::SUCCESS", "GSL::FAILURE", "GSL::CONTINUE", "GSL::EDOM", 
          "GSL::ERANGE", "GSL::EFAULT", "GSL::EINVAL", "GSL::EFAILED", 
          "GSL::EFACTOR", "GSL::ESANITY", "GSL::ENOMEM", "GSL::EBADFUNC", 
          "GSL::ERUNAWAY", "GSL::EMAXITER", "GSL::EZERODIV", "GSL::EBADTOL", 
          "GSL::ETOL", "GSL::EUNDRFLW", "GSL::EOVRFLW", "GSL::ELOSS", 
          "GSL::EROUND", "GSL::EBADLEN", "GSL::ENOTSQR", "GSL::ESING",
          "GSL::EDIVERGE", "GSL::EUNSUP", "GSL::EUNIMPL", "GSL::ECACHE",
          "GSL::ETABLE", "GSL::ENOPROG", "GSL::ENOPROGJ", "GSL::ETOLF", 
          "GSL::ETOLX", "GSL::ETOLG", "GSL::EOF"]

Errors = Array.new(errors.size)

errors.each_index do |i|
  Errors[i] = ErrCode.new(eval("#{errors[i]}"), errors[i])
end

for i in 0...Errors.size
  printf("%s = %d\n", Errors[i].name, Errors[i].number)
end

for i in 0...Errors.size
  status = 0
  for j in 0...Errors.size
    if j != i
      status |= (Errors[i].number == Errors[j].number) ? 1 : 0
    end
    GSL::Test::test(status, "#{Errors[i].name} is distinct from other error values")
  end
end

for i in 0...Errors.size
  status = 0
  e1 = Errors[i].number
  for j in 0...Errors.size
    if j != i
      e2 = Errors[j].number
      status |= (GSL::strerror(e1) == GSL::strerror(e2)) ? 1 : 0
    end
  end
  GSL::Test::test(status, "#{Errors[i].name} has a distinct error message")
end
