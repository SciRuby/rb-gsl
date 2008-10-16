=begin
= Error Handling

== Error codes
The GSL routines report an error whenever they cannot perform the task 
requested of them. For example, a root-finding function would return a 
non-zero error code if could not converge to the requested accuracy, 
or exceeded a limit on the number of iterations. Situations like this 
are a normal occurrence when using any mathematical library and 
you should check the return status of the functions that you call.

Whenever a routine reports an error the return value specifies the type of error. 
The return value is analogous to the value of the variable errno in the C library. 
The caller can examine the return code and decide what action to take, including 
ignoring the error if it is not considered serious.

The error code numbers in GSL as (({GSL_EDOM})) are defined in Ruby/GSL 
as Ruby constants under the (({GSL})) module. Here are some of them:
  * (({GSL::EDOM})) - Domain error; used by mathematical functions when an 
    argument value does not fall into the domain over which the function is 
    defined (like (({EDOM})) in the C library)
  * (({GSL::ERANGE})) - Range error; used by mathematical functions when the 
    result value is not representable because of overflow or underflow 
    (like (({ERANGE})) in the C library)
  * (({GSL::ENOMEM})) - No memory available. The system cannot allocate more 
    virtual memory because its capacity is full (like (({ENOMEM})) in the 
    C library). This error is reported when a GSL routine encounters problems 
    when trying to allocate memory with malloc.
  * (({GSL::EINVAL})) - Invalid argument. This is used to indicate various 
    kinds of problems with passing the wrong argument to a library function 
    (like (({EINVAL})) in the C library).

== Error handler
In Ruby/GSL, the default GSL error handler is replaced by an other one which calls
(({rb_raise()})). Thus whenever a GSL routine reports a fatal error, 
a Ruby Exception is generated.

--- GSL::set_error_handler(proc)
--- GSL::set_error_handler { |reason, file, line, errno| ... }
    This replaces the Ruby/GSL default error handler by a user-defined handler
    given by a Proc object ((|proc|)) or a block.

((<prev|URL:use.html>))
((<next|URL:math.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))
=end
