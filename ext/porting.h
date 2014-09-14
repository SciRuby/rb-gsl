#include <ruby/version.h>

#if RUBY_API_VERSION_CODE < 20100
  #define rb_obj_reveal(o,k) RBASIC(o)->klass = k
#endif
