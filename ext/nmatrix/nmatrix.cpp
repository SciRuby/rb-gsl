/////////////////////////////////////////////////////////////////////
// = NMatrix
//
// A linear algebra library for scientific computation in Ruby.
// NMatrix is part of SciRuby.
//
// NMatrix was originally inspired by and derived from NArray, by
// Masahiro Tanaka: http://narray.rubyforge.org
//
// == Copyright Information
//
// SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
// NMatrix is Copyright (c) 2012, Ruby Science Foundation
//
// Please see LICENSE.txt for additional copyright notices.
//
// == Contributing
//
// By contributing source code to SciRuby, you agree to be bound by
// our Contributor Agreement:
//
// * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
//
// == nmatrix.c
//

/////////////////////////////////////////////////////////////////////
// = NMatrix
//
// A linear algebra library for scientific computation in Ruby.
// NMatrix is part of SciRuby.
//
// NMatrix was originally inspired by and derived from NArray, by
// Masahiro Tanaka: http://narray.rubyforge.org
//
// == Copyright Information
//
// SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
// NMatrix is Copyright (c) 2012, Ruby Science Foundation
//
// Please see LICENSE.txt for additional copyright notices.
//
// == Contributing
//
// By contributing source code to SciRuby, you agree to be bound by
// our Contributor Agreement:
//
// * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
//
// == data.cpp
//
// Functions and data for dealing the data types.

/*
 * Standard Includes
 */

#include <ruby.h>

/*
 * Project Includes
 */

#include "nmatrix.h"
#include "ruby_symbols.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

/*
 * Forward Declarations
 */
 
static dtype_t	dtype_from_string(VALUE str);
static dtype_t	dtype_from_symbol(VALUE sym);
static dtype_t	dtype_guess(VALUE v);
static double		get_time(void);
static dtype_t	interpret_dtype(int argc, VALUE* argv, stype_t stype);
static stype_t	interpret_stype(VALUE arg);
static stype_t	stype_from_string(VALUE str);
static stype_t	stype_from_symbol(VALUE sym);

/*
 * Functions
 */

///////////////////
// Ruby Bindings //
///////////////////

void Init_nmatrix() {
	
	///////////////////////
	// Class Definitions //
	///////////////////////
	
	cNMatrix = rb_define_class("NMatrix", rb_cObject);
	cNVector = rb_define_class("NVector", cNMatrix);
	
	// Special exceptions
	nm_eDataTypeError    = rb_define_class("DataTypeError",			rb_eStandardError);
	nm_eStorageTypeError = rb_define_class("StorageTypeError",	rb_eStandardError);

	///////////////////
	// Class Methods //
	///////////////////
	
	rb_define_alloc_func(cNMatrix, nm_alloc);
	
	/*
	 * FIXME: These need to be bound in a better way.
	rb_define_singleton_method(cNMatrix, "__cblas_gemm__", nm_cblas_gemm, 13);
	rb_define_singleton_method(cNMatrix, "__cblas_gemv__", nm_cblas_gemv, 11);
	rb_define_singleton_method(cNMatrix, "upcast", nm_upcast, 2);
	 */
	
	//////////////////////
	// Instance Methods //
	//////////////////////
	
	rb_define_method(cNMatrix, "initialize", nm_init, -1);
	
	rb_define_method(cNMatrix, "initialize_copy", nm_init_copy, 1);
	rb_define_method(cNMatrix, "initialize_cast_copy", nm_init_cast_copy, 2);
	rb_define_method(cNMatrix, "as_dtype", nm_cast_copy, 1);
	
	rb_define_method(cNMatrix, "dtype", nm_dtype, 0);
	rb_define_method(cNMatrix, "stype", nm_stype, 0);
	rb_define_method(cNMatrix, "cast",  nm_scast_copy, 2);

	rb_define_method(cNMatrix, "[]", nm_mref, -1);
	rb_define_method(cNMatrix, "[]=", nm_mset, -1);
	rb_define_method(cNMatrix, "rank", nm_rank, 0);
	rb_define_method(cNMatrix, "shape", nm_shape, 0);
	rb_define_method(cNMatrix, "transpose", nm_transpose_new, 0);
	rb_define_method(cNMatrix, "det_exact", nm_det_exact, 0);
	rb_define_method(cNMatrix, "transpose!", nm_transpose_self, 0);
	rb_define_method(cNMatrix, "complex_conjugate!", nm_complex_conjugate_bang, 0);

	rb_define_method(cNMatrix, "each", nm_each, 0);

	rb_define_method(cNMatrix, "*", nm_ew_multiply, 1);
	rb_define_method(cNMatrix, "/", nm_ew_divide, 1);
	rb_define_method(cNMatrix, "+", nm_ew_add, 1);
	rb_define_method(cNMatrix, "-", nm_ew_subtract, 1);
	rb_define_method(cNMatrix, "%", nm_ew_mod, 1);
	rb_define_method(cNMatrix, "eql?", nm_eqeq, 1);
	rb_define_method(cNMatrix, "dot", nm_multiply, 1);
	
	/*
	 * TODO: Write new elementwise code for boolean operations
	rb_define_method(cNMatrix, "==", nm_ew_eqeq, 1);
	rb_define_method(cNMatrix, "!=", nm_ew_neq, 1);
	rb_define_method(cNMatrix, "<=", nm_ew_leq, 1);
	rb_define_method(cNMatrix, ">=", nm_ew_geq, 1);
	rb_define_method(cNMatrix, "<", nm_ew_lt, 1);
	rb_define_method(cNMatrix, ">", nm_ew_gt, 1);
	 */

	rb_define_method(cNMatrix, "symmetric?", nm_symmetric, 0);
	rb_define_method(cNMatrix, "hermitian?", nm_hermitian, 0);

	rb_define_method(cNMatrix, "capacity", nm_capacity, 0);

	/*
	 * FIXME: I don't think these should actually be exposed to the Ruby class.
	rb_define_method(cNMatrix, "__yale_ija__", nm_yale_ija, 0);
	rb_define_method(cNMatrix, "__yale_a__", nm_yale_a, 0);
	rb_define_method(cNMatrix, "__yale_size__", nm_yale_size, 0);
	rb_define_method(cNMatrix, "__yale_ia__", nm_yale_ia, 0);
	rb_define_method(cNMatrix, "__yale_ja__", nm_yale_ja, 0);
	rb_define_method(cNMatrix, "__yale_d__", nm_yale_d, 0);
	rb_define_method(cNMatrix, "__yale_lu__", nm_yale_lu, 0);
	rb_define_const(cNMatrix, "YALE_GROWTH_CONSTANT", rb_float_new(YALE_GROWTH_CONSTANT));
	*/
	
	/////////////
	// Aliases //
	/////////////
	
	rb_define_alias(cNMatrix, "dim", "rank");
	rb_define_alias(cNMatrix, "equal?", "eql?");
	
	///////////////////////
	// Symbol Generation //
	///////////////////////
	
	ruby_symbols_init();
}

/*
 * Create a new NMatrix.
 *
 * There are several ways to do this. At a minimum, dimensions and either a dtype or initial values are needed, e.g.,
 *
 *     NMatrix.new(3, :int64)       # square 3x3 dense matrix
 *     NMatrix.new([3,4], :float32) # 3x4 matrix
 *     NMatrix.new(3, 0)            # 3x3 dense matrix initialized to all zeros
 *     NMatrix.new([3,3], [1,2,3])  # [[1,2,3],[1,2,3],[1,2,3]]
 *
 * NMatrix will try to guess the dtype from the first value in the initial values array.
 *
 * You can also provide the stype prior to the dimensions. However, non-dense matrices cannot take initial values, and
 * require a dtype (e.g., :int64):
 *
 *     NMatrix.new(:yale, [4,3], :int64)
 *     NMatrix.new(:list, 5, :rational128)
 *
 * For Yale, you can also give an initial size for the non-diagonal component of the matrix:
 *
 *     NMatrix.new(:yale, [4,3], 2, :int64)
 *
 * Finally, you can be extremely specific, and define a matrix very exactly:
 *
 *     NMatrix.new(:dense, [2,2,2], [0,1,2,3,4,5,6,7], :int8)
 *
 * There is one additional constructor for advanced users, which takes seven arguments and is only for creating Yale matrices
 * with known IA, JA, and A arrays. This is used primarily internally for IO, e.g., reading Matlab matrices, which are
 * stored in old Yale format.
 *
 * Just be careful! There are no overflow warnings in NMatrix.
 */
static VALUE nm_init(int argc, VALUE* argv, VALUE nm) {
  char    ZERO = 0;
  VALUE   QNIL = Qnil;
  dtype_t dtype;
  stype_t stype;
  size_t  rank, offset = 0;
  size_t* shape;
  size_t  init_cap = 0, init_val_len = 0;
  void*   init_val = NULL;
  NMATRIX* nmatrix;

  /*
   * READ ARGUMENTS
   */

  if (argc < 2) {
  	rb_raise(rb_eArgError, "Expected 2-4 arguments (or 8 for internal Yale creation)");
  	return Qnil;
  }

  if (!SYMBOL_P(argv[0]) && !IS_STRING(argv[0])) {
    stype = DENSE_STORE;
    
  } else {
  	// 0: String or Symbol
    stype  = interpret_stype(argv[0]);
    offset = 1;
  }

  // If there are 7 arguments and Yale, refer to a different init function with fewer sanity checks.
  if (argc == 8) {
    if (stype == YALE_STORE) {
    	return nm_init_yale_from_old_yale(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], nm);
    	
		} else {
		  rb_raise(rb_eArgError, "Expected 2-4 arguments (or 7 for internal Yale creation)");
		  return Qnil;
		}
  }
	
	// 1: Array or Fixnum
  shape = interpret_shape(argv[offset], &rank);
  // 2-3: dtype
  dtype = interpret_dtype(argc-1-offset, argv+offset+1, stype);

  if (IS_NUMERIC(argv[1+offset]) || TYPE(argv[1+offset]) == T_ARRAY) {
  	// Initial value provided (could also be initial capacity, if yale).
  	
    if (stype == YALE_STORE) {
      init_cap = FIX2UINT(argv[1+offset]);
      
    } else {
    	// 4: initial value / dtype
      init_val = interpret_initial_value(argv[1+offset], dtype);
      
      if (TYPE(argv[1+offset]) == T_ARRAY) {
      	init_val_len = RARRAY_LEN(argv[1+offset]);
      	
      } else {
      	init_val_len = 1;
      }
    }
    
  } else {
  	// DType is RUBYOBJ.
  	
    if (stype == DENSE_STORE) {
    	/*
    	 * No need to initialize dense with any kind of default value unless it's
    	 * an RUBYOBJ matrix.
    	 */
      if (dtype == RUBYOBJ) {
      	// Pretend [nil] was passed for RUBYOBJ.
      	init_val = ALLOC(VALUE);
        *(VALUE*)init_val = QNIL;
        
        init_val_len = 1;
        
      } else {
      	init_val = NULL;
      }
      
    } else if (stype == LIST_STORE) {
    	init_val = ALLOC_N(int8_t, DTYPE_SIZES[dtype]);
      memset(init_val, 0, DTYPE_SIZES[dtype]);
    }
  }
	
  // TODO: Update to allow an array as the initial value.
	
  UnwrapNMatrix(nm, nmatrix);

  nmatrix->stype = stype;
  
  switch (stype) {
  	case DENSE_STORE:
  		nmatrix->storage = (STORAGE*)dense_storage_create(dtype, shape, rank, init_val, init_val_len);
  		break;
  		
  	case LIST_STORE:
  		nmatrix->storage = (STORAGE*)list_storage_create(dtype, shape, rank, init_val);
  		break;
  		
  	case YALE_STORE:
  		nmatrix->storage = (STORAGE*)yale_storage_create(dtype, shape, rank, init_cap);
  		yale_storage_init(nmatrix->storage);
  		
  		// Do we not need to free the initial value when using other stypes?
  		free(init_val);
  		
  		if (!nmatrix->storage) {
  			rb_raise(rb_eNoMemError, "Yale allocation failed.");
  		}
  		
  		break;
  }

  return nm;
}

///////////////////////
// Utility Functions //
///////////////////////

/*
 * Converts a string to a data type.
 */
dtype_t dtype_from_string(VALUE str) {
  size_t index;
  
  for (index = 0; index < NUM_DTYPES; ++index) {
  	if (!strncmp(RSTRING_PTR(str), DTYPE_NAMES[index], RSTRING_LEN(str))) {
  		return index;
  	}
  }
  
  rb_raise(rb_eArgError, "Invalid data type speicified.");
}

/*
 * Converts a symbol to a data type.
 */
dtype_t dtype_from_symbol(VALUE sym) {
  size_t index;
  
  for (index = 0; index < NUM_DTYPES; ++index) {
    if (SYM2ID(sym) == rb_intern(nm_dtypestring[i])) {
    	return index;
    }
  }
  
  rb_raise(rb_eArgError, "Invalid data type speicified.");
}

/*
 * Guess the data type given a value.
 *
 * TODO: Probably needs some work for Bignum.
 */
dtype_t dtype_guess(VALUE v) {
  switch(TYPE(v)) {
  case T_TRUE:
  case T_FALSE:
    return BYTE;
    
  case T_STRING:
    if (RSTRING_LEN(v) == 1) {
    	return BYTE;
    	
    } else {
    	rb_raise(rb_eArgError, "Strings of length > 1 may not be stored in a matrix.");
    }

#if SIZEOF_INT == 8
  case T_FIXNUM:
    return INT64;
    
  case T_RATIONAL:
    return RATIONAL128;
    
#else
# if SIZEOF_INT == 4
  case T_FIXNUM:
    return INT32;
    
  case T_RATIONAL:
    return RATIONAL64;
    
#else
  case T_FIXNUM:
    return INT16;
    
  case T_RATIONAL:
    return RATIONAL32;
# endif
#endif

  case T_BIGNUM:
    return INT64;

#if SIZEOF_FLOAT == 4
  case T_COMPLEX:
    return COMPLEX128;
    
  case T_FLOAT:
    return FLOAT64;
    
#else
# if SIZEOF_FLOAT == 2
  case T_COMPLEX:
    return COMPLEX64;
    
  case T_FLOAT:
    return FLOAT32;
# endif
#endif

  case T_ARRAY:
  	/*
  	 * May be passed for dense -- for now, just look at the first element.
  	 * 
  	 * TODO: Look at entire array for most specific type.
  	 */
  	
    return dtype_guess(RARRAY_PTR(v)[0]);

  case T_NIL:
  default:
    rb_raise(rb_eArgError, "Unable to guess a data type from provided parameters; data type must be specified manually.");
  }
}

#ifdef BENCHMARK
/*
 * A simple function used when benchmarking NMatrix.
 */
static double get_time(void) {
  struct timeval t;
  struct timezone tzp;
  
  gettimeofday(&t, &tzp);
  
  return t.tv_sec + t.tv_usec*1e-6;
}
#endif

/*
 * The argv parameter will be either 1 or 2 elements.  If 1, could be either
 * initial or dtype.  If 2, is initial and dtype. This function returns the
 * dtype.
 */
static dtype_t interpret_dtype(int argc, VALUE* argv, stype_t stype) {
  int offset;
  
  switch (argc) {
  	case 1:
  		offset = 0;
  		break;
  	
  	case 2:
  		offset = 1;
  		break;
  		
  	default:
  		rb_raise(rb_eArgError, "Need an initial value or a dtype.");
  		break;
  }

  if (SYMBOL_P(argv[offset])) {
  	return dtype_from_symbol(argv[offset]);
  	
  } else if (IS_STRING(argv[offset])) {
  	return dtype_from_string(StringValue(argv[offset]));
  	
  } else if (stype == S_YALE) {
  	rb_raise(rb_eArgError, "Yale storage class requires a dtype.");
  	
  } else {
  	return dtype_guess(argv[0]);
  }
}

/*
 * Convert an Ruby value or an array of Ruby values into initial C values.
 *
 * FIXME: Remove the references to SetFuncs.
 */
void* interpret_initial_value(VALUE arg, dtype_t dtype) {
  void* init_val;

  if (TYPE(arg) == T_ARRAY) {
  	// Array
    
    init_val = ALLOC_N(int8_t, DTYPE_SIZES[dtype] * RARRAY_LEN(arg));
    SetFuncs[dtype][RUBYOBJ](RARRAY_LEN(arg), init_val, DTYPE_SIZES[dtype], RARRAY_PTR(arg), DTYPE_SIZES[RUBYOBJ]);
    
  } else {
  	// Single value
  	
    init_val = ALLOC_N(int8_t, DTYPE_SIZES[dtype]);
    SetFuncs[dtype][RUBYOBJ](1, init_val, 0, &arg, 0);
  }

  return init_val;
}

/*
 * Convert the shape argument, which may be either a Ruby value or an array of
 * Ruby values, into C values.  The second argument is where the rank of the
 * matrix will be stored.  The function itself returns a pointer to the array
 * describing the shape, which must be freed manually.
 */
size_t* interpret_shape(VALUE arg, size_t* rank) {
  size_t index;
  size_t* shape;

  if (TYPE(arg) == T_ARRAY) {
    *rank = RARRAY_LEN(arg);
    shape = ALLOC_N(size_t, *rank);
    
    for (index = *rank; index-- > 0;) {
      shape[index] = (size_t)(FIX2UINT(RARRAY_PTR(arg)[index]));
    }
    
  } else if (FIXNUM_P(arg)) {
    *rank = 2;
    shape = ALLOC_N(size_t, *rank);
    
    shape[0] = (size_t)FIX2UINT(arg);
    shape[1] = (size_t)FIX2UINT(arg);
    
  } else {
    rb_raise(rb_eArgError, "Expected an array of numbers or a single fixnum for matrix shape");
  }
	
  return shape;
}

/*
 * Convert a Ruby symbol or string into an storage type.
 */
static stype_t interpret_stype(VALUE arg) {
  if (SYMBOL_P(arg)) {
  	return stype_from_symbol(arg);
  	
  } else if (IS_STRING(arg)) {
  	return stype_from_string(StringValue(arg));
  	
  } else {
  	rb_raise(rb_eArgError, "Expected storage type");
  }
}

/*
 * Converts a string to a storage type. Only looks at the first three
 * characters.
 */
static stype_t stype_from_string(VALUE str) {
  size_t index;
  
  for (index = 0; index < NUM_STYPES; ++index) {
    if (!strncmp(RSTRING_PTR(str), STYPE_NAMES[index], 3)) {
    	return index;
    }
  }
  
  return DENSE_STORE;
}

/*
 * Converts a symbol to a storage type.
 */
static stype_t stype_from_symbol(VALUE sym) {
  
  size_t index;
  for (index = 0; index < NUM_STYPES; ++index) {
    if (SYM2ID(sym) == rb_intern(STYPE_NAMES[index])) {
    	return index;
    }
  }
  
  return DENSE_STORE;
}

