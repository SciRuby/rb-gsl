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
#include "ruby_constants.h"

/*
 * Macros
 */

typedef VALUE (*METHOD)(...);

/*
 * Global Variables
 */

/*
 * Forward Declarations
 */

static VALUE nm_init(int argc, VALUE* argv, VALUE nm);
static VALUE nm_init_cast_copy(VALUE copy, VALUE original, VALUE new_stype, VALUE new_dtype);
static VALUE nm_init_yale_from_old_yale(VALUE shape, VALUE dtype, VALUE ia, VALUE ja, VALUE a, VALUE from_dtype, VALUE from_itype, VALUE nm);
static VALUE nm_dtype(VALUE self);
static VALUE nm_stype(VALUE self);
static VALUE nm_rank(VALUE self);
static VALUE nm_shape(VALUE self);
static VALUE nm_capacity(VALUE self);
static VALUE nm_each(VALUE nmatrix);
static VALUE nm_symmetric(VALUE self);
static VALUE nm_hermitian(VALUE self);

static dtype_t	dtype_from_rbstring(VALUE str);
static dtype_t	dtype_from_rbsymbol(VALUE sym);
static itype_t  itype_from_rbsymbol(VALUE sym);
static dtype_t	dtype_guess(VALUE v);
static dtype_t	interpret_dtype(int argc, VALUE* argv, stype_t stype);
static void*		interpret_initial_value(VALUE arg, dtype_t dtype);
static size_t*	interpret_shape(VALUE arg, size_t* rank);
static stype_t	interpret_stype(VALUE arg);
static stype_t	stype_from_rbstring(VALUE str);
static stype_t	stype_from_rbsymbol(VALUE sym);

#ifdef BENCHMARK
static double get_time(void);
#endif

/*
 * Functions
 */

///////////////////
// Ruby Bindings //
///////////////////

extern "C" {

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
	
//	rb_define_alloc_func(cNMatrix, nm_alloc);
	
	/*
	 * FIXME: These need to be bound in a better way.
	rb_define_singleton_method(cNMatrix, "__cblas_gemm__", nm_cblas_gemm, 13);
	rb_define_singleton_method(cNMatrix, "__cblas_gemv__", nm_cblas_gemv, 11);
	rb_define_singleton_method(cNMatrix, "upcast", nm_upcast, 2);
	 */
	
	//////////////////////
	// Instance Methods //
	//////////////////////
	
	rb_define_method(cNMatrix, "initialize", (VALUE(*)(...))nm_init, -1);
	
	rb_define_method(cNMatrix, "initialize", (METHOD)nm_init, -1);
	
	//rb_define_method(cNMatrix, "initialize_copy", (METHOD)nm_init_copy, 1);
	//rb_define_method(cNMatrix, "initialize_cast_copy", (METHOD)nm_init_cast_copy, 2);
	//rb_define_method(cNMatrix, "as_dtype", (METHOD)nm_cast_copy, 1);
	
	rb_define_method(cNMatrix, "dtype", (METHOD)nm_dtype, 0);
	rb_define_method(cNMatrix, "stype", (METHOD)nm_stype, 0);
	rb_define_method(cNMatrix, "cast",  (METHOD)nm_init_cast_copy, 2);

	rb_define_method(cNMatrix, "[]", (METHOD)nm_mref, -1);
	rb_define_method(cNMatrix, "slice", (METHOD)nm_mget, -1);
	rb_define_method(cNMatrix, "[]=", (METHOD)nm_mset, -1);
	rb_define_method(cNMatrix, "is_ref?", (METHOD)nm_is_ref, 0);
	rb_define_method(cNMatrix, "rank", (METHOD)nm_rank, 0);
	rb_define_method(cNMatrix, "shape", (METHOD)nm_shape, 0);
	rb_define_method(cNMatrix, "transpose", (METHOD)nm_transpose_new, 0);
	rb_define_method(cNMatrix, "det_exact", (METHOD)nm_det_exact, 0);
	rb_define_method(cNMatrix, "transpose!", (METHOD)nm_transpose_self, 0);
	rb_define_method(cNMatrix, "complex_conjugate!", (METHOD)nm_complex_conjugate_bang, 0);

	rb_define_method(cNMatrix, "each", (METHOD)nm_each, 0);

	rb_define_method(cNMatrix, "*", (METHOD)nm_ew_multiply, 1);
	rb_define_method(cNMatrix, "/", (METHOD)nm_ew_divide, 1);
	rb_define_method(cNMatrix, "+", (METHOD)nm_ew_add, 1);
	rb_define_method(cNMatrix, "-", (METHOD)nm_ew_subtract, 1);
	rb_define_method(cNMatrix, "%", (METHOD)nm_ew_mod, 1);
	rb_define_method(cNMatrix, "eql?", (METHOD)nm_eqeq, 1);
	rb_define_method(cNMatrix, "dot", (METHOD)nm_multiply, 1);
	
	/*
	 * TODO: Write new elementwise code for boolean operations
	rb_define_method(cNMatrix, "==", (METHOD)nm_ew_eqeq, 1);
	rb_define_method(cNMatrix, "!=", (METHOD)nm_ew_neq, 1);
	rb_define_method(cNMatrix, "<=", (METHOD)nm_ew_leq, 1);
	rb_define_method(cNMatrix, ">=", (METHOD)nm_ew_geq, 1);
	rb_define_method(cNMatrix, "<", (METHOD)nm_ew_lt, 1);
	rb_define_method(cNMatrix, ">", (METHOD)nm_ew_gt, 1);
	 */

	rb_define_method(cNMatrix, "symmetric?", (METHOD)nm_symmetric, 0);
	rb_define_method(cNMatrix, "hermitian?", (METHOD)nm_hermitian, 0);

	rb_define_method(cNMatrix, "capacity", (METHOD)nm_capacity, 0);

	/*
	 * FIXME: I don't think these should actually be exposed to the Ruby class.
	rb_define_method(cNMatrix, "__yale_ija__", (METHOD)nm_yale_ija, 0);
	rb_define_method(cNMatrix, "__yale_a__", (METHOD)nm_yale_a, 0);
	rb_define_method(cNMatrix, "__yale_size__", (METHOD)nm_yale_size, 0);
	rb_define_method(cNMatrix, "__yale_ia__", (METHOD)nm_yale_ia, 0);
	rb_define_method(cNMatrix, "__yale_ja__", (METHOD)nm_yale_ja, 0);
	rb_define_method(cNMatrix, "__yale_d__", (METHOD)nm_yale_d, 0);
	rb_define_method(cNMatrix, "__yale_lu__", (METHOD)nm_yale_lu, 0);
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
 * End of the Ruby binding functions in the `extern "C" {}` block.
 */
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

  if (!SYMBOL_P(argv[0]) && !RUBYVAL_IS_STRING(argv[0])) {
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

  if (RUBYVAL_IS_NUMERIC(argv[1+offset]) || RUBYVAL_IS_ARRAY(argv[1+offset])) {
  	// Initial value provided (could also be initial capacity, if yale).
  	
    if (stype == YALE_STORE) {
      init_cap = FIX2UINT(argv[1+offset]);
      
    } else {
    	// 4: initial value / dtype
      init_val = interpret_initial_value(argv[1+offset], dtype);
      
      if (RUBYVAL_IS_ARRAY(argv[1+offset])) {
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
  		yale_storage_init((YALE_STORAGE*)(nmatrix->storage));
  		
  		// Do we not need to free the initial value when using other stypes?
  		free(init_val);
  		
  		if (!nmatrix->storage) {
  			rb_raise(rb_eNoMemError, "Yale allocation failed.");
  		}
  		
  		break;
  }

  return nm;
}

/*
 * Create a new NMatrix helper for handling internal ia, ja, and a arguments.
 *
 * This constructor is only called by Ruby code, so we can skip most of the
 * checks.
 */
static VALUE nm_init_yale_from_old_yale(VALUE shape, VALUE dtype, VALUE ia, VALUE ja, VALUE a, VALUE from_dtype, VALUE from_itype, VALUE nm) {
  size_t rank     = 2;
  size_t* shape_  = interpret_shape(shape, &rank);
  dtype_t dtype_  = dtype_from_rbsymbol(dtype);
  char *ia_       = RSTRING_PTR(ia),
       *ja_       = RSTRING_PTR(ja),
       *a_        = RSTRING_PTR(a);
  dtype_t from_dtype_ = dtype_from_rbsymbol(from_dtype);
  itype_t from_itype_ = itype_from_rbsymbol(from_itype);
  NMATRIX* nmatrix;

  UnwrapNMatrix( nm, nmatrix );

  nmatrix->stype   = YALE_STORE;
  nmatrix->storage = (STORAGE*)yale_storage_create_from_old_yale(dtype_, shape_, ia_, ja_, a_, from_dtype_, from_itype_);

  return nm;
}



/*
 * Copy constructor for changing dtypes and stypes.
 */
static VALUE nm_init_cast_copy(VALUE copy, VALUE original, VALUE new_stype_symbol, VALUE new_dtype_symbol) {
  NMATRIX *lhs, *rhs;

  dtype_t new_dtype = dtype_from_rbsymbol(new_dtype_symbol);
  stype_t new_stype = stype_from_rbsymbol(new_stype_symbol);

  CheckNMatrixType(original);

  if (copy == original) return copy;

  UnwrapNMatrix( original, rhs );
  UnwrapNMatrix( copy,     lhs );
  //lhs = ALLOC(NMATRIX); // FIXME: If this fn doesn't work, try switching comments between this line and the above.
  lhs->stype = new_stype;

  // Copy the storage
  static STORAGE* (*ttable[NUM_STYPES][NUM_STYPES])(const STORAGE*, dtype_t) = {
    { dense_storage_cast_copy,  dense_storage_from_list,  dense_storage_from_yale },
    { list_storage_from_dense,  list_storage_cast_copy,   list_storage_from_yale  },
    { yale_storage_from_dense,  yale_storage_from_list,   yale_storage_cast_copy  }
  };

  lhs->storage = ttable[new_stype][rhs->stype](rhs, new_dtype);

  STYPE_MARK_TABLE(mark_table);

  return Data_Wrap_Struct(cNMatrix, mark_table[new_stype], nm_delete, copy);
}


/*
 * Allocator.
 */
static VALUE nm_alloc(VALUE klass) {
  NMATRIX* mat = ALLOC(NMATRIX);
  mat->storage = NULL;
  mat->stype   = NUM_STYPES;

  STYPE_MARK_TABLE(mark_table);

  return Data_Wrap_Struct(klass, mark_table[mat->stype], nm_delete, mat);
}


/*
 * Get the data type (dtype) of a matrix, e.g., :byte, :int8, :int16, :int32,
 * :int64, :float32, :float64, :complex64, :complex128, :rational32,
 * :rational64, :rational128, or :object (the last is a Ruby object).
 */
static VALUE nm_dtype(VALUE self) {
  ID dtype = rb_intern(DTYPE_NAMES[NM_DTYPE(self)]);
  return ID2SYM(dtype);
}


/*
 * Get the storage type (stype) of a matrix, e.g., :yale, :dense, or :list.
 */
static VALUE nm_stype(VALUE self) {
  ID stype = rb_intern(STYPE_NAMES[NM_STYPE(self)]);
  return ID2SYM(stype);
}


/*
 * Find the capacity of an NMatrix. The capacity only differs from the size for Yale matrices, which occasionally
 * allocate more space than they need. For list and dense, capacity gives the number of elements in the matrix.
 */
static VALUE nm_capacity(VALUE self) {
  VALUE cap;

  switch(NM_STYPE(self)) {
  case YALE_STORE:
    cap = UINT2NUM(((YALE_STORAGE*)(NM_STORAGE(self)))->capacity);
    break;

  case DENSE_STORE:
    cap = UINT2NUM(storage_count_max_elements( NM_DENSE_STORAGE(self)->rank, NM_DENSE_STORAGE(self)->shape ));
    break;

  case LIST_STORE:
    cap = UINT2NUM(list_storage_count_elements( NM_LIST_STORAGE(self) ));
    break;

  default:
    rb_raise(nm_eStorageTypeError, "unrecognized stype in nm_capacity()");
  }

  return cap;
}


/*
 * Get the rank of an NMatrix (the number of dimensions).
 *
 * In other words, if you set your matrix to be 3x4, the rank is 2. If the matrix was initialized as 3x4x3, the rank
 * is 3.
 *
 * This function may lie slightly for NVectors, which are internally stored as rank 2 (and have an orientation), but
 * act as if they're rank 1.
 */
static VALUE nm_rank(VALUE self) {
  return rubyobj_from_cval(&(NM_STORAGE(self)->rank), INT64).rval;
}


/*
 * Get the shape (dimensions) of a matrix.
 */
static VALUE nm_shape(VALUE self) {
  STORAGE* s   = NM_STORAGE(self);
  size_t index;

  // Copy elements into a VALUE array and then use those to create a Ruby array with rb_ary_new4.
  VALUE* shape = ALLOCA_N(VALUE, s->rank);
  for (index = 0; index < s->rank; ++index)
    shape[index] = rubyobj_from_cval( s->shape + sizeof(size_t)*index, SIZE_T ).rval;
  //SetFuncs[NM_ROBJ][NM_SIZE_T]( s->rank, shape, sizeof(VALUE), s->shape, sizeof(size_t));

  return rb_ary_new4(s->rank, shape);
}


// Helper function for nm_symmetric and nm_hermitian.
static VALUE is_symmetric(VALUE self, bool hermitian) {
  NMATRIX* m;
  UnwrapNMatrix(self, m);

  if (m->storage->shape[0] == m->storage->shape[1] && m->storage->rank == 2) {

    if (NM_STYPE(self) == DENSE_STORE) {
      if (!hermitian) {
        if (dense_storage_is_symmetric((DENSE_STORAGE*)(m->storage), m->storage->shape[0])) return Qtrue;
      } else if (dense_storage_is_hermitian((DENSE_STORAGE*)(m->storage), m->storage->shape[0])) return Qtrue;
    } else {
      // TODO: Implement, at the very least, yale_is_symmetric. Model it after yale/transp.template.c.
      rb_raise(rb_eNotImpError, "symmetric? and hermitian? only implemented for dense currently");
    }

  }

  return Qfalse;
}


/*
 * Is this matrix symmetric?
 */
static VALUE nm_symmetric(VALUE self) {
  return is_symmetric(self, false);
}

/*
 * Is this matrix hermitian?
 *
 * Definition: http://en.wikipedia.org/wiki/Hermitian_matrix
 *
 * For non-complex matrices, this function should return the same result as symmetric?.
 */
static VALUE nm_hermitian(VALUE self) {
  return is_symmetric(self, true);
}


// Borrowed this function from NArray. Handles 'each' iteration on a dense matrix.
//
// Additionally, handles separately matrices containing VALUEs and matrices containing
// other types of data.
static VALUE nm_dense_each(VALUE nmatrix) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)(NM_STORAGE(nmatrix));
  VALUE v;
  size_t i;

  //void (*copy)();

  if (NM_DTYPE(nmatrix) == RUBYOBJ) {

    // matrix of Ruby objects -- yield those objects directly
    for (i = 0; i < storage_count_max_elements(s->rank, s->shape); ++i)
      rb_yield( *((VALUE*)((char*)(s->elements) + i*DTYPE_SIZES[NM_DTYPE(nmatrix)])) );

  } else {
    // We're going to copy the matrix element into a Ruby VALUE and then operate on it. This way user can't accidentally
    // modify it and cause a seg fault.
    //copy = SetFuncs[NM_ROBJ][NM_DTYPE(nmatrix)];

    for (i = 0; i < storage_count_max_elements(s->rank, s->shape); ++i) {
      v = rubyobj_from_cval((char*)(s->elements) + i*DTYPE_SIZES[NM_DTYPE(nmatrix)], NM_DTYPE(nmatrix)).rval;
      // (*copy)(1, &v, 0, (char*)(s->elements) + i*nm_sizeof[NM_DTYPE(nmatrix)], 0);
      rb_yield(v); // yield to the copy we made
    }
  }

  return nmatrix;
}


/*
 * Iterate over the matrix as you would an Enumerable (e.g., Array).
 *
 * Currently only works for dense.
 */
static VALUE nm_each(VALUE nmatrix) {
  volatile VALUE nm = nmatrix; // not sure why we do this, but it gets done in ruby's array.c.

  switch(NM_STYPE(nm)) {
  case DENSE_STORE:
    return nm_dense_each(nm);
  default:
    rb_raise(rb_eNotImpError, "only dense matrix's each method works right now");
  }
}


/*
 * Check to determine whether matrix is a reference to another matrix.
 */
VALUE nm_is_ref(VALUE self) {
  if (NM_STYPE(self) == DENSE_STORE) // refs only allowed for dense matrices.
    return (NM_DENSE_SRC(self) == NM_STORAGE(self)) ? Qfalse : Qtrue;

  return Qfalse;
}

///////////////////////
// Utility Functions //
///////////////////////


/*
 * Converts a string to a data type.
 */
dtype_t dtype_from_rbstring(VALUE str) {
  size_t index;
  
  for (index = 0; index < NUM_DTYPES; ++index) {
  	if (!strncmp(RSTRING_PTR(str), DTYPE_NAMES[index], RSTRING_LEN(str))) {
  		return static_cast<dtype_t>(index);
  	}
  }
  
  rb_raise(rb_eArgError, "Invalid data type specified.");
}

/*
 * Converts a symbol to a data type.
 */
static dtype_t dtype_from_rbsymbol(VALUE sym) {
  size_t index;
  
  for (index = 0; index < NUM_DTYPES; ++index) {
    if (SYM2ID(sym) == rb_intern(DTYPE_NAMES[index])) {
    	return static_cast<dtype_t>(index);
    }
  }
  
  rb_raise(rb_eArgError, "Invalid data type specified.");
}


/*
 * Converts a symbol to an index type.
 */
static itype_t itype_from_rbsymbol(VALUE sym) {
  size_t index;

  for (index = 0; index < NUM_ITYPES; ++index) {
    if (SYM2ID(sym) == rb_intern(ITYPE_NAMES[index])) {
    	return static_cast<itype_t>(index);
    }
  }

  rb_raise(rb_eArgError, "Invalid index type specified.");
}


/*
 * Guess the data type given a value.
 *
 * TODO: Probably needs some work for Bignum.
 */
static dtype_t dtype_guess(VALUE v) {
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
  	return dtype_from_rbsymbol(argv[offset]);
  	
  } else if (RUBYVAL_IS_STRING(argv[offset])) {
  	return dtype_from_rbstring(StringValue(argv[offset]));
  	
  } else if (stype == YALE_STORE) {
  	rb_raise(rb_eArgError, "Yale storage class requires a dtype.");
  	
  } else {
  	return dtype_guess(argv[0]);
  }
}

/*
 * Convert an Ruby value or an array of Ruby values into initial C values.
 */
static void* interpret_initial_value(VALUE arg, dtype_t dtype) {
  unsigned int index;
  void* init_val;
  
  if (RUBYVAL_IS_ARRAY(arg)) {
  	// Array
    
    init_val = ALLOC_N(int8_t, DTYPE_SIZES[dtype] * RARRAY_LEN(arg));
    for (index = 0; index < RARRAY_LEN(arg); ++index) {
    	rubyval_to_cval(RARRAY_PTR(arg)[index], dtype, (char*)init_val + (index * DTYPE_SIZES[dtype]));
    }
    
  } else {
  	// Single value
  	
    init_val = rubyobj_to_cval(arg, dtype);
  }

  return init_val;
}

/*
 * Convert the shape argument, which may be either a Ruby value or an array of
 * Ruby values, into C values.  The second argument is where the rank of the
 * matrix will be stored.  The function itself returns a pointer to the array
 * describing the shape, which must be freed manually.
 */
static size_t* interpret_shape(VALUE arg, size_t* rank) {
  size_t index;
  size_t* shape;

  if (RUBYVAL_IS_ARRAY(arg)) {
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
  	return stype_from_rbsymbol(arg);
  	
  } else if (RUBYVAL_IS_STRING(arg)) {
  	return stype_from_rbstring(StringValue(arg));
  	
  } else {
  	rb_raise(rb_eArgError, "Expected storage type");
  }
}

/*
 * Converts a string to a storage type. Only looks at the first three
 * characters.
 */
static stype_t stype_from_rbstring(VALUE str) {
  size_t index;
  
  for (index = 0; index < NUM_STYPES; ++index) {
    if (!strncmp(RSTRING_PTR(str), STYPE_NAMES[index], 3)) {
    	return static_cast<stype_t>(index);
    }
  }
  
  return DENSE_STORE;
}

/*
 * Converts a symbol to a storage type.
 */
static stype_t stype_from_rbsymbol(VALUE sym) {
  size_t index;
  
  for (index = 0; index < NUM_STYPES; ++index) {
    if (SYM2ID(sym) == rb_intern(STYPE_NAMES[index])) {
    	return static_cast<stype_t>(index);
    }
  }
  
  return DENSE_STORE;
}


/////////////////
// Exposed API //
/////////////////

extern "C" {

/*
 * Create a dense matrix. Used by the NMatrix GSL fork. Unlike nm_create, this one copies all of the
 * arrays and such passed in -- so you don't have to allocate and pass a new shape object for every
 * matrix you want to create, for example. Same goes for elements.
 *
 * Returns a properly-wrapped Ruby object as a VALUE.
 *
 * TODO: Add a column-major option for libraries that use column-major matrices.
 */
VALUE rb_nmatrix_dense_create(dtype_t dtype, size_t* shape, size_t rank, void* elements, size_t length) {
  NMATRIX* nm;
  VALUE klass;
  size_t nm_rank;
  size_t* shape_copy;

  // Do not allow a rank of 1; if rank == 1, this should probably be an NVector instead, but that still has rank 2.
  if (rank == 1) {
    klass					= cNVector;
    nm_rank				= 2;
    shape_copy		= ALLOC_N(size_t, nm_rank);
    shape_copy[0]	= shape[0];
    shape_copy[1]	= 1;
    
  } else {
    klass				= cNMatrix;
    nm_rank			= rank;
    shape_copy	= ALLOC_N(size_t, nm_rank);
    memcpy(shape_copy, shape, sizeof(size_t)*nm_rank);
  }

  // Copy elements
  void* elements_copy = ALLOC_N(char, DTYPE_SIZES[dtype]*length);
  memcpy(elements_copy, elements, DTYPE_SIZES[dtype]*length);

  // allocate and create the matrix and its storage
  nm = nm_create(DENSE_STORE, dense_storage_create(dtype, shape_copy, rank, elements_copy, length));

  // tell Ruby about the matrix and its storage, particularly how to garbage collect it.
  return Data_Wrap_Struct(klass, dense_storage_mark, dense_storage_delete, nm);
}


/*
 * Create a dense vector. Used by the NMatrix GSL fork.
 *
 * Basically just a convenience wrapper for rb_nmatrix_dense_create().
 *
 * Returns a properly-wrapped Ruby NVector object as a VALUE.
 *
 * TODO: Add a transpose option for setting the orientation of the vector?
 */
VALUE rb_nvector_dense_create(dtype_t dtype, void* elements, size_t length) {
  size_t rank = 1, shape = length;
  return rb_nmatrix_dense_create(dtype, &shape, rank, elements, length);
}

}

