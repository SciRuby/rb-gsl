/*
  gsl_nmatrix.c

  Written by John Woods     Jul/2012
 */

#include "rb_gsl_config.h"
//#ifdef HAVE_NMATRIX_H
#include "rb_gsl_array.h"
#include "rb_gsl_with_nmatrix.h"

// comment out masa 20130221
/*
typedef enum {
	NM_BYTE				  =  0, // unsigned char
	NM_INT8				  =  1, // char
	NM_INT16				=  2, // short
	NM_INT32				=  3, // int
	NM_INT64				=  4, // long
	NM_FLOAT32			=  5, // float
	NM_FLOAT64			=  6, // double
	NM_COMPLEX64		=  7, // Complex64 class
	NM_COMPLEX128	  =  8, // Complex128 class
	NM_RATIONAL32	  =  9, // Rational32 class
	NM_RATIONAL64	  = 10, // Rational64 class
	NM_RATIONAL128	= 11, // Rational128 class
	NM_RUBYOBJ			= 12  // Ruby VALUE type
} nm_dtype_t;


typedef struct nm_storage {
	// Common elements found in all storage types. Should not be re-arranged.
	nm_dtype_t	dtype;
	size_t	rank;
	size_t*	shape;
	size_t*	offset;

	//virtual void empty(void) = 0;
} STORAGE;

typedef struct nm_dense_storage {
	nm_dtype_t	dtype;
	size_t	rank;
	size_t*	shape;
	size_t*	offset;

	size_t*	stride;
	int			count;
	void*		src;
	void*		elements;
} DENSE_STORAGE;

typedef struct nm_yale_storage {
	nm_dtype_t	dtype;
	size_t	rank;
	size_t*	shape;
	size_t*	offset;

	// Yale storage specific elements.
	void* a;

	// Strictly non-diagonal non-zero count!
	size_t ndnz;

	size_t	capacity;
	itype_t	itype;
	void*		ija;
} YALE_STORAGE;


struct nm_list_node {
  size_t key;
  void*  val;
  struct nm_list_node* next;
};


struct nm_list {
  struct nm_list_node* first;
};


typedef struct nm_list_storage {
	// List storage specific elements.
	void* default_val;
	struct nm_list* rows;
} LIST_STORAGE;


typedef enum {
  NM_DENSE_STORE,
  NM_LIST_STORE,
  NM_YALE_STORE
} nm_stype_t;


typedef struct nmatrix {
	// Method of storage (csc, dense, etc).
	nm_stype_t		stype;
	// Pointer to storage struct.
	STORAGE*	storage;
} NMATRIX;

*/
//end masa

// re-define so it doesn't use reinterpret_cast<>
//#define NM_STRUCT(val)          ((NMATRIX*)(DATA_PTR(val)))
#define NM_STORAGE(val)         (NM_STRUCT(val)->storage)
#define NM_LIST_STORAGE(val)    ((LIST_STORAGE*)(NM_STORAGE(val)))
#define NM_YALE_STORAGE(val)    ((YALE_STORAGE*)(NM_STORAGE(val)))
#define NM_DENSE_STORAGE(val)   ((DENSE_STORAGE*)(NM_STORAGE(val)))

//#define NM_DENSE_SRC(val)       (NM_DENSE_STORAGE(val)->src)
#define NM_RANK(val)            (NM_STORAGE(val)->rank)
#define NM_DTYPE(val)           (NM_STORAGE(val)->dtype)
//#define NM_ITYPE(val)           (NM_YALE_STORAGE(val)->itype)
#define NM_STYPE(val)           (NM_STRUCT(val)->stype)
#define NM_SHAPE(val,i)         (NM_STORAGE(val)->shape[(i)])
#define NM_SHAPE0(val)          (NM_STORAGE(val)->shape[0])
#define NM_SHAPE1(val)          (NM_STORAGE(val)->shape[1])


// External API
extern VALUE rb_nmatrix_dense_create(nm_dtype_t dtype, size_t* shape, size_t rank, void* elements, size_t length);
extern VALUE rb_nvector_dense_create(nm_dtype_t dtype, void* elements, size_t length);

extern VALUE nm_eDataTypeError, nm_eStorageTypeError;

/* GSL::Vector -> NMatrix */

static VALUE rb_gsl_vector_to_nvector(VALUE obj, VALUE klass) {
  gsl_vector *v = NULL;
  Data_Get_Struct(obj, gsl_vector, v);

  return rb_nvector_dense_create(NM_FLOAT64, v->data, v->size);
}


static VALUE rb_gsl_vector_complex_to_nvector(VALUE obj, VALUE klass) {
  gsl_vector *v = NULL;
  Data_Get_Struct(obj, gsl_vector, v);

  return rb_nvector_dense_create(NM_COMPLEX128, v->data, v->size);
}


static VALUE rb_gsl_vector_int_to_nvector(VALUE obj, VALUE klass) {
  gsl_vector *v = NULL;
  Data_Get_Struct(obj, gsl_vector, v);

  return rb_nvector_dense_create(NM_INT64, v->data, v->size);
}


static VALUE rb_gsl_matrix_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_matrix *m = NULL;
  Data_Get_Struct(obj, gsl_matrix, m);

  return rb_nmatrix_dense_create(NM_FLOAT64, &(m->size1), 2, m->data, m->size1 * m->size2);
}


static VALUE rb_gsl_matrix_complex_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_matrix *m = NULL;
  Data_Get_Struct(obj, gsl_matrix, m);

  return rb_nmatrix_dense_create(NM_COMPLEX128, &(m->size1), 2, m->data, m->size1 * m->size2);
}


static VALUE rb_gsl_matrix_int_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_matrix *m = NULL;
  Data_Get_Struct(obj, gsl_matrix, m);

  return rb_nmatrix_dense_create(NM_INT64, &(m->size1), 2, m->data, m->size1 * m->size2);
}


gsl_vector* nv_to_gv(VALUE nm) {
  DENSE_STORAGE* s = NM_DENSE_STORAGE(nm);
  gsl_vector* v = gsl_vector_alloc( s->count );

  if (s->dtype != NM_FLOAT64) {
    rb_raise(nm_eDataTypeError, "requires dtype of :float64 to convert to a GSL vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(double));

  return v;
}


gsl_vector_complex* nv_to_gv_complex(VALUE nm) {
  DENSE_STORAGE* s = NM_DENSE_STORAGE(nm);
  gsl_vector_complex* v = gsl_vector_complex_alloc( s->count );

  if (s->dtype != NM_COMPLEX128) {
    rb_raise(nm_eDataTypeError, "requires dtype of :complex128 to convert to a GSL complex vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(double)*2);

  return v;
}


gsl_vector_int* nv_to_gv_int(VALUE nm) {
  DENSE_STORAGE* s = NM_DENSE_STORAGE(nm);
  gsl_vector_int* v = gsl_vector_int_alloc( s->count );

  if (s->dtype != NM_INT64) {
    rb_raise(nm_eDataTypeError, "requires dtype of :int64 to convert to a GSL int vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(int));

  return v;
}


gsl_matrix* nm_to_gm(VALUE nm) {
  DENSE_STORAGE* s = NM_DENSE_STORAGE(nm);
  gsl_matrix* m = gsl_matrix_alloc( s->shape[0], s->shape[1] );

  if (s->dtype != NM_FLOAT64) {
    rb_raise(nm_eDataTypeError, "requires dtype of :float64 to convert from a GSL double vector");
  }

  memcpy(m->data, s->elements, s->count);
  return m;
}

gsl_matrix_complex* nm_to_gm_complex(VALUE nm) {
  DENSE_STORAGE* s = NM_DENSE_STORAGE(nm);
  gsl_matrix_complex* m = gsl_matrix_complex_alloc( s->shape[0], s->shape[1] );

  if (s->dtype != NM_COMPLEX128) {
    rb_raise(nm_eDataTypeError, "requires dtype of :complex128 to convert from a GSL complex vector");
  }

  memcpy(m->data, s->elements, s->count);
  return m;
}


gsl_matrix_int* nm_to_gm_int(VALUE nm) {
  DENSE_STORAGE* s = NM_DENSE_STORAGE(nm);
  gsl_matrix_int* m = gsl_matrix_int_alloc( s->shape[0], s->shape[1] );

  if (s->dtype != NM_INT64) {
    rb_raise(nm_eDataTypeError, "requires dtype of :int64 to convert from a GSL int vector");
  }

  memcpy(m->data, s->elements, s->count);
  return m;
}


VALUE rb_gsl_nv_to_gsl_vector(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, nv_to_gv(n));
}

VALUE rb_gsl_nv_to_gsl_vector_complex(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, nv_to_gv_complex(n));
}

VALUE rb_gsl_nv_to_gsl_vector_int(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, nv_to_gv_int(n));
}

VALUE rb_gsl_nm_to_gsl_matrix(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, nm_to_gm(n));
}

VALUE rb_gsl_nm_to_gsl_matrix_complex(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, nm_to_gm_complex(n));
}

VALUE rb_gsl_nm_to_gsl_matrix_int(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_matrix_int, 0, gsl_matrix_int_free, nm_to_gm_int(n));
}

VALUE rb_gsl_nv_to_gsl_vector_method(VALUE nv) {
  VALUE v;

  if (NM_DTYPE(nv) == NM_COMPLEX64 || NM_DTYPE(nv) == NM_COMPLEX128) {
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, nv_to_gv_complex(nv));
  } else {
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, nv_to_gv(nv));
  }

  return v;
}



void Init_gsl_nmatrix(VALUE module) {
  rb_define_method(cgsl_vector, "to_nv", rb_gsl_vector_to_nvector, 0);
  rb_define_alias(cgsl_vector, "to_nm", "to_nv");

  rb_define_singleton_method(cgsl_vector, "to_gslv",    rb_gsl_nv_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "to_gv",      rb_gsl_nv_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "nv_to_gslv", rb_gsl_nv_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "nv_to_gv",   rb_gsl_nv_to_gsl_vector, 1);

//  rb_define_method(cgsl_vector_complex, "to_nv", rb_gsl_vector_to_nvector_complex, 0);
//  rb_define_alias(cgsl_vector_complex, "to_nm", "to_nv");

  rb_define_singleton_method(cgsl_vector_complex, "to_gslv",    rb_gsl_nv_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "to_gv",      rb_gsl_nv_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "nv_to_gslv", rb_gsl_nv_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "nv_to_gv",   rb_gsl_nv_to_gsl_vector_complex, 1);

//  rb_define_method(cgsl_vector_int, "to_nv", rb_gsl_vector_to_nvector_int, 0);
//  rb_define_alias(cgsl_vector_int, "to_nm", "to_nv");

  rb_define_singleton_method(cgsl_vector_int, "to_gslv",    rb_gsl_nv_to_gsl_vector_int, 1);
  rb_define_singleton_method(cgsl_vector_int, "to_gv",      rb_gsl_nv_to_gsl_vector_int, 1);
  rb_define_singleton_method(cgsl_vector_int, "nv_to_gslv", rb_gsl_nv_to_gsl_vector_int, 1);
  rb_define_singleton_method(cgsl_vector_int, "nv_to_gv",   rb_gsl_nv_to_gsl_vector_int, 1);

  rb_define_method(cgsl_matrix, "to_nm", rb_gsl_matrix_to_nmatrix, 0);
  rb_define_singleton_method(cgsl_matrix, "nm_to_gslm",  rb_gsl_nm_to_gsl_matrix, 1);

  rb_define_method(cgsl_matrix_complex, "to_nm", rb_gsl_matrix_complex_to_nmatrix, 0);
  rb_define_singleton_method(cgsl_matrix_complex, "nm_to_gslm",  rb_gsl_nm_to_gsl_matrix_complex, 1);

  rb_define_method(cgsl_matrix_int, "to_nm", rb_gsl_matrix_int_to_nmatrix, 0);
  rb_define_singleton_method(cgsl_matrix_int, "nm_to_gslm",  rb_gsl_nm_to_gsl_matrix_int, 1);

  rb_define_method(cNMatrix, "to_gslv",    rb_gsl_nv_to_gsl_vector_method, 0);
  rb_define_alias(cNMatrix, "to_gv", "to_gslv");
}

//#endif // HAVE_NMATRIX_H
