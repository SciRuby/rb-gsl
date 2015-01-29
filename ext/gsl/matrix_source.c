/*
  matrix_source.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#ifdef BASE_DOUBLE
#define NUMCONV(x) NUM2DBL(x)
#define NUMCONV2(x) NUM2DBL(x)
#define PRINTF_FORMAT "%4.3e "
#define MAT_ROW_COL MATRIX_ROW_COL
#define MAT_P MATRIX_P
#define MAT_ROW_P MATRIX_ROW_P
#define MAT_COL_P MATRIX_COL_P
#define C_TO_VALUE rb_float_new
#define C_TO_VALUE2 rb_float_new
#define CHECK_MAT CHECK_MATRIX
#define MAT_VIEW_P MATRIX_VIEW_P
#define VEC_ROW_COL VECTOR_ROW_COL
#define VEC_P VECTOR_P
#define VEC_ROW_P VECTOR_ROW_P
#define VEC_COL_P VECTOR_COL_P
#define CHECK_VEC CHECK_VECTOR
#define VEC_VIEW_P VECTOR_VIEW_P
#elif defined(BASE_INT)
#define NUMCONV(x) FIX2INT(x)
#define NUMCONV2(x) NUM2INT(x)
#define PRINTF_FORMAT "%d "
#define MAT_ROW_COL MATRIX_INT_ROW_COL
#define MAT_P MATRIX_INT_P
#define C_TO_VALUE INT2FIX
#define C_TO_VALUE2 INT2NUM
#define MAT_ROW_P MATRIX_INT_ROW_P
#define MAT_COL_P MATRIX_INT_COL_P
#define CHECK_MAT CHECK_MATRIX_INT
#define MAT_VIEW_P MATRIX_INT_VIEW_P
#define VEC_ROW_COL VECTOR_INT_ROW_COL
#define VEC_P VECTOR_INT_P
#define VEC_ROW_P VECTOR_INT_ROW_P
#define VEC_COL_P VECTOR_INT_COL_P
#define CHECK_VEC CHECK_VECTOR_INT
#define VEC_VIEW_P VECTOR_INT_VIEW_P
#endif

// From ext/vector_source.c
void FUNCTION(get_range,beg_en_n)(VALUE range, BASE *beg, BASE *en, size_t *n, int *step);

// From ext/vector_source.c
void get_range_beg_en_n_for_size(VALUE range,
    int *beg, int *en, size_t *n, int *step, size_t size);

void parse_submatrix_args(int argc, VALUE *argv, size_t size1, size_t size2,
    size_t *i, size_t *j, size_t *n1, size_t *n2);

#ifdef BASE_DOUBLE
void parse_submatrix_args(int argc, VALUE *argv, size_t size1, size_t size2,
    size_t *i, size_t *j, size_t *n1, size_t *n2)
{
  int ii, ij, in1, in2, end, step;

  switch (argc) {
  // no args -> same as submatrix(0, 0, size1, size2)
  case 0:
    *i = 0; *j = 0;
    *n1 = size1; *n2 = size2;
    break;

  // Fixnum -> Same as submatrix(i/size2, i%size2, 1, 1)
  case 1:
    CHECK_FIXNUM(argv[0]);
    ii = FIX2INT(argv[0]);
    *n1 = size1 * size2;
    if(ii < 0) ii += *n1;
    // TODO Bounds check?
    *i = ii / size2; *j = ii % size2;
    *n1 = 1; *n2 = 1;
    break;


  // nil, nil -> all rows, all cols (Matrix::View)
  // nil, Range -> all rows, Range cols (Matrix::View)
  // nil, Fixnum -> all rows, single col (Vector::Col::View)
  // Range, nil -> Range rows, all cols (Matrix::View)
  // Range, Range -> Range rows, Range cols (Matrix::View)
  // Range, Fixnum -> Range rows, single col (Vector::Col::View)
  // Fixnum, nil -> single row, all cols (Vector::View)
  // Fixnum, Range -> single row, Range cols (Vector::View)
  // Fixnum, Fixnum -> single row, single col (Matrix::View)
  case 2: 
    // nil, ...
    if(NIL_P(argv[0])) {
      // Parse second arg
      if(NIL_P(argv[1])) {
        // nil, nil -> all rows, all cols (Matrix::View)
        *i = 0; *j = 0;
        *n1 = size1; *n2 = size2;
      } else if(rb_obj_is_kind_of(argv[1], rb_cRange)) {
        // nil, Range -> all rows, Range cols (Matrix::View)
        *i = 0; *n1 = size1;
        get_range_beg_en_n_for_size(argv[1], &ij, &end, n2, &step, size2);
        if(step < 0 || *n2 <=0) {
          rb_raise(rb_eRangeError, "begin > end");
        }
        *j = (size_t)ij;
      } else {
        // nil, Fixnum -> all rows, single col (Vector::Col::View)
        ij = NUM2INT(argv[1]);
        if(ij < 0) ij += size2;
        *i = 0; *j = (size_t)ij;
        *n1 = size1; *n2 = 0; // *n2 == 0 tells #submatrix to return Vector::Col::View
      }
    // Range, ...
    } else if(rb_obj_is_kind_of(argv[0], rb_cRange)) {
      get_range_beg_en_n_for_size(argv[0], &ii, &end, n1, &step, size1);
      if(step < 0 || *n1 <= 0) {
        rb_raise(rb_eRangeError, "arg0: begin > end");
      }
      *i = (size_t)ii;
      // Parse second arg
      if(NIL_P(argv[1])) {
        // Range, nil -> Range rows, all cols (Matrix::View)
        *j = 0; *n2 = size2;
      } else if(rb_obj_is_kind_of(argv[1], rb_cRange)) {
        // Range, Range -> Range rows, Range cols (Matrix::View)
        get_range_beg_en_n_for_size(argv[1], &ij, &end, n2, &step, size2);
        if(step < 0 || *n2 <= 0) {
          rb_raise(rb_eRangeError, "arg1: begin > end");
        }
        *j = (size_t)ij;
      } else {
        // Range, Fixnum -> Range rows, single col (Vector::Col::View)
        ij = NUM2INT(argv[1]);
        if(ij < 0) ij += size2;
        *j = (size_t) ij;
        *n2 = 0; // *n2 == 0 tells #submatrix to return Vector::Col::View
      }
    // Fixnum, ...
    } else {
      ii = NUM2INT(argv[0]);
      if(ii < 0) ii += size1;
      if(NIL_P(argv[1])) {
        // Fixnum, nil -> single row, all cols (Vector::View)
        *i = (size_t)ii; *j = 0;
        *n1 = 0; *n2 = size2; // *n1 == 0 tells #submatrix to return Vector::View
      } else if(rb_obj_is_kind_of(argv[1], rb_cRange)) {
        // Fixnum, Range -> single row, Range cols (Vector::View)
        get_range_beg_en_n_for_size(argv[1], &ij, &end, n2, &step, size2);
        if(step < 0 || *n2 <= 0) {
          rb_raise(rb_eRangeError, "arg1: begin > end");
        }
        *i = (size_t)ii;
        *j = (size_t)ij;
        *n1 = 0; // *n1 == 0 tells #submatrix to return Vector::View
      } else {
        // Fixnum, Fixnum -> single row, single col (Matrix::View)
        ij = NUM2INT(argv[1]);
        if(ij < 0) ij += size2;
        *i = (size_t)ii; *j = (size_t)ij;
        *n1 = 1; *n2 = 1;
      }
    }
    break;

  // nil, Fixnum, Fixnum -> All rows, some cols
  // Range, Fixnum, Fixnum -> Range rows, some cols
  // Fixnum, Fixnum, nil -> Some rows, all cols
  // Fixnum, Fixnum, Range -> Some rows, Range cols
  case 3:
    // nil, Fixnum, Fixnum
    if(NIL_P(argv[0])) {
      // nil, Fixnum, Fixnum -> All rows, some cols
      CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
      *i = 0;
      ij = FIX2INT(argv[1]);
      *n1 = size1;
      in2 = FIX2INT(argv[2]);
      if(ij < 0) ij += size2;
      *j = (size_t) ij;
      *n2 = (size_t) in2;
    // Range, Fixnum, Fixnum
    } else if(rb_obj_is_kind_of(argv[0], rb_cRange)) {
      // Range, Fixnum, Fixnum -> Range rows, some cols
      CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
      get_range_beg_en_n_for_size(argv[0], &ii, &end, n1, &step, size1);
      if(step < 0 || *n1 <= 0) {
        rb_raise(rb_eRangeError, "arg0: begin > end");
      }
      ij = FIX2INT(argv[1]);
      in2 = FIX2INT(argv[2]);
      if(ij < 0) ij += size2;
      *i = (size_t)ii;
      *j = (size_t) ij;
      *n2 = (size_t) in2;
    // Fixnum, Fixnum, ...
    } else {
      CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
      ii = FIX2INT(argv[0]);
      if(ii < 0) ii += size1; 
      in1 = FIX2INT(argv[1]);
      *i = (size_t)ii;
      *n1 = (size_t)in1;
      // Parse arg2
      if(NIL_P(argv[2])) {
        // Fixnum, Fixnum, nil -> Some rows, all cols
        *j = 0;
        *n2 = size2;
      } else if(rb_obj_is_kind_of(argv[2], rb_cRange)) {
        // Fixnum, Fixnum, Range -> Some rows, Range cols
        get_range_beg_en_n_for_size(argv[2], &ij, &end, n2, &step, size2);
        if(step < 0 || *n2 <= 0) {
          rb_raise(rb_eRangeError, "arg2: begin > end");
        }
        *j = (size_t)ij;
      } else {
        rb_raise(rb_eArgError, "expected third argument to be nil or Range, not %s",
           rb_class2name(CLASS_OF(argv[2])));
      }
    }
    break;

  case 4:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    CHECK_FIXNUM(argv[2]); CHECK_FIXNUM(argv[3]);
    ii = FIX2INT(argv[0]);  ij = FIX2INT(argv[1]);
    in1 = FIX2INT(argv[2]);  in2 = FIX2INT(argv[3]);
    if(ii < 0) ii += size1;
    if(ij < 0) ij += size2;
    // TODO Bounds check?
    *i = (size_t)ii; *j = (size_t)ij;
    *n1 = (size_t)in1; *n2 = (size_t)in2;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 to 4)", argc);
    break;
  }
}
#endif

VALUE FUNCTION(rb_gsl_matrix,do_something)(VALUE obj, void (*f)(GSL_TYPE(gsl_matrix) *))
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  (*f)(m);
  return obj;
}

static VALUE FUNCTION(create_matrix,from_range_shape)(VALUE range, VALUE nn1, VALUE nn2);
static VALUE FUNCTION(create_matrix,from_ranges)(int argc, VALUE *argv);
GSL_TYPE(gsl_matrix)* FUNCTION(gsl_matrix,alloc_from_colvectors)(int argc, VALUE *argv);

static VALUE FUNCTION(rb_gsl_matrix,alloc)(int argc, VALUE *argv, VALUE klass)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t n1, n2;
#ifdef HAVE_NARRAY_H
  size_t n;
  VALUE ary;
  struct NARRAY *na;
#endif
  if (argc < 1) rb_raise(rb_eArgError, 
       "wrong number of arguments (%d for >= 1)", argc);

#ifdef HAVE_NARRAY_H
  if (NA_IsNArray(argv[0])) {
    GetNArray(argv[0], na);
    n = na->shape[0]*na->shape[1];
    m = FUNCTION(gsl_matrix,alloc)(na->shape[1], na->shape[0]);
    if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
#ifdef BASE_DOUBLE
    ary = na_change_type(argv[0], NA_DFLOAT);
#else
    ary = na_change_type(argv[0], NA_LINT);
#endif
    memcpy(m->data, NA_PTR_TYPE(ary,BASE*), n*sizeof(BASE));
    return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_matrix,free), m);
  }
#endif
  switch (TYPE(argv[0])) {
  case T_FIXNUM:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)", 
          argc);
    CHECK_FIXNUM(argv[1]);
    n1 = FIX2INT(argv[0]); n2 = FIX2INT(argv[1]);
    m = FUNCTION(gsl_matrix,calloc)(n1, n2);
    break;
  case T_ARRAY:
    if (argc == 1) {
      m = FUNCTION(gsl_matrix,alloc_from_arrays)(argc, argv);
      break;
    }
    if (CLASS_OF(argv[1]) == rb_cRange) argv[1] = rb_gsl_range2ary(argv[1]);
    switch (TYPE(argv[1])) {
    case T_ARRAY:
      m = FUNCTION(gsl_matrix,alloc_from_arrays)(argc, argv);
      break;
    case T_FIXNUM:
      if (argc != 3) rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)",
            argc);
      CHECK_FIXNUM(argv[2]);
      m = FUNCTION(gsl_matrix,alloc_from_array_sizes)(argv[0], argv[1], argv[2]);
      break;
    default:
      rb_raise(rb_eTypeError, 
         "wrong argument type %s\nUsage: new(n1, n2), "
         "new([], [], [], ...), new([], n1, n2)", 
         rb_class2name(CLASS_OF(argv[1])));
      break;
    }
    break;
  default:
    if (CLASS_OF(argv[0]) == rb_cRange) {
      if (argc==3 && TYPE(argv[1]) == T_FIXNUM && TYPE(argv[2]) == T_FIXNUM)
  return FUNCTION(create_matrix,from_range_shape)(argv[0], argv[1], argv[2]);
      else
  return FUNCTION(create_matrix,from_ranges)(argc, argv);
    } else if (VEC_P(argv[0])) {
      if (argc == 3 && FIXNUM_P(argv[1]) && FIXNUM_P(argv[2])) {
  m = FUNCTION(gsl_matrix,alloc_from_vector_sizes)(argv[0], argv[1], argv[2]);
      } else {
  if (VEC_COL_P(argv[0])) {
    m = FUNCTION(gsl_matrix,alloc_from_colvectors)(argc, argv);
  } else {
    m = FUNCTION(gsl_matrix,alloc_from_vectors)(argc, argv);
  }
      }
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s\n"
         "Usage: new(n1, n2), new([], [], [], ...), new([], n1, n2)", 
         rb_class2name(CLASS_OF(argv[0])));
    }
    break;
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

void FUNCTION(set_ptr_data,by_range)(BASE *ptr, size_t n, VALUE range);

static GSL_TYPE(gsl_matrix)* FUNCTION(cr_matrix,from_ranges)(int argc, VALUE *argv)
{
  GSL_TYPE(gsl_matrix) *m;
  BASE beg, en;
  size_t i, n;
  int step;
  FUNCTION(get_range,beg_en_n)(argv[0], &beg, &en, &n, &step);
  m = FUNCTION(gsl_matrix,calloc)(argc, n);
  FUNCTION(set_ptr_data,by_range)(m->data, n, argv[0]);
  for (i = 1; (int) i < argc; i++) {
    if (CLASS_OF(argv[i]) != rb_cRange) 
      rb_raise(rb_eTypeError, "wrong argument type %s (Range expected)",
         rb_class2name(CLASS_OF(argv[i])));
    FUNCTION(set_ptr_data,by_range)(m->data+i*n, n, argv[i]);
  }
  return m;
}

static VALUE FUNCTION(create_matrix,from_ranges)(int argc, VALUE *argv)
{
  GSL_TYPE(gsl_matrix) *m;
  m = FUNCTION(cr_matrix,from_ranges)(argc, argv);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

static GSL_TYPE(gsl_matrix)* FUNCTION(cr_matrix,from_range_shape)(VALUE range, VALUE nn1, VALUE nn2)
{
  size_t n1, n2;
  GSL_TYPE(gsl_matrix) *m;
  n1 = FIX2INT(nn1);  n2 = FIX2INT(nn2);
  m = FUNCTION(gsl_matrix,alloc)(n1, n2);
  FUNCTION(set_ptr_data,by_range)(m->data, n1*n2, range);
  return m;
}

static VALUE FUNCTION(create_matrix,from_range_shape)(VALUE range, VALUE nn1, VALUE nn2)
{
  GSL_TYPE(gsl_matrix) *m;
  m = FUNCTION(cr_matrix,from_range_shape)(range, nn1, nn2);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

/* create a matrix from arrays */
GSL_TYPE(gsl_matrix)* FUNCTION(gsl_matrix,alloc_from_arrays)(int argc, VALUE *argv)
{
  size_t n, i, j;
  GSL_TYPE(gsl_matrix) *m = NULL;
  if (CLASS_OF(argv[0]) == rb_cRange) argv[0] = rb_gsl_range2ary(argv[0]);
  else Check_Type(argv[0], T_ARRAY);
  n = RARRAY_LEN(argv[0]);
  m = FUNCTION(gsl_matrix,alloc)(argc, n);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
  for (i = 0; (int) i < argc; i++) {
    if (CLASS_OF(argv[i]) == rb_cRange) argv[i] = rb_gsl_range2ary(argv[i]);
    else Check_Type(argv[i], T_ARRAY);
    for (j = 0; j < n; j++) {
      if ((int) j >= RARRAY_LEN(argv[i])) FUNCTION(gsl_matrix,set)(m, i, j, 0);
      else FUNCTION(gsl_matrix,set)(m, i, j, NUMCONV2(rb_ary_entry(argv[i], j)));
    }
  }
  return m;
}

GSL_TYPE(gsl_matrix)* FUNCTION(gsl_matrix,alloc_from_vectors)(int argc, VALUE *argv)
{
  size_t i;
  GSL_TYPE(gsl_matrix) *m = NULL;
  GSL_TYPE(gsl_vector) *v = NULL;
  if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
  CHECK_VEC(argv[0]);
  Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), v);
  m = FUNCTION(gsl_matrix,alloc)(argc, v->size);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
  for (i = 0; (int) i < argc; i++) {
    CHECK_VEC(argv[i]);
    Data_Get_Struct(argv[i], GSL_TYPE(gsl_vector), v);
    FUNCTION(gsl_matrix,set_row)(m, i, v);
  }
  return m;
}

GSL_TYPE(gsl_matrix)* FUNCTION(gsl_matrix,alloc_from_colvectors)(int argc, VALUE *argv)
{
  size_t i;
  GSL_TYPE(gsl_matrix) *m = NULL;
  GSL_TYPE(gsl_vector) *v = NULL;
  if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
  CHECK_VEC(argv[0]);
  Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), v);
//  m = FUNCTION(gsl_matrix,alloc)(argc, v->size);
  m = FUNCTION(gsl_matrix,alloc)(v->size, argc);  
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
  for (i = 0; (int) i < argc; i++) {
    CHECK_VEC(argv[i]);
    Data_Get_Struct(argv[i], GSL_TYPE(gsl_vector), v);
    FUNCTION(gsl_matrix,set_col)(m, i, v);
  }
  return m;
}

/* create a matrix from two sizes and an array */
GSL_TYPE(gsl_matrix)* FUNCTION(gsl_matrix,alloc_from_array_sizes)(VALUE ary,
                  VALUE nn1, VALUE nn2)
{
  size_t n1, n2, len;
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t i, j, k;
  CHECK_FIXNUM(nn1); CHECK_FIXNUM(nn2);
  Check_Type(ary, T_ARRAY);
  n1 = FIX2INT(nn1);
  n2 = FIX2INT(nn2);
  m = FUNCTION(gsl_matrix,alloc)(n1, n2);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
  k = 0;
  len = RARRAY_LEN(ary);
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++, k++) {
      if (k >= len)
  FUNCTION(gsl_matrix,set)(m, i, j, (BASE) 0);
      else
  FUNCTION(gsl_matrix,set)(m, i, j, (BASE) NUMCONV2(rb_ary_entry(ary, k)));
      
    }
  }
  return m;
}

GSL_TYPE(gsl_matrix)* FUNCTION(gsl_matrix,alloc_from_vector_sizes)(VALUE ary,
                   VALUE nn1, VALUE nn2)
{
  size_t n1, n2;
  GSL_TYPE(gsl_matrix) *m = NULL;
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t i, j, k;
  CHECK_VEC(ary);
  CHECK_FIXNUM(nn1); CHECK_FIXNUM(nn2);
  Data_Get_Struct(ary, GSL_TYPE(gsl_vector), v);
  n1 = FIX2INT(nn1); n2 = FIX2INT(nn2);
  m = FUNCTION(gsl_matrix,alloc)(n1, n2);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
  k = 0;
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++, k++) {
      if (k >= v->size) 
  FUNCTION(gsl_matrix,set)(m, i, j, (BASE) 0);
      else
  FUNCTION(gsl_matrix,set)(m, i, j, FUNCTION(gsl_vector,get)(v, k));
    }
  }
  return m;
}

static VALUE FUNCTION(rb_gsl_matrix,calloc)(VALUE klass, VALUE nn1, VALUE nn2)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  CHECK_FIXNUM(nn1); CHECK_FIXNUM(nn2);
  m = FUNCTION(gsl_matrix,calloc)(FIX2INT(nn1), FIX2INT(nn2));
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_calloc failed");
  return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_matrix,free), m);
}

static VALUE FUNCTION(rb_gsl_matrix,diagonal_singleton)(int argc, VALUE *argv, VALUE klass)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  GSL_TYPE(gsl_vector) *v = NULL;
  VALUE ary, tmp;
  size_t len, i;
  switch (argc) {
  case 1:
    switch (TYPE(argv[0])) {
    case T_FIXNUM:
    case T_FLOAT:
      len = FIX2INT(argv[0]);
      m = FUNCTION(gsl_matrix,alloc)(len, len);
      for (i = 0; i < len; i++)
  FUNCTION(gsl_matrix,set)(m, i, i, 1);
      return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_matrix,free), m);
      break;
    default:
      /* do next */
      break;
    }
    if (rb_obj_is_kind_of(argv[0], rb_cRange))
      ary = rb_gsl_range2ary(argv[0]);
    else ary = argv[0];
    switch (TYPE(ary)) {
    case T_ARRAY:
      len = RARRAY_LEN(ary);
      m = FUNCTION(gsl_matrix,calloc)(len, len);
      for (i = 0; i < len; i++) {
  tmp = rb_ary_entry(ary, i);
  FUNCTION(gsl_matrix,set)(m, i, i, NUMCONV2(tmp));
      }
      break;
    default:
      CHECK_VEC(ary);
      Data_Get_Struct(ary, GSL_TYPE(gsl_vector), v);
      len = v->size;
      m = FUNCTION(gsl_matrix,calloc)(len, len);
      for (i = 0; i < len; i++) {
  FUNCTION(gsl_matrix,set)(m, i, i, FUNCTION(gsl_vector,get)(v, i));
      }
      break;
    }
    break;
  default:
    m = FUNCTION(gsl_matrix,calloc)(argc, argc);
    for (i = 0; (int) i < argc; i++) {
      FUNCTION(gsl_matrix,set)(m, i, i, NUMCONV2(argv[i]));
    }
    break;
  }
  return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_matrix,free), m);
}

static VALUE FUNCTION(rb_gsl_matrix,eye)(int argc, VALUE *argv, VALUE klass)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t n1, n2, n, i;
  switch (argc) {
  case 1:
    CHECK_FIXNUM(argv[0]);
    n = FIX2INT(argv[0]);
    n1 = n2 = n;
    break;
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    n1 = FIX2INT(argv[0]); n2 = FIX2INT(argv[1]);
    n = GSL_MIN_INT(n1, n2);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  m = FUNCTION(gsl_matrix,calloc)(n1, n2);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_calloc failed");
  for (i = 0; i < n; i++) {
    FUNCTION(gsl_matrix,set)(m, i, i, 1);
  }
  return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_matrix,free), m);
}

static VALUE FUNCTION(rb_gsl_matrix,ones)(int argc, VALUE *argv, VALUE klass)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t n1, n2, n, i, j;
  switch (argc) {
  case 1:
    CHECK_FIXNUM(argv[0]);
    n = FIX2INT(argv[0]);
    n1 = n2 = n;
    break;
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    n1 = FIX2INT(argv[0]); n2 = FIX2INT(argv[1]);
    n = GSL_MIN_INT(n1, n2);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  m = FUNCTION(gsl_matrix,calloc)(n1, n2);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_calloc failed");
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) {
      FUNCTION(gsl_matrix,set)(m, i, j, (BASE) 1);
    }
  }
  return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_matrix,free), m);
}

static VALUE FUNCTION(rb_gsl_matrix,zeros)(int argc, VALUE *argv, VALUE klass)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t n1, n2, n, i, j;
  switch (argc) {
  case 1:
    CHECK_FIXNUM(argv[0]);
    n = FIX2INT(argv[0]);
    n1 = n2 = n;
    break;
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    n1 = FIX2INT(argv[0]); n2 = FIX2INT(argv[1]);
    n = GSL_MIN_INT(n1, n2);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  m = FUNCTION(gsl_matrix,calloc)(n1, n2);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_calloc failed");
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) {
      FUNCTION(gsl_matrix,set)(m, i, j, (BASE) 0);
    }
  }
  return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_matrix,free), m);
}

static VALUE FUNCTION(rb_gsl_matrix,identity)(VALUE klass, VALUE nn)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t n, i;
  CHECK_FIXNUM(nn);
  n = FIX2INT(nn);
  m = FUNCTION(gsl_matrix,calloc)(n, n);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_calloc failed");
  for (i = 0; i < n; i++) {
    FUNCTION(gsl_matrix,set)(m, i, i, 1);
  }
  return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_matrix,free), m);
}

static VALUE FUNCTION(rb_gsl_matrix,size1)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  return INT2FIX(m->size1);
}

static VALUE FUNCTION(rb_gsl_matrix,size2)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  return INT2FIX(m->size2);
}

static VALUE FUNCTION(rb_gsl_matrix,shape)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  return rb_ary_new3(2, INT2FIX(m->size1), INT2FIX(m->size2));
}

static VALUE FUNCTION(rb_gsl_matrix,submatrix)(int argc, VALUE *argv, VALUE obj);
static VALUE FUNCTION(rb_gsl_matrix,get)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  VALUE retval;
  int ii, ij;

  if(argc == 2 && TYPE(argv[0]) == T_FIXNUM && TYPE(argv[1]) == T_FIXNUM) {
    // m[i,j]
    Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
    ii = FIX2INT(argv[0]);
    ij = FIX2INT(argv[1]);
    if(ii < 0) ii += m->size1;
    if(ij < 0) ij += m->size2;
    retval = C_TO_VALUE2(FUNCTION(gsl_matrix,get)(m, (size_t)ii, (size_t)ij));
  } else if(argc == 1 && TYPE(argv[0]) == T_FIXNUM) {
    // m[i]
    Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
    ii = FIX2INT(argv[0]);
    if(ii < 0) ii += m->size1 * m->size2;
    retval = C_TO_VALUE2(FUNCTION(gsl_matrix,get)(m, (size_t)(ii / m->size2), (size_t)(ii % m->size2)));
  } else if(argc == 1 && TYPE(argv[0]) == T_ARRAY) {
    // m[[i,j]], to support m[m.max_index]
    if(RARRAY_LEN(argv[0]) == 2) {
      Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
      ii = FIX2INT(RARRAY_PTR(argv[0])[0]);
      ij = FIX2INT(RARRAY_PTR(argv[0])[1]);
      if(ii < 0) ii += m->size1;
      if(ij < 0) ij += m->size2;
      retval = C_TO_VALUE2(FUNCTION(gsl_matrix,get)(m, (size_t)ii, (size_t)ij));
    } else {
      rb_raise(rb_eArgError, "Array index must have length 2, not %d", (int) RARRAY_LEN(argv[0]));
    }
  } else {
    retval = FUNCTION(rb_gsl_matrix,submatrix)(argc, argv, obj);
  }

  return retval;
}

void FUNCTION(rb_gsl_vector,set_subvector)(int argc, VALUE *argv, GSL_TYPE(gsl_vector) *v, VALUE other);
static VALUE FUNCTION(rb_gsl_matrix,set)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m, *mother;
  QUALIFIED_VIEW(gsl_matrix,view) mv;
  QUALIFIED_VIEW(gsl_vector,view) vv;
  VALUE other, row, row_set_argv[2];
  int ii, ij, step;
  size_t i, j, k, n1, n2, nother;
  BASE beg, end;

  if(argc < 1 || argc > 5) {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1-5)", argc);
  }

  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  other = argv[argc-1];

  if(argc == 1 && TYPE(argv[0]) == T_ARRAY) {
    // m.set([row0,row1,...])
    n1 = RARRAY_LEN(argv[0]);
    if(n1 > m->size1) n1 = m->size1;
    row_set_argv[0] = INT2FIX(0);
    // Each given row must have as manay elements as m has columns.
    // The bounds check happens inside rb_gsl_vector*set_subvector().
    row_set_argv[1] = INT2FIX(m->size2);
    for(k = 0; k < n1 && k < m->size1; k++) {
      vv = FUNCTION(gsl_matrix,row)(m, k);
      FUNCTION(rb_gsl_vector,set_subvector)(2, row_set_argv, &vv.vector, rb_ary_entry(argv[0], k));
    }
  } else if(argc == 1) {
    // m[] = x
    FUNCTION(gsl_matrix,set_all)(m, NUMCONV2(other));
  } else if(argc==2 && TYPE(argv[0]) == T_ARRAY && TYPE(argv[1]) != T_ARRAY) {
    // m.set([i, j], x) or m[[i,j]] = x
    ii = FIX2INT(rb_ary_entry(argv[0], 0));
    ij = FIX2INT(rb_ary_entry(argv[0], 1));
    if(ii < 0) ii += m->size1;
    if(ij < 0) ij += m->size2;
    FUNCTION(gsl_matrix,set)(m, (size_t)ii, (size_t)ij, NUMCONV2(argv[1]));
  } else if(argc == 3 && TYPE(argv[0]) == T_FIXNUM && TYPE(argv[1]) == T_FIXNUM) {
    // m[i,j] = x
    ii = FIX2INT(argv[0]);
    ij = FIX2INT(argv[1]);
    if(ii < 0) ii += m->size1;
    if(ij < 0) ij += m->size2;
    FUNCTION(gsl_matrix,set)(m, (size_t)ii, (size_t)ij, NUMCONV2(other));
  } else if(TYPE(argv[0]) == T_ARRAY) {
    // m.set(row0,row1,...)
    n1 = argc;
    if(n1 > m->size1) n1 = m->size1;
    row_set_argv[0] = INT2FIX(0);
    row_set_argv[1] = INT2FIX(m->size2);
    for(k = 0; k < n1 && k < m->size1; k++) {
      vv = FUNCTION(gsl_matrix,row)(m, k);
      FUNCTION(rb_gsl_vector,set_subvector)(2, row_set_argv, &vv.vector, argv[k]);
    }
  } else {
    // x -> assignment to m.submatrix(i...)
    parse_submatrix_args(argc-1, argv, m->size1, m->size2, &i, &j, &n1, &n2);
    if(n1 == 0) n1 = 1;
    if(n2 == 0) n2 = 1;
    mv = FUNCTION(gsl_matrix,submatrix)(m, i, j, n1, n2);
    if(rb_obj_is_kind_of(other, GSL_TYPE(cgsl_matrix))) {
      // m[...] = m_other
      Data_Get_Struct(other, GSL_TYPE(gsl_matrix), mother);
      if(n1 * n2 != mother->size1 * mother->size2) {
        rb_raise(rb_eRangeError, "sizes do not match (%d x %d != %d x %d)",
     (int) n1, (int) n2, (int) mother->size1, (int) mother->size2);
      }
      // TODO Change to gsl_matrix_memmove if/when GSL has such a function
      // because gsl_matrix_memcpy does not handle overlapping regions (e.g.
      // Views) well.
      FUNCTION(gsl_matrix,memcpy)(&mv.matrix, mother);
    } else if(rb_obj_is_kind_of(other, rb_cArray)) {
      row_set_argv[0] = INT2FIX(0);
      row_set_argv[1] = INT2FIX(n2);

      if(n1 == 1) {
        // m[...] = [col0, ...] # single row
        vv = FUNCTION(gsl_matrix,row)(&mv.matrix, 0);
        FUNCTION(rb_gsl_vector,set_subvector)(2, row_set_argv, &vv.vector, other);
      } else {
        // m[...] = [[row0], [row1], ...] # multiple rows
        if((int) n1 != RARRAY_LEN(other)) {
          rb_raise(rb_eRangeError, "row counts do not match (%d != %d)",
       (int) n1, (int) RARRAY_LEN(other));
        }
        for(k = 0; k < n1; k++) {
          vv = FUNCTION(gsl_matrix,row)(&mv.matrix, k);
          row = rb_ary_entry(other, k);
          FUNCTION(rb_gsl_vector,set_subvector)(2, row_set_argv, &vv.vector, row);
        }
      }
    } else if(rb_obj_is_kind_of(other, rb_cRange)) {
      // m[...] = beg..end
      FUNCTION(get_range,beg_en_n)(other, &beg, &end, &nother, &step);
      if(n1 * n2 != nother) {
        rb_raise(rb_eRangeError, "sizes do not match (%d x %d != %d)", (int) n1, (int) n2, (int) nother);
      }
      for(k = 0; k < nother; k++) {
        FUNCTION(gsl_matrix,set)(&mv.matrix, k / n2, k % n2, beg);
        beg += step;
      }
    } else {
      // m[...] = x
      FUNCTION(gsl_matrix,set_all)(&mv.matrix, NUMCONV2(other));
    }
  }

  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,set_all)(VALUE obj, VALUE x)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,set_all)(m, NUMCONV2(x));
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,set_zero)(VALUE obj)
{
  return FUNCTION(rb_gsl_matrix,do_something)(obj, FUNCTION(gsl_matrix,set_zero));
}

static VALUE FUNCTION(rb_gsl_matrix,set_identity)(VALUE obj)
{
  return FUNCTION(rb_gsl_matrix,do_something)(obj,  FUNCTION(gsl_matrix,set_identity));
}

static VALUE FUNCTION(rb_gsl_matrix,print)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  printf("[ ");
  for (i = 0; i < m->size1; i++) {
    if (i != 0) printf("  ");
    for (j = 0; j < m->size2; j++) {
      printf(PRINTF_FORMAT, FUNCTION(gsl_matrix,get)(m, i, j));
    }
    if (i == m->size1 - 1) printf("]\n");
    else printf("\n");
  }
  return Qnil;
}

#ifdef BASE_DOUBLE
#define SHOW_ELM 6
#else
#define SHOW_ELM 12
#endif

static VALUE FUNCTION(rb_gsl_matrix,to_s)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  char buf[32], format[32], format2[32];
  size_t i, j;
  VALUE str;
  BASE x;
  int dig = 8;
#ifdef BASE_INT
  BASE min;
  BASE max;
#endif
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
#ifdef BASE_INT
  min = FUNCTION(gsl_matrix,min)(m);
  max = gsl_matrix_int_max(m);
  dig = (int) GSL_MAX(fabs(max),fabs(min));
  if (dig > 0) dig = ceil(log10(dig+1e-10));
  else dig = 1;
  if (min < 0) dig += 1;
  sprintf(format, "%%%dd ", (int) dig);
  strcpy(format2, format);
#else
  strcpy(format, PRINTF_FORMAT);
  strcpy(format2, " "PRINTF_FORMAT);
#endif
  str = rb_str_new2("[ ");
  for (i = 0; i < m->size1; i++) {
    if (i != 0) {
      strcpy(buf, "  ");
      rb_str_cat(str, buf, strlen(buf));
    }
    for (j = 0; j < m->size2; j++) {
      x = FUNCTION(gsl_matrix,get)(m, i, j);
      if (x < 0)
  sprintf(buf, format, x); 
      else
  sprintf(buf, format2, x); 
      rb_str_cat(str, buf, strlen(buf));
      if ((int) j >= (55/dig)) {
        strcpy(buf, "... ");
        rb_str_cat(str, buf, strlen(buf));
        break;
      }
    }
    if (i >= 20) {
      strcpy(buf, "\n  ... ]");
      rb_str_cat(str, buf, strlen(buf));
      break;
    }
    if (i == m->size1 - 1) {
      strcpy(buf, "]");
      rb_str_cat(str, buf, strlen(buf));
    } else {
      strcpy(buf, "\n");
      rb_str_cat(str, buf, strlen(buf));
    }
  }
  return str;
}
#undef SHOW_ELM

static VALUE FUNCTION(rb_gsl_matrix,inspect)(VALUE obj)
{
  VALUE str;
  char buf[64];
  sprintf(buf, "%s\n", rb_class2name(CLASS_OF(obj)));
  str = rb_str_new2(buf);
  return rb_str_concat(str, FUNCTION(rb_gsl_matrix,to_s)(obj));
}

#ifdef BASE_DOUBLE
#define PRINTF_FORMAT2 "%g"
#elif defined(BASE_INT)
#define PRINTF_FORMAT2 "%d"
#endif
static VALUE FUNCTION(rb_gsl_matrix,fprintf)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  if (argc != 1 && argc != 2) {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), h);
  fp = rb_gsl_open_writefile(argv[0], &flag);
  if (argc == 2) {
    Check_Type(argv[1], T_STRING);
    status = FUNCTION(gsl_matrix,fprintf)(fp, h, STR2CSTR(argv[1]));
  } else {
    status = FUNCTION(gsl_matrix,fprintf)(fp, h, PRINTF_FORMAT2);
  }
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_gsl_matrix,printf)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *h = NULL;
  int status;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), h);
  if (argc == 1) {
    Check_Type(argv[0], T_STRING);
    status = FUNCTION(gsl_matrix,fprintf)(stdout, h, STR2CSTR(argv[0]));
  } else {
    status = FUNCTION(gsl_matrix,fprintf)(stdout, h, PRINTF_FORMAT2);
  }
  return INT2FIX(status);
}

#undef PRINTF_FORMAT2
static VALUE FUNCTION(rb_gsl_matrix,fscanf)(VALUE obj, VALUE io)
{
  GSL_TYPE(gsl_matrix) *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), h);
  fp = rb_gsl_open_readfile(io, &flag);
  status = FUNCTION(gsl_matrix,fscanf)(fp, h);
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_gsl_matrix,set_diagonal)(VALUE obj, VALUE diag)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  GSL_TYPE(gsl_vector) *v;
  size_t i, len;
  BASE x;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  switch (TYPE(diag)) {
  case T_FIXNUM: case T_BIGNUM: case T_FLOAT:
    x = (BASE) NUMCONV2(diag);
    for (i = 0; i < m->size1; i++) FUNCTION(gsl_matrix,set)(m, i, i, x);
    break;
  case T_ARRAY:
    len = GSL_MIN_INT((int) m->size1, RARRAY_LEN(diag));
    for (i = 0; i < len; i++) {
      FUNCTION(gsl_matrix,set)(m, i, i, NUMCONV2(rb_ary_entry(diag, i)));
    }
    break;
  default:
    if (VEC_P(diag)) {
      Data_Get_Struct(diag, GSL_TYPE(gsl_vector), v);
      len = GSL_MIN_INT(m->size1, v->size);
      for (i = 0; i < len; i++) {
  FUNCTION(gsl_matrix,set)(m, i, i, FUNCTION(gsl_vector,get)(v, i));
      }
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Vector or Array expected)",
         rb_class2name(CLASS_OF(diag)));
    }
    break;
  }
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,get_row)(VALUE obj, VALUE i)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  GSL_TYPE(gsl_vector) *v = NULL;
  CHECK_FIXNUM(i);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  v = FUNCTION(gsl_vector,alloc)(m->size1);
  if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
  FUNCTION(gsl_matrix,get_row)(v, m, FIX2INT(i));
  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), v);
}

static VALUE FUNCTION(rb_gsl_matrix,get_col)(VALUE obj, VALUE i)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  GSL_TYPE(gsl_vector) *v = NULL;
  CHECK_FIXNUM(i);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  v = FUNCTION(gsl_vector,alloc)(m->size2);
  if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
  FUNCTION(gsl_matrix,get_col)(v, m, FIX2INT(i));
  // TODO This is NOT returning a view!  Is there a macro more appropriate than
  // QUALIFIED_VIEW?
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,col), 0, FUNCTION(gsl_vector,free), v);
}

static VALUE FUNCTION(rb_gsl_matrix,set_row)(VALUE obj, VALUE i, VALUE vv)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  GSL_TYPE(gsl_vector) *v = NULL;
  int flag = 0;
  size_t j;
  CHECK_FIXNUM(i);
  if (CLASS_OF(vv) == rb_cRange) vv = rb_gsl_range2ary(vv);
  if (TYPE(vv) == T_ARRAY) {
    v = FUNCTION(gsl_vector,alloc)(RARRAY_LEN(vv));
    for (j = 0; (int) j < RARRAY_LEN(vv); j++) {
      FUNCTION(gsl_vector,set)(v, j, NUMCONV2(rb_ary_entry(vv, j)));
    }
    flag = 1;
  } else {
    CHECK_VEC(vv);
    Data_Get_Struct(vv, GSL_TYPE(gsl_vector), v);
  }
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,set_row)(m, FIX2INT(i), v);
  if (flag == 1) FUNCTION(gsl_vector,free)(v);
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,set_col)(VALUE obj, VALUE j, VALUE vv)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  GSL_TYPE(gsl_vector) *v = NULL;
  int flag = 0;
  size_t i;
  CHECK_FIXNUM(j);
  if (CLASS_OF(vv) == rb_cRange) vv = rb_gsl_range2ary(vv);
  if (TYPE(vv) == T_ARRAY) {
    v = FUNCTION(gsl_vector,alloc)(RARRAY_LEN(vv));
    for (i = 0; (int) i < RARRAY_LEN(vv); i++) {
      FUNCTION(gsl_vector,set)(v, i, NUMCONV2(rb_ary_entry(vv, i)));
    }
    flag = 1;
  } else {
    CHECK_VECTOR(vv);
    Data_Get_Struct(vv, GSL_TYPE(gsl_vector), v);
  }
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,set_col)(m, FIX2INT(j), v);
  if (flag == 1) FUNCTION(gsl_vector,free)(v);
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,clone)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  FUNCTION(gsl_matrix,memcpy)(mnew, m);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,memcpy)(VALUE obj, VALUE mm1, VALUE mm2)
{
  GSL_TYPE(gsl_matrix) *m1 = NULL, *m2 = NULL;
  CHECK_MAT(mm1);  CHECK_MAT(mm2);
  Data_Get_Struct(mm1, GSL_TYPE(gsl_matrix), m1);
  Data_Get_Struct(mm2, GSL_TYPE(gsl_matrix), m2);
  FUNCTION(gsl_matrix,memcpy)(m1, m2);
  return mm1;
}

static VALUE FUNCTION(rb_gsl_matrix,isnull)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  return INT2FIX(FUNCTION(gsl_matrix,isnull)(m));
}

static VALUE FUNCTION(rb_gsl_matrix,isnull2)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  if (FUNCTION(gsl_matrix,isnull)(m)) return Qtrue;
  else return Qfalse;
}


/* singleton */
static VALUE FUNCTION(rb_gsl_matrix,swap)(VALUE obj, VALUE mm1, VALUE mm2)
{
  GSL_TYPE(gsl_matrix) *m1 = NULL, *m2 = NULL;
  CHECK_MAT(mm1);  CHECK_MAT(mm2);
  Data_Get_Struct(mm1, GSL_TYPE(gsl_matrix), m1);
  Data_Get_Struct(mm2, GSL_TYPE(gsl_matrix), m2);
  FUNCTION(gsl_matrix,swap)(m1, m2);
  return mm1;
}

static VALUE FUNCTION(rb_gsl_matrix,swap_rows_bang)(VALUE obj, VALUE i, VALUE j)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  CHECK_FIXNUM(i);   CHECK_FIXNUM(j);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,swap_rows)(m, FIX2INT(i), FIX2INT(j));
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,swap_rows)(VALUE obj, VALUE i, VALUE j)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mnew;
  CHECK_FIXNUM(i);   CHECK_FIXNUM(j);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(make_matrix,clone)(m);
  FUNCTION(gsl_matrix,swap_rows)(mnew, FIX2INT(i), FIX2INT(j));
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,swap_columns_bang)(VALUE obj, VALUE i, VALUE j)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  CHECK_FIXNUM(i);   CHECK_FIXNUM(j);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,swap_columns)(m, FIX2INT(i), FIX2INT(j));
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,swap_columns)(VALUE obj, VALUE i, VALUE j)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mnew;
  CHECK_FIXNUM(i);   CHECK_FIXNUM(j);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(make_matrix,clone)(m);
  FUNCTION(gsl_matrix,swap_columns)(mnew, FIX2INT(i), FIX2INT(j));
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,swap_rowcol_bang)(VALUE obj, VALUE i, VALUE j)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  CHECK_FIXNUM(i);   CHECK_FIXNUM(j);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,swap_rowcol)(m, FIX2INT(i), FIX2INT(j));
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,swap_rowcol)(VALUE obj, VALUE i, VALUE j)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mnew = NULL;
  CHECK_FIXNUM(i);   CHECK_FIXNUM(j);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(make_matrix,clone)(m);
  FUNCTION(gsl_matrix,swap_rowcol)(mnew, FIX2INT(i), FIX2INT(j));
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,transpose_memcpy)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size2, m->size1);
  FUNCTION(gsl_matrix,transpose_memcpy)(mnew, m);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,transpose_bang)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,transpose)(m);
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,max)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  return C_TO_VALUE2(FUNCTION(gsl_matrix,max)(m));
}

static VALUE FUNCTION(rb_gsl_matrix,min)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  return C_TO_VALUE2(FUNCTION(gsl_matrix,min)(m));
}

static VALUE FUNCTION(rb_gsl_matrix,minmax)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  BASE min, max;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,minmax)(m, &min, &max);
  return rb_ary_new3(2, C_TO_VALUE2(min), C_TO_VALUE2(max));
}

static VALUE FUNCTION(rb_gsl_matrix,max_index)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t imax, jmax;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,max_index)(m, &imax, &jmax);
  return rb_ary_new3(2, INT2FIX(imax), INT2FIX(jmax));
}

static VALUE FUNCTION(rb_gsl_matrix,min_index)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t imin, jmin;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,min_index)(m, &imin, &jmin);
  return rb_ary_new3(2, INT2FIX(imin), INT2FIX(jmin));
}

static VALUE FUNCTION(rb_gsl_matrix,minmax_index)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t imin, jmin, imax, jmax;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,minmax_index)(m, &imin, &jmin, &imax, &jmax);
  return rb_ary_new3(2, rb_ary_new3(2, INT2FIX(imin), INT2FIX(jmin)),
         rb_ary_new3(2, INT2FIX(imax), INT2FIX(jmax)));
}

static VALUE FUNCTION(rb_gsl_matrix,fwrite)(VALUE obj, VALUE io)
{
  GSL_TYPE(gsl_matrix) *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), h);
  f = rb_gsl_open_writefile(io, &flag);
  status = FUNCTION(gsl_matrix,fwrite)(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_gsl_matrix,fread)(VALUE obj, VALUE io)
{
  GSL_TYPE(gsl_matrix) *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), h);
  f = rb_gsl_open_readfile(io, &flag);
  status = FUNCTION(gsl_matrix,fread)(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_gsl_matrix,trace)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t i;
  BASE trace = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  for (i = 0; i < m->size1; i++) {
    trace += FUNCTION(gsl_matrix,get)(m, i, i);
  }
  return C_TO_VALUE2(trace);
}


static VALUE FUNCTION(rb_gsl_matrix,uplus)(VALUE obj)
{
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,uminus)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mnew = NULL;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      FUNCTION(gsl_matrix,set)(mnew, i, j, -FUNCTION(gsl_matrix,get)(m, i, j));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
}

VALUE FUNCTION(rb_gsl_matrix,power)(VALUE obj, VALUE bb)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mtmp = NULL, *mnew = NULL;
  size_t i, b;
  CHECK_FIXNUM(bb);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  b = FIX2INT(bb); 
  mtmp = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  FUNCTION(gsl_matrix,memcpy)(mnew, m);
  for (i = 1; i < b; i++) {
    FUNCTION(gsl_matrix,memcpy)(mtmp, mnew);
#ifdef BASE_DOUBLE
    gsl_linalg_matmult(mtmp, m, mnew);
#else
    gsl_linalg_matmult_int(mtmp, m, mnew);
#endif
  }
  FUNCTION(gsl_matrix,free)(mtmp);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,submatrix)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  QUALIFIED_VIEW(gsl_matrix,view) *mv = NULL;
  QUALIFIED_VIEW(gsl_vector,view) *vv = NULL;
  size_t i, j, n1, n2;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  parse_submatrix_args(argc, argv, m->size1, m->size2, &i, &j, &n1, &n2);
  if(n1 == 0) {
    vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
    *vv = FUNCTION(gsl_matrix,subrow)(m, i, j, n2);
    return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,view), 0, free, vv);
  }
  else if(n2 == 0) {
    vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
    *vv = FUNCTION(gsl_matrix,subcolumn)(m, j, i, n1);
    return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,col_view), 0, free, vv);
  } else {
    mv = ALLOC(QUALIFIED_VIEW(gsl_matrix,view));
    *mv = FUNCTION(gsl_matrix,submatrix)(m, i, j, n1, n2);
    return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_matrix,view), 0, free, mv);
  }
}

static VALUE FUNCTION(rb_gsl_matrix,return_vector_view)(VALUE obj, VALUE index,
                QUALIFIED_VIEW(gsl_vector,view) (*f)(GSL_TYPE(gsl_matrix)*,
                   size_t))
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  QUALIFIED_VIEW(gsl_vector,view) *vv = NULL;
  CHECK_FIXNUM(index);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
  *vv = (*f)(m, FIX2INT(index));
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,view), 0, free, vv);
}

static VALUE FUNCTION(rb_gsl_matrix,row)(VALUE obj, VALUE i)
{
  return FUNCTION(rb_gsl_matrix,return_vector_view)(obj, i, FUNCTION(gsl_matrix,row));
}

static VALUE FUNCTION(rb_gsl_matrix,column)(VALUE obj, VALUE j)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  QUALIFIED_VIEW(gsl_vector,view) *vv = NULL;
  CHECK_FIXNUM(j);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
  *vv = FUNCTION(gsl_matrix,column)(m, FIX2INT(j));
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,col_view), 0, free, vv);
}

#ifdef GSL_1_10_LATER
static VALUE FUNCTION(rb_gsl_matrix,subrow)(VALUE obj, VALUE i, VALUE offset,
  VALUE n)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  QUALIFIED_VIEW(gsl_vector,view) *vv = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
  *vv = FUNCTION(gsl_matrix,subrow)(m, FIX2INT(i), FIX2INT(offset), FIX2INT(n));
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,view), 0, free, vv);
}
static VALUE FUNCTION(rb_gsl_matrix,subcolumn)(VALUE obj, VALUE j, VALUE offset,
  VALUE n)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  QUALIFIED_VIEW(gsl_vector,view) *vv = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
  *vv = FUNCTION(gsl_matrix,subcolumn)(m, FIX2INT(j), FIX2INT(offset), FIX2INT(n));
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,col_view), 0, free, vv);
}
#endif

static VALUE FUNCTION(rb_gsl_matrix,diagonal)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  QUALIFIED_VIEW(gsl_vector,view) *vv = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
  *vv = FUNCTION(gsl_matrix,diagonal)(m);
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,view), 0, free, vv);
}

static VALUE FUNCTION(rb_gsl_matrix,subdiagonal)(VALUE obj, VALUE k)
{
  return rb_gsl_matrix_return_vector_view(obj, k, gsl_matrix_subdiagonal);
}

static VALUE FUNCTION(rb_gsl_matrix,superdiagonal)(VALUE obj, VALUE k)
{
  return rb_gsl_matrix_return_vector_view(obj, k, gsl_matrix_superdiagonal);
}

static VALUE FUNCTION(rb_gsl_matrix,vector_view)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  QUALIFIED_VIEW(gsl_vector,view) *vv;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
  vv->vector.size = m->size1*m->size2;
  vv->vector.owner = 0;
  vv->vector.stride = 1;
  vv->vector.data = m->data;
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,view), 0, free, vv);
}


static VALUE FUNCTION(rb_gsl_matrix,each_row)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  QUALIFIED_VIEW(gsl_vector,view) *vv;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  for (i = 0; i < m->size1; i++) {
    vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
    *vv = FUNCTION(gsl_matrix,row)(m, i);
#ifdef BASE_DOUBLE
    rb_yield(Data_Wrap_Struct(cgsl_vector_view, 0, free, vv));
#else
    rb_yield(Data_Wrap_Struct(cgsl_vector_int_view, 0, free, vv));
#endif
  }
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,each_col)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  QUALIFIED_VIEW(gsl_vector,view) *vv;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  for (i = 0; i < m->size2; i++) {
    vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
    *vv = FUNCTION(gsl_matrix,column)(m, i);
#ifdef BASE_DOUBLE
    rb_yield(Data_Wrap_Struct(cgsl_vector_col_view, 0, free, vv));
#else
    rb_yield(Data_Wrap_Struct(cgsl_vector_int_col_view, 0, free, vv));
#endif
  }
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,scale_bang)(VALUE obj, VALUE x)
{
  GSL_TYPE(gsl_matrix) *m;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,scale)(m, NUMCONV(x));
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,scale)(VALUE obj, VALUE b)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mnew;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(make_matrix,clone)(m);
  FUNCTION(gsl_matrix,scale)(mnew, NUMCONV(b));
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,add_constant_bang)(VALUE obj, VALUE x)
{
  GSL_TYPE(gsl_matrix) *m;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(gsl_matrix,add_constant)(m, NUMCONV(x));
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,add_constant)(VALUE obj, VALUE b)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mnew;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(make_matrix,clone)(m);
  FUNCTION(gsl_matrix,add_constant)(mnew, NUMCONV(b));
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static int FUNCTION(mygsl_matrix,equal)(GSL_TYPE(gsl_matrix) *a, GSL_TYPE(gsl_matrix) *b, double eps)
{
  size_t i, j;
  BASE x, y;
  if (a->size1 != b->size1 || a->size2 != b->size2) return 0;
  for (i = 0; i < a->size1; i++) {
    for (j = 0; j < a->size2; j++) {
      x = FUNCTION(gsl_matrix,get)(a, i, j);
      y = FUNCTION(gsl_matrix,get)(b, i, j);
      if (fabs(x-y) > eps) return 0;
    }
  }
  return 1;
}
 
#ifdef HAVE_TENSOR_TENSOR_H
EXTERN VALUE cgsl_tensor, cgsl_tensor_int;
VALUE rb_gsl_tensor_equal(int argc, VALUE *argv, VALUE obj);
VALUE rb_gsl_tensor_int_equal(int argc, VALUE *argv, VALUE obj);
#ifdef BASE_DOUBLE
#define TEN_P(x) TENSOR_P(x)
#else
#define TEN_P(x) TENSOR_INT_P(x)
#endif
#endif

static VALUE FUNCTION(rb_gsl_matrix,equal)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *a, *b;
  double eps = 1e-10;
  VALUE bb;
  switch (argc) {
  case 2:
    bb = argv[0];
    eps = NUM2DBL(argv[1]);
    break;
  case 1:
    bb = argv[0];
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
#ifdef HAVE_TENSOR_TENSOR_H
  if (TEN_P(bb)) {
    return FUNCTION(rb_gsl_tensor,equal)(argc, argv, obj);
  }
#endif
  CHECK_MAT(bb);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), a);
  Data_Get_Struct(bb, GSL_TYPE(gsl_matrix), b);
  if (FUNCTION(mygsl_matrix,equal)(a, b, eps) == 1) return Qtrue;
  else return Qfalse;
}

#ifdef HAVE_TENSOR_TENSOR_H
#ifdef TEN_P
#undef TEN_P
#endif
#endif

static VALUE FUNCTION(rb_gsl_matrix,equal_singleton)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *a, *b;
  VALUE aa, bb;
  double eps = 1e-10;
  BASE x, y;
  size_t i, j;
  switch (argc) {
  case 3:
    aa = argv[0];
    bb = argv[1];
    eps = NUM2DBL(argv[2]);
    break;
  case 2:
    aa = argv[0];
    bb = argv[1];
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  CHECK_MAT(aa);  CHECK_MAT(bb);
  Data_Get_Struct(aa, GSL_TYPE(gsl_matrix), a);
  Data_Get_Struct(bb, GSL_TYPE(gsl_matrix), b);
  if (a->size1 != b->size1 || a->size2 != b->size2) return Qfalse;
  for (i = 0; i < a->size1; i++) {
    for (j = 0; j < a->size2; j++) {
      x = FUNCTION(gsl_matrix,get)(a, i, j);
      y = FUNCTION(gsl_matrix,get)(b, i, j);
      if (fabs(x-y) > eps) return Qfalse;
    }
  }
  return Qtrue;
}

#ifdef HAVE_TENSOR_TENSOR_H
#include "include/rb_gsl_tensor.h"
static VALUE FUNCTION(rb_gsl_matrix,to_tensor)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  GSL_TYPE(rbgsl_tensor) *t;
  unsigned int rank;
  size_t dim;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  if (m->size1 != m->size2) rb_raise(rb_eRuntimeError, "matrix must have equal dimensions");
  rank = 2;
  dim = m->size1;
  t = FUNCTION(rbgsl_tensor,alloc)(rank, dim);
  memcpy(t->tensor->data, m->data, sizeof(BASE)*t->tensor->size);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), t);
}
#endif

static VALUE FUNCTION(rb_gsl_matrix,collect)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mnew;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      FUNCTION(gsl_matrix,set)(mnew, i, j, NUMCONV(rb_yield(C_TO_VALUE(FUNCTION(gsl_matrix,get)(m, i, j)))));
    }
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,collect_bang)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      FUNCTION(gsl_matrix,set)(m, i, j, NUMCONV(rb_yield(C_TO_VALUE(FUNCTION(gsl_matrix,get)(m, i, j)))));
    }
  }
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,upper)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mnew;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(make_matrix,clone)(m);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < i; j++) {
      FUNCTION(gsl_matrix,set)(mnew, i, j, 0);
    }
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,lower)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m = NULL, *mnew;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(make_matrix,clone)(m);
  for (i = 0; i < m->size1; i++) {
    for (j = i+1; j < m->size2; j++) {
      FUNCTION(gsl_matrix,set)(mnew, i, j, 0);
    }
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

/* singleton method, creates Pascal matrix of dimension n */
static VALUE FUNCTION(rb_gsl_matrix,pascal1)(VALUE obj, VALUE n)
{
  GSL_TYPE(gsl_matrix) *m;
  BASE x;
  size_t i, j, dim;
  CHECK_FIXNUM(n);
  dim = (size_t) FIX2INT(n);
  m = FUNCTION(gsl_matrix,alloc)(dim, dim);
  for (j = 0; j < dim; j++) FUNCTION(gsl_matrix,set)(m, 0, j, (BASE) 1);
  for (i = 1; i < dim; i++) {
    FUNCTION(gsl_matrix,set)(m, i, 0, (BASE) 1);
    for (j = 1; j < dim; j++) {
      x=FUNCTION(gsl_matrix,get)(m,i-1,j)+FUNCTION(gsl_matrix,get)(m,i,j-1);
      FUNCTION(gsl_matrix,set)(m, i, j, x);
    }
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

#ifdef BASE_DOUBLE
static VALUE FUNCTION(rb_gsl_matrix,hilbert)(VALUE obj, VALUE n)
{
  GSL_TYPE(gsl_matrix) *m;
  double x;
  size_t i, j, dim;
  CHECK_FIXNUM(n);
  dim = (size_t) FIX2INT(n);
  m = FUNCTION(gsl_matrix,alloc)(dim, dim);
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      x = 1.0/(i + j + 1);
      FUNCTION(gsl_matrix,set)(m, i, j, (BASE) x);
    }
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

double mygsl_binomial_coef(unsigned int n, unsigned int k);
static VALUE FUNCTION(rb_gsl_matrix,invhilbert)(VALUE obj, VALUE n)
{
  GSL_TYPE(gsl_matrix) *m;
  double x, y;
  size_t i, j, dim;
  CHECK_FIXNUM(n);
  dim = (size_t) FIX2INT(n);
  m = FUNCTION(gsl_matrix,alloc)(dim, dim);
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if ((i+j)%2 == 0) x = 1;
      else x = -1;
      x *= (i + j + 1);
      x *= mygsl_binomial_coef(dim + i, dim - j - 1);
      x *= mygsl_binomial_coef(dim + j, dim - i - 1);
      y = mygsl_binomial_coef(i + j, i);
      x *= y*y;
      FUNCTION(gsl_matrix,set)(m, i, j, (BASE) x);
    }
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}
#endif

static void FUNCTION(mygsl_matrix,vandermonde)(GSL_TYPE(gsl_matrix) *m, GSL_TYPE(gsl_vector) *v)
{
  size_t i, j;

  for (i = 0; i < v->size; i++) {
    for (j = 0; j < v->size; j++) {
      FUNCTION(gsl_matrix,set)(m, i, j, (BASE) gsl_pow_int(FUNCTION(gsl_vector,get)(v, i), v->size-j-1));
    }
  }
}

/* singleton */
static VALUE FUNCTION(rb_gsl_matrix,vandermonde)(VALUE obj, VALUE vv)
{
  GSL_TYPE(gsl_vector) *v;
  GSL_TYPE(gsl_matrix) *m;
  int flag = 0;
  if (TYPE(vv) == T_ARRAY) {
    v = FUNCTION(make_cvector,from_rarray)(vv);
    flag = 1;
  } else if (VEC_P(vv)) {
    Data_Get_Struct(vv, GSL_TYPE(gsl_vector), v);
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Array or Vector expected)",
       rb_class2name(CLASS_OF(vv)));
  }
  m = FUNCTION(gsl_matrix,alloc)(v->size, v->size);
  FUNCTION(mygsl_matrix,vandermonde)(m, v);
  if (flag == 1) FUNCTION(gsl_vector,free)(v);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

static void FUNCTION(mygsl_matrix,toeplitz)(GSL_TYPE(gsl_matrix) *m, GSL_TYPE(gsl_vector) *v)
{
  size_t i, j;
  for (i = 0; i < v->size; i++) {
    for (j = 0; j < v->size; j++) {
      if (j >= i) 
  FUNCTION(gsl_matrix,set)(m, i, j, FUNCTION(gsl_vector,get)(v, j-i));
      else
  FUNCTION(gsl_matrix,set)(m, i, j, FUNCTION(gsl_vector,get)(v, i-j));
    }
  }
}

static VALUE FUNCTION(rb_gsl_matrix,toeplitz)(VALUE obj, VALUE vv)
{
  GSL_TYPE(gsl_vector) *v;
  GSL_TYPE(gsl_matrix) *m;
  int flag = 0;
  if (TYPE(vv) == T_ARRAY) {
    v = FUNCTION(make_cvector,from_rarray)(vv);
    flag = 1;
  } else if (VEC_P(vv)) {
    Data_Get_Struct(vv, GSL_TYPE(gsl_vector), v);
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Array or Vector expected)",
       rb_class2name(CLASS_OF(vv)));
  }
  m = FUNCTION(gsl_matrix,alloc)(v->size, v->size);
  FUNCTION(mygsl_matrix,toeplitz)(m, v);
  if (flag == 1) FUNCTION(gsl_vector,free)(v);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

void FUNCTION(mygsl_vector,to_m_circulant)(GSL_TYPE(gsl_matrix) *m, GSL_TYPE(gsl_vector) *v);
static VALUE FUNCTION(rb_gsl_matrix,circulant)(VALUE obj, VALUE vv)
{
  GSL_TYPE(gsl_vector) *v;
  GSL_TYPE(gsl_matrix) *m;
  int flag = 0;
  if (TYPE(vv) == T_ARRAY) {
    v = FUNCTION(make_cvector,from_rarray)(vv);
    flag = 1;
  } else if (VEC_P(vv)) {
    Data_Get_Struct(vv, GSL_TYPE(gsl_vector), v);
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Array or Vector expected)",
       rb_class2name(CLASS_OF(vv)));
  }
  m = FUNCTION(gsl_matrix,alloc)(v->size, v->size);
  FUNCTION(mygsl_vector,to_m_circulant)(m, v);
  if (flag == 1) FUNCTION(gsl_vector,free)(v);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

static void FUNCTION(mygsl_matrix,indgen)(GSL_TYPE(gsl_matrix) *m,
            BASE start, BASE step)
{
  size_t i, j;
  BASE n;
  n = start;
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      FUNCTION(gsl_matrix,set)(m, i, j, n);
      n += step;
    }
  }
}

static VALUE FUNCTION(rb_gsl_matrix,indgen_singleton)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  size_t n1, n2;
  BASE start = 0, step = 1;
  switch (argc) {
  case 4:
    step = NUMCONV2(argv[3]);
    /* no break */
  case 3:
    start = NUMCONV2(argv[2]);
    /* no break */
  case 2:
    n1 = NUM2INT(argv[0]);
    n2 = NUM2INT(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2-4)", argc);
    break;
  }
  m = FUNCTION(gsl_matrix,alloc)(n1, n2);
  FUNCTION(mygsl_matrix,indgen)(m, start, step);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

static VALUE FUNCTION(rb_gsl_matrix,indgen)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m, *mnew;
  BASE start = 0, step = 1;
  switch (argc) {
  case 2:
    step = NUMCONV2(argv[1]);
    /* no break */
  case 1:
    start = NUMCONV2(argv[0]);
    break;
  case 0:
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0-2)", argc);
    break;
  }
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  FUNCTION(mygsl_matrix,indgen)(mnew, start, step);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,indgen_bang)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  BASE start = 0, step = 1;
  switch (argc) {
  case 2:
    step = NUMCONV2(argv[1]);
    /* no break */
  case 1:
    start = NUMCONV2(argv[0]);
    break;
  case 0:
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0-2)", argc);
    break;
  }
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  FUNCTION(mygsl_matrix,indgen)(m, start, step);
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,to_a)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  VALUE ma, ra;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  ma = rb_ary_new2(m->size1);
  for(i=0; i < m->size1; i++) {
    ra = rb_ary_new2(m->size2);
    rb_ary_store(ma, i, ra);
    for(j=0; j < m->size2; j++) {
      rb_ary_store(ra, j, C_TO_VALUE(FUNCTION(gsl_matrix,get)(m, i, j)));
    }
  }
  return ma;
}

static VALUE FUNCTION(rb_gsl_matrix,to_v)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  GSL_TYPE(gsl_vector) *v;
  size_t i, j, k;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  v = FUNCTION(gsl_vector,alloc)(m->size1*m->size2);
  //  memcpy(v->data, m->data, sizeof(BASE)*v->size);
  for (i = 0, k = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++, k++) {
      FUNCTION(gsl_vector,set)(v, k, FUNCTION(gsl_matrix,get)(m, i, j));
    }
  }

  if(m->size1 > 1 && m->size2 == 1) {
    return Data_Wrap_Struct(CONCAT2(GSL_TYPE(cgsl_vector),col), 0, FUNCTION(gsl_vector,free), v);
  } else {
    return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), v);
  }
}

static VALUE FUNCTION(rb_gsl_matrix,to_vview)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  QUALIFIED_VIEW(gsl_vector,view) *v;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  v = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
  v->vector.size = m->size1*m->size2;
  v->vector.stride = 1;
  v->vector.owner = 0;
  v->vector.data = m->data;
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,view), 0, FUNCTION(gsl_vector,free), v);
}

static VALUE FUNCTION(rb_gsl_matrix,norm)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  size_t i, n;
  BASE x = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  n = m->size1*m->size2;
  for (i = 0; i < n; i++) x += m->data[i]*m->data[i];
  return rb_float_new(sqrt(x));
}

static int FUNCTION(mygsl_matrix,reverse_columns)(GSL_TYPE(gsl_matrix) *dst,
               GSL_TYPE(gsl_matrix) *src)
{
  size_t j;
  QUALIFIED_VIEW(gsl_vector,view) col;
  if (dst->size1 != src->size1 || dst->size2 != src->size2)
    rb_raise(rb_eRuntimeError, "matrix sizes are different.");
  for (j = 0; j < src->size2; j++) {
    col = FUNCTION(gsl_matrix,column)(src, j);
    FUNCTION(gsl_matrix,set_col)(dst, dst->size2-1-j, &col.vector);
  }
  return 0;
}

static int FUNCTION(mygsl_matrix,reverse_rows)(GSL_TYPE(gsl_matrix) *dst,
               GSL_TYPE(gsl_matrix) *src)
{
  size_t i;
  QUALIFIED_VIEW(gsl_vector,view) row;
  if (dst->size1 != src->size1 || dst->size2 != src->size2)
    rb_raise(rb_eRuntimeError, "matrix sizes are different.");
  for (i = 0; i < src->size1; i++) {
    row = FUNCTION(gsl_matrix,row)(src, i);
    FUNCTION(gsl_matrix,set_row)(dst, dst->size1-1-i, &row.vector);
  }
  return 0;
}

static VALUE FUNCTION(rb_gsl_matrix,reverse_columns)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m, *mnew;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  FUNCTION(mygsl_matrix,reverse_columns)(mnew, m);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,reverse_columns_bang)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m, *mnew;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  FUNCTION(mygsl_matrix,reverse_columns)(mnew, m);
  FUNCTION(gsl_matrix,memcpy)(m, mnew);
  FUNCTION(gsl_matrix,free)(mnew);
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,reverse_rows)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m, *mnew;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  FUNCTION(mygsl_matrix,reverse_rows)(mnew, m);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,reverse_rows_bang)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m, *mnew;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  FUNCTION(mygsl_matrix,reverse_rows)(mnew, m);
  FUNCTION(gsl_matrix,memcpy)(m, mnew);
  FUNCTION(gsl_matrix,free)(mnew);
  return obj;
}

static VALUE FUNCTION(rb_gsl_matrix,block)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_block), 0, NULL, m->block);
}

static VALUE FUNCTION(rb_gsl_matrix,info)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  char buf[256];
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  sprintf(buf, "Class:      %s\n", rb_class2name(CLASS_OF(obj)));
  sprintf(buf, "%sSuperClass: %s\n", buf, rb_class2name(RCLASS_SUPER(CLASS_OF(obj))));
  sprintf(buf, "%sDimension:  %dx%d\n", buf, (int) m->size1, (int) m->size2);
  sprintf(buf, "%sSize:       %d\n", buf, (int) (m->size1*m->size2));
  return rb_str_new2(buf);
}

static VALUE FUNCTION(rb_gsl_matrix,any)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  QUALIFIED_VIEW(gsl_vector,view) vv;
  GSL_TYPE(gsl_vector) *v;
  gsl_vector_int *vnew;
  size_t j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  vnew = gsl_vector_int_alloc(m->size2);
  for (j = 0; j < m->size2; j++) {
    vv = FUNCTION(gsl_matrix,column)(m, j);
    v = &vv.vector;
    if (FUNCTION(gsl_vector,isnull(v))) gsl_vector_int_set(vnew, j, 0);
    else gsl_vector_int_set(vnew, j, 1);
  }
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, vnew);
}

static VALUE FUNCTION(rb_gsl_matrix,all)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  QUALIFIED_VIEW(gsl_vector,view) vv;
  GSL_TYPE(gsl_vector) *v;
  gsl_vector_int *vnew;
  size_t i, j;
  int flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  vnew = gsl_vector_int_alloc(m->size2);
  for (j = 0; j < m->size2; j++) {
    vv = FUNCTION(gsl_matrix,column)(m, j);
    v = &vv.vector;
    /*    if (FUNCTION(gsl_vector,isnull(v))) gsl_vector_int_set(vnew, j, 0);
    else gsl_vector_int_set(vnew, j, 1);*/
    for (i = 0; i < v->size; i++) {
      if (FUNCTION(gsl_vector,get)(v, i) == (BASE) 0) {
  gsl_vector_int_set(vnew, j, 0);
  flag = 0;
  break;
      } else {
  flag = 1;
      }
    }
    if (flag == 1) gsl_vector_int_set(vnew, j, 1);
  }
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, vnew);
}

static VALUE FUNCTION(rb_gsl_matrix,rot90)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m, *mnew, *mtmp;
  int p;
  switch (argc) {
  case 0:
    p = 1;
    break;
  case 1:
    p = FIX2INT(argv[0])%4;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    break;
  }
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  switch (p) {
  case 0:
    mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
    //    FUNCTION(gsl_matrix,memcpy(mnew, m));
    FUNCTION(gsl_matrix,memcpy)(mnew, m);
    break;
  case 1:
  case -3:
    mtmp = FUNCTION(gsl_matrix,alloc)(m->size2, m->size1);
    FUNCTION(gsl_matrix,transpose_memcpy(mtmp, m));
    mnew = FUNCTION(gsl_matrix,alloc)(m->size2, m->size1);
    FUNCTION(mygsl_matrix,reverse_rows)(mnew, mtmp);
    FUNCTION(gsl_matrix,free)(mtmp);
    break;
  case 2:
  case -2:
    mtmp = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
    FUNCTION(mygsl_matrix,reverse_rows)(mtmp,m);
    mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
    FUNCTION(mygsl_matrix,reverse_columns)(mnew, mtmp);
    FUNCTION(gsl_matrix,free)(mtmp);
    break;
  case 3:
  case -1:
    mtmp = FUNCTION(gsl_matrix,alloc)(m->size2, m->size1);
    FUNCTION(gsl_matrix,transpose_memcpy(mtmp, m));
    mnew = FUNCTION(gsl_matrix,alloc)(m->size2, m->size1);
    FUNCTION(mygsl_matrix,reverse_columns)(mnew, mtmp);
    FUNCTION(gsl_matrix,free)(mtmp);
    break;
  default:
    return Qnil;
    break;
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

void FUNCTION(mygsl_vector,diff)(GSL_TYPE(gsl_vector) *vdst,
         GSL_TYPE(gsl_vector) *vsrc, size_t n);

static VALUE FUNCTION(rb_gsl_matrix,diff)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m, *mnew;
  QUALIFIED_VIEW(gsl_vector,view) v1, v2;
  size_t n, j;
  switch (argc) {
  case 0:
    n = 1;
    break;
  case 1:
    n = FIX2INT(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    break;
  }
  if (n <= 0) return obj;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1-n, m->size2);
  if (m->size1 <= n) return obj;
  for (j = 0; j < m->size2; j++) {
    v1 = FUNCTION(gsl_matrix,column)(m, j);
    v2 = FUNCTION(gsl_matrix,column)(mnew, j);
    FUNCTION(mygsl_vector,diff)(&v2.vector, &v1.vector, n);
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,test)(VALUE obj, int (*f)(const double))
{
  GSL_TYPE(gsl_matrix) *m;
  gsl_matrix_int *mi;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mi = gsl_matrix_int_alloc(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      gsl_matrix_int_set(mi, i, j, (*f)(FUNCTION(gsl_matrix,get)(m, i, j)));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix_int, 0, gsl_matrix_int_free, mi);
}

static VALUE FUNCTION(rb_gsl_matrix,isnan)(VALUE obj)
{
  return FUNCTION(rb_gsl_matrix,test)(obj, gsl_isnan);
}

static VALUE FUNCTION(rb_gsl_matrix,isinf)(VALUE obj)
{
  return FUNCTION(rb_gsl_matrix,test)(obj, gsl_isinf);
}

static VALUE FUNCTION(rb_gsl_matrix,finite)(VALUE obj)
{
  return FUNCTION(rb_gsl_matrix,test)(obj, gsl_finite);
}

static VALUE FUNCTION(rb_gsl_matrix,sgn)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m, *mnew;
  BASE x;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      x = FUNCTION(gsl_matrix,get)(m, i, j);
      FUNCTION(gsl_matrix,set)(mnew, i, j, (BASE)(x>0 ? 1 : (x<0 ? -1 : 0)));
    }
  }

  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,abs)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m, *mnew;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      FUNCTION(gsl_matrix,set)(mnew, i, j, (BASE) fabs(FUNCTION(gsl_matrix,get)(m, i, j)));
    }
  }

  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,horzcat)(VALUE obj, VALUE mm2)
{
  GSL_TYPE(gsl_matrix) *m, *m2, *mnew;
  QUALIFIED_VIEW(gsl_vector,view) v;
  size_t j, k;
  CHECK_MAT(mm2);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  Data_Get_Struct(mm2, GSL_TYPE(gsl_matrix), m2);
  if (m->size1 != m2->size1) 
    rb_raise(rb_eRuntimeError, "Different number of rows (%d and %d).",
       (int) m->size1, (int) m2->size1);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2 + m2->size2);
  for (j = 0, k = 0; j < m->size2; j++, k++) {
    v = FUNCTION(gsl_matrix,column)(m, j);
    FUNCTION(gsl_matrix,set_col)(mnew, k, &v.vector);
  }
  for (j = 0; j < m2->size2; j++, k++) {
    v = FUNCTION(gsl_matrix,column)(m2, j);
    FUNCTION(gsl_matrix,set_col)(mnew, k, &v.vector);
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,horzcat_singleton)(VALUE klass, VALUE mm, VALUE mm2)
{
  CHECK_MAT(mm);
  return FUNCTION(rb_gsl_matrix,horzcat)(mm, mm2);
}


static VALUE FUNCTION(rb_gsl_matrix,vertcat)(VALUE obj, VALUE mm2)
{
  GSL_TYPE(gsl_matrix) *m, *m2, *mnew;
  QUALIFIED_VIEW(gsl_vector,view) v;
  size_t i, k;
  CHECK_MAT(mm2);
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  Data_Get_Struct(mm2, GSL_TYPE(gsl_matrix), m2);
  if (m->size2 != m2->size2) 
    rb_raise(rb_eRuntimeError, "Different number of columns (%d and %d).",
       (int) m->size2, (int) m2->size2);
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1 + m2->size1, m->size2);
  for (i = 0, k = 0; i < m->size1; i++, k++) {
    v = FUNCTION(gsl_matrix,row)(m, i);
    FUNCTION(gsl_matrix,set_row)(mnew, k, &v.vector);
  }
  for (i = 0; i < m2->size1; i++, k++) {
    v = FUNCTION(gsl_matrix,row)(m2, i);
    FUNCTION(gsl_matrix,set_row)(mnew, k, &v.vector);
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,vertcat_singleton)(VALUE klass, VALUE mm, VALUE mm2)
{
  CHECK_MAT(mm);
  return FUNCTION(rb_gsl_matrix,vertcat)(mm, mm2);
}

#ifdef GSL_1_9_LATER
static VALUE FUNCTION(rb_gsl_matrix,property)(VALUE obj,
  int (*f)(const GSL_TYPE(gsl_matrix)*)) {
  GSL_TYPE(gsl_matrix) *m;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  return INT2FIX((*f)(m));  
}

static VALUE FUNCTION(rb_gsl_matrix,property2)(VALUE obj,
  int (*f)(const GSL_TYPE(gsl_matrix) *)) {
  GSL_TYPE(gsl_matrix) *m;  
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  if ((*f)(m)) return Qtrue;
  else return Qfalse;
}
static VALUE FUNCTION(rb_gsl_matrix,ispos)(VALUE obj)
{
  return FUNCTION(rb_gsl_matrix,property)(obj, FUNCTION(gsl_matrix,ispos));
}
static VALUE FUNCTION(rb_gsl_matrix,ispos2)(VALUE obj)
{
  return FUNCTION(rb_gsl_matrix,property2)(obj, FUNCTION(gsl_matrix,ispos));  
}
static VALUE FUNCTION(rb_gsl_matrix,isneg)(VALUE obj)
{
  return FUNCTION(rb_gsl_matrix,property)(obj, FUNCTION(gsl_matrix,isneg));  
}
static VALUE FUNCTION(rb_gsl_matrix,isneg2)(VALUE obj)
{
  return FUNCTION(rb_gsl_matrix,property2)(obj, FUNCTION(gsl_matrix,isneg));  
}
#endif

#ifdef GSL_1_10_LATER
static VALUE FUNCTION(rb_gsl_matrix,isnonneg)(VALUE obj)
{
  return FUNCTION(rb_gsl_matrix,property)(obj, FUNCTION(gsl_matrix,isnonneg));  
}
static VALUE FUNCTION(rb_gsl_matrix,isnonneg2)(VALUE obj)
{
  return FUNCTION(rb_gsl_matrix,property2)(obj, FUNCTION(gsl_matrix,isnonneg));  
}
#endif

static VALUE FUNCTION(rb_gsl_matrix,symmetrize)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m, *mnew;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  if (m->size1 != m->size2)
    rb_raise(rb_eRuntimeError, "symmetrize: not a square matrix.\n");
  mnew = FUNCTION(gsl_matrix,alloc)(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = i; j < m->size2; j++) {
      FUNCTION(gsl_matrix,set)(mnew, i, j, FUNCTION(gsl_matrix,get)(m, i, j));
    }
    for (j = 0; j < i; j++) {
      FUNCTION(gsl_matrix,set)(mnew, i, j, FUNCTION(gsl_matrix,get)(m, j, i));
    }
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), mnew);
}

static VALUE FUNCTION(rb_gsl_matrix,symmetrize_bang)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  size_t i, j;
  Data_Get_Struct(obj, GSL_TYPE(gsl_matrix), m);
  if (m->size1 != m->size2)
    rb_raise(rb_eRuntimeError, "symmetrize: not a square matrix.\n");
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < i; j++) {
      FUNCTION(gsl_matrix,set)(m, i, j, FUNCTION(gsl_matrix,get)(m, j, i));
    }
  }
  return obj;
}

void FUNCTION(Init_gsl_matrix,init)(VALUE module)
{
  /*  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "new", FUNCTION(rb_gsl_matrix,alloc), -1);*/
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "[]", FUNCTION(rb_gsl_matrix,alloc), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "alloc", FUNCTION(rb_gsl_matrix,alloc), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "calloc", FUNCTION(rb_gsl_matrix,calloc), 2);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "eye", FUNCTION(rb_gsl_matrix,eye), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "ones", FUNCTION(rb_gsl_matrix,ones), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "zeros", FUNCTION(rb_gsl_matrix,zeros), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "diagonal", 
           FUNCTION(rb_gsl_matrix,diagonal_singleton), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "diag", 
           FUNCTION(rb_gsl_matrix,diagonal_singleton), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "identity", 
           FUNCTION(rb_gsl_matrix,identity), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "scalar", 
           FUNCTION(rb_gsl_matrix,identity), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "unit", 
           FUNCTION(rb_gsl_matrix,identity), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "I", 
           FUNCTION(rb_gsl_matrix,identity), 1);

  /*****/
  rb_define_method(GSL_TYPE(cgsl_matrix), "size1", 
       FUNCTION(rb_gsl_matrix,size1), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "size2", 
       FUNCTION(rb_gsl_matrix,size2), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "shape", 
       FUNCTION(rb_gsl_matrix,shape), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "size", "shape");

  rb_define_method(GSL_TYPE(cgsl_matrix), "get", FUNCTION(rb_gsl_matrix,get), -1);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "[]", "get");
  rb_define_method(GSL_TYPE(cgsl_matrix), "set", FUNCTION(rb_gsl_matrix,set), -1);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "[]=", "set");

  rb_define_method(GSL_TYPE(cgsl_matrix), "set_all", 
       FUNCTION(rb_gsl_matrix,set_all), 1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "set_zero",  
       FUNCTION(rb_gsl_matrix,set_zero), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "set_identity",  
       FUNCTION(rb_gsl_matrix,set_identity), 0);

  rb_define_method(GSL_TYPE(cgsl_matrix), "print", 
       FUNCTION(rb_gsl_matrix,print), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "inspect",  
       FUNCTION(rb_gsl_matrix,inspect), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "to_s", 
       FUNCTION(rb_gsl_matrix,to_s), 0);
  
  rb_define_method(GSL_TYPE(cgsl_matrix), "set_diagonal", 
       FUNCTION(rb_gsl_matrix,set_diagonal), 1);
  
  rb_define_method(GSL_TYPE(cgsl_matrix), "get_row", 
       FUNCTION(rb_gsl_matrix,get_row), 1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "get_column", 
       FUNCTION(rb_gsl_matrix,get_col), 1);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "get_col", "get_column");
  rb_define_method(GSL_TYPE(cgsl_matrix), "set_column", 
       FUNCTION(rb_gsl_matrix,set_col), 2);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "set_col", "set_column");
  rb_define_method(GSL_TYPE(cgsl_matrix), "set_row", 
       FUNCTION(rb_gsl_matrix,set_row), 2);
  
  rb_define_method(GSL_TYPE(cgsl_matrix), "clone", 
       FUNCTION(rb_gsl_matrix,clone), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "duplicate", "clone");
  rb_define_alias(GSL_TYPE(cgsl_matrix), "dup", "clone");
  rb_define_method(GSL_TYPE(cgsl_matrix), "isnull", 
       FUNCTION(rb_gsl_matrix,isnull), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "isnull?", 
       FUNCTION(rb_gsl_matrix,isnull2), 0);
  
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "memcpy", 
           FUNCTION(rb_gsl_matrix,memcpy), 2);
  
  rb_define_method(GSL_TYPE(cgsl_matrix), "swap_rows", 
       FUNCTION(rb_gsl_matrix,swap_rows), 2);
  rb_define_method(GSL_TYPE(cgsl_matrix), "swap_rows!",
       FUNCTION(rb_gsl_matrix,swap_rows_bang), 2);
  rb_define_method(GSL_TYPE(cgsl_matrix), "swap_columns",
       FUNCTION(rb_gsl_matrix,swap_columns), 2);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "swap_cols", "swap_columns");
  rb_define_method(GSL_TYPE(cgsl_matrix), "swap_columns!",
       FUNCTION(rb_gsl_matrix,swap_columns_bang), 2);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "swap_cols!", "swap_columns!");
  rb_define_method(GSL_TYPE(cgsl_matrix), "swap_rowcol", 
       FUNCTION(rb_gsl_matrix,swap_rowcol), 2);
  rb_define_method(GSL_TYPE(cgsl_matrix), "swap_rowcol!", 
       FUNCTION(rb_gsl_matrix,swap_rowcol_bang), 2);
  rb_define_method(GSL_TYPE(cgsl_matrix), "transpose_memcpy",
       FUNCTION(rb_gsl_matrix,transpose_memcpy), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "transpose", "transpose_memcpy");
  rb_define_alias(GSL_TYPE(cgsl_matrix), "trans", "transpose_memcpy");
  rb_define_method(GSL_TYPE(cgsl_matrix), "transpose!", 
       FUNCTION(rb_gsl_matrix,transpose_bang), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "trans!", "transpose!");
  rb_define_method(GSL_TYPE(cgsl_matrix), "reverse_columns", 
       FUNCTION(rb_gsl_matrix,reverse_columns), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "fliplr", "reverse_columns");
  rb_define_method(GSL_TYPE(cgsl_matrix), "reverse_columns!", 
       FUNCTION(rb_gsl_matrix,reverse_columns_bang), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "reverse_rows", 
       FUNCTION(rb_gsl_matrix,reverse_rows), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "flipud", "reverse_rows");
  rb_define_method(GSL_TYPE(cgsl_matrix), "reverse_rows!", 
       FUNCTION(rb_gsl_matrix,reverse_rows_bang), 0);
  /*****/
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "swap", 
           FUNCTION(rb_gsl_matrix,swap), 2);

  rb_define_method(GSL_TYPE(cgsl_matrix), "max", FUNCTION(rb_gsl_matrix,max), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "min", FUNCTION(rb_gsl_matrix,min), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "minmax",
       FUNCTION(rb_gsl_matrix,minmax), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "max_index",
       FUNCTION(rb_gsl_matrix,max_index), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "min_index", 
       FUNCTION(rb_gsl_matrix,min_index), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "minmax_index",
       FUNCTION(rb_gsl_matrix,minmax_index), 0);

  rb_define_method(GSL_TYPE(cgsl_matrix), "fwrite", 
       FUNCTION(rb_gsl_matrix,fwrite), 1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "fread", 
       FUNCTION(rb_gsl_matrix,fread), 1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "fprintf",
       FUNCTION(rb_gsl_matrix,fprintf), -1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "printf", 
       FUNCTION(rb_gsl_matrix,printf), -1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "fscanf",
       FUNCTION(rb_gsl_matrix,fscanf), 1);

  rb_define_method(GSL_TYPE(cgsl_matrix), "trace", 
       FUNCTION(rb_gsl_matrix,trace), 0);

  rb_define_method(GSL_TYPE(cgsl_matrix), "-@",  
       FUNCTION(rb_gsl_matrix,uminus), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "+@",  
       FUNCTION(rb_gsl_matrix,uplus), 0);


/*****/
  rb_define_method(GSL_TYPE(cgsl_matrix), "submatrix", 
       FUNCTION(rb_gsl_matrix,submatrix), -1);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "view", "submatrix");

  rb_define_method(GSL_TYPE(cgsl_matrix), "row", FUNCTION(rb_gsl_matrix,row), 1);
  /*  rb_define_alias(GSL_TYPE(cgsl_matrix), "[]", "row");*/
  rb_define_method(GSL_TYPE(cgsl_matrix), "column", 
       FUNCTION(rb_gsl_matrix,column), 1);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "col", "column");

#ifdef GSL_1_10_LATER
 rb_define_method(GSL_TYPE(cgsl_matrix), "subrow", 
       FUNCTION(rb_gsl_matrix,subrow), 3);
 rb_define_method(GSL_TYPE(cgsl_matrix), "subcolumn", 
       FUNCTION(rb_gsl_matrix,subcolumn), 3);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "subcol", "subcolumn");       
#endif

  rb_define_method(GSL_TYPE(cgsl_matrix), "diagonal", 
       FUNCTION(rb_gsl_matrix,diagonal), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "diag", "diagonal");

  rb_define_method(GSL_TYPE(cgsl_matrix), "subdiagonal", 
       FUNCTION(rb_gsl_matrix,subdiagonal), 1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "superdiagonal", 
       FUNCTION(rb_gsl_matrix,superdiagonal), 1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "vector_view", 
       FUNCTION(rb_gsl_matrix,vector_view), 0);

  rb_define_method(GSL_TYPE(cgsl_matrix), "each_row", 
       FUNCTION(rb_gsl_matrix,each_row), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "each_col",
       FUNCTION(rb_gsl_matrix,each_col), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "each_column", "each_col");

  rb_define_method(GSL_TYPE(cgsl_matrix), "scale", 
       FUNCTION(rb_gsl_matrix,scale), 1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "scale!",
       FUNCTION(rb_gsl_matrix,scale_bang), 1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "add_constant",
       FUNCTION(rb_gsl_matrix,add_constant), 1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "add_constant!", 
       FUNCTION(rb_gsl_matrix,add_constant_bang), 1);

  rb_define_method(GSL_TYPE(cgsl_matrix), "equal?", 
       FUNCTION(rb_gsl_matrix,equal), -1);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "==", "equal?");
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "equal?", 
           FUNCTION(rb_gsl_matrix,equal_singleton), -1);

  rb_define_method(GSL_TYPE(cgsl_matrix), "power", 
       FUNCTION(rb_gsl_matrix,power), 1);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "**", "power");
  rb_define_alias(GSL_TYPE(cgsl_matrix), "^", "power");

  rb_define_method(GSL_TYPE(cgsl_matrix), "collect", 
       FUNCTION(rb_gsl_matrix,collect), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "collect!", 
       FUNCTION(rb_gsl_matrix,collect_bang), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "map", "collect");
  rb_define_alias(GSL_TYPE(cgsl_matrix), "map!", "collect!");
#ifdef HAVE_TENSOR_TENSOR_H
  rb_define_method(GSL_TYPE(cgsl_matrix), "to_tensor", 
       FUNCTION(rb_gsl_matrix,to_tensor), 0);
#endif

  /*****/

  rb_define_method(GSL_TYPE(cgsl_matrix), "to_a", 
       FUNCTION(rb_gsl_matrix,to_a), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "to_v", 
       FUNCTION(rb_gsl_matrix,to_v), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "to_vview", 
       FUNCTION(rb_gsl_matrix,to_vview), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "data", "to_vview");
  rb_define_method(GSL_TYPE(cgsl_matrix), "norm",
       FUNCTION(rb_gsl_matrix,norm), 0);

  rb_define_method(GSL_TYPE(cgsl_matrix), "upper", 
       FUNCTION(rb_gsl_matrix,upper), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "lower", 
       FUNCTION(rb_gsl_matrix,lower), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "block", 
       FUNCTION(rb_gsl_matrix,block), 0);

  /***** Special matrices *****/
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "pascal",
           FUNCTION(rb_gsl_matrix,pascal1), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "vandermonde", 
           FUNCTION(rb_gsl_matrix,vandermonde), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "vander",
           FUNCTION(rb_gsl_matrix,vandermonde), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "toeplitz", 
           FUNCTION(rb_gsl_matrix,toeplitz), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "circulant", 
           FUNCTION(rb_gsl_matrix,circulant), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "indgen",
           FUNCTION(rb_gsl_matrix,indgen_singleton), -1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "indgen", 
       FUNCTION(rb_gsl_matrix,indgen), -1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "indgen!", 
       FUNCTION(rb_gsl_matrix,indgen_bang), -1);

  rb_define_method(GSL_TYPE(cgsl_matrix), "info", 
       FUNCTION(rb_gsl_matrix,info), 0);
#ifdef BASE_DOUBLE
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "hilbert", 
           FUNCTION(rb_gsl_matrix,hilbert), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "hilb", 
           FUNCTION(rb_gsl_matrix,hilbert), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "invhilbert", 
           FUNCTION(rb_gsl_matrix,invhilbert), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "invhilb", 
           FUNCTION(rb_gsl_matrix,invhilbert), 1);
#endif

  rb_define_method(GSL_TYPE(cgsl_matrix), "any", 
       FUNCTION(rb_gsl_matrix,any), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "all", 
       FUNCTION(rb_gsl_matrix,all), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "rot90", 
       FUNCTION(rb_gsl_matrix,rot90), -1);

  rb_define_method(GSL_TYPE(cgsl_matrix), "diff", 
       FUNCTION(rb_gsl_matrix,diff), -1);
  rb_define_method(GSL_TYPE(cgsl_matrix), "isnan", FUNCTION(rb_gsl_matrix,isnan), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "isinf", FUNCTION(rb_gsl_matrix,isinf), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "finite", FUNCTION(rb_gsl_matrix,finite), 0);

  rb_define_method(GSL_TYPE(cgsl_matrix), "sgn", FUNCTION(rb_gsl_matrix,sgn), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "signum", "sgn");
  rb_define_method(GSL_TYPE(cgsl_matrix), "abs", FUNCTION(rb_gsl_matrix,abs), 0);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "fabs", "abs");

  rb_define_method(GSL_TYPE(cgsl_matrix), "horzcat", FUNCTION(rb_gsl_matrix,horzcat), 1);
  rb_define_alias(GSL_TYPE(cgsl_matrix), "cat", "horzcat");
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "horzcat", FUNCTION(rb_gsl_matrix,horzcat_singleton), 2);

  rb_define_method(GSL_TYPE(cgsl_matrix), "vertcat", FUNCTION(rb_gsl_matrix,vertcat), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_matrix), "vertcat", FUNCTION(rb_gsl_matrix,vertcat_singleton), 2);

#ifdef GSL_1_9_LATER
  rb_define_method(GSL_TYPE(cgsl_matrix), "ispos", FUNCTION(rb_gsl_matrix,ispos), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "ispos?", FUNCTION(rb_gsl_matrix,ispos2), 0); 
  rb_define_method(GSL_TYPE(cgsl_matrix), "isneg", FUNCTION(rb_gsl_matrix,isneg), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "isneg?", FUNCTION(rb_gsl_matrix,isneg2), 0);    
#endif

#ifdef GSL_1_10_LATER
  rb_define_method(GSL_TYPE(cgsl_matrix), "isnonneg", FUNCTION(rb_gsl_matrix,isnonneg), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "isnonneg?", FUNCTION(rb_gsl_matrix,isnonneg2), 0);
#endif

  rb_define_method(GSL_TYPE(cgsl_matrix), "symmetrize", FUNCTION(rb_gsl_matrix,symmetrize), 0);
  rb_define_method(GSL_TYPE(cgsl_matrix), "symmetrize!", FUNCTION(rb_gsl_matrix,symmetrize_bang), 0);
}

#undef NUMCONV
#undef NUMCONV2
#undef PRINTF_FORMAT
#undef VEC_ROW_COL
#undef VEC_P
#undef C_TO_VALUE
#undef C_TO_VALUE2
#undef VEC_COL_P
#undef VEC_ROW_P
#undef CHECK_VEC
#undef VEC_VIEW_P

#undef MAT_ROW_COL
#undef MAT_P
#undef MAT_COL_P
#undef MAT_ROW_P
#undef CHECK_MAT
#undef MAT_VIEW_P
