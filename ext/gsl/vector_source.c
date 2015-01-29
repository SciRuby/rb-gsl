/*
  vector_source.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2005 by Yoshiki Tsunesada
                               Cameron McBride

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#ifdef BASE_DOUBLE
#define NUMCONV(x) NUM2DBL(x)
#define NUMCONV2(x) NUM2DBL(x)
#define PRINTF_FORMAT "%5.3e "
#define VEC_ROW_COL VECTOR_ROW_COL
#define VEC_P VECTOR_P
#define VEC_ROW_P VECTOR_ROW_P
#define VEC_COL_P VECTOR_COL_P
#define C_TO_VALUE rb_float_new
#define C_TO_VALUE2 rb_float_new
#define CHECK_VEC CHECK_VECTOR
#define VEC_VIEW_P VECTOR_VIEW_P
#elif defined(BASE_INT)
#define NUMCONV(x) FIX2INT(x)
#define NUMCONV2(x) NUM2INT(x)
#define PRINTF_FORMAT "%d "
#define VEC_ROW_COL VECTOR_INT_ROW_COL
#define VEC_P VECTOR_INT_P
#define C_TO_VALUE INT2FIX
#define C_TO_VALUE2 INT2NUM
#define VEC_ROW_P VECTOR_INT_ROW_P
#define VEC_COL_P VECTOR_INT_COL_P
#define CHECK_VEC CHECK_VECTOR_INT
#define VEC_VIEW_P VECTOR_INT_VIEW_P
#endif

void FUNCTION(get_range,beg_en_n)(VALUE range, BASE *beg, BASE *en, size_t *n, int *step);

void get_range_beg_en_n_for_size(VALUE range,
    int *beg, int *en, size_t *n, int *step, size_t size);

void parse_subvector_args(int argc, VALUE *argv, size_t size,
    size_t *offset, size_t *stride, size_t *n);

void FUNCTION(get_range,beg_en_n)(VALUE range, BASE *beg, BASE *en, size_t *n, int *step)
{
  *beg = NUMCONV2(rb_funcall3(range, rb_gsl_id_beg, 0, NULL));
  *en = NUMCONV2(rb_funcall3(range, rb_gsl_id_end, 0, NULL));
  *n = (size_t) fabs(*en - *beg);
  if (!RTEST(rb_funcall3(range, rb_gsl_id_excl, 0, NULL))) *n += 1;
  if (*en < *beg) *step = -1; else *step = 1;
}

#ifdef BASE_INT
void get_range_beg_en_n_for_size(VALUE range, int *beg, int *en, size_t *n, int *step, size_t size)
{
  *beg = NUM2INT(rb_funcall3(range, rb_gsl_id_beg, 0, NULL));
  if(*beg < 0) *beg += size;
  *en = NUM2INT(rb_funcall3(range, rb_gsl_id_end, 0, NULL));
  if(*en < 0) *en += size;
  *n = (size_t) fabs(*en - *beg);
  if (!RTEST(rb_funcall3(range, rb_gsl_id_excl, 0, NULL))) *n += 1;
  if (*en < *beg) *step = -1; else *step = 1;
}

void parse_subvector_args(int argc, VALUE *argv, size_t size,
    size_t *offset, size_t *stride, size_t *n)
{
  int begin = 0, end, step, length;
  *stride = 1;
  switch (argc) {
  case 0:
    *n = size;
    break;
  case 1:
    if(rb_obj_is_kind_of(argv[0], rb_cRange)) {
      get_range_beg_en_n_for_size(argv[0], &begin, &end, n, &step, size);
      // TODO Should we do bounds checking or risk letting GSL do it?
      // On one hand, it seems like we should do as little as possible to stay as
      // thin and fast as possible.  On the other hand, it seems like we don't
      // want to let Ruby crash if GSL does not have bounds checking enabled.
      if(begin < 0 || (size_t)begin >= size) {
        rb_raise(rb_eRangeError,
            "begin value %d is out of range for Vector of length %d",
     begin, (int) size);
      }
      if(end < 0 || (size_t)end >= size) {
        rb_raise(rb_eRangeError,
            "end value %d is out of range for Vector of length %d",
     end, (int) size);
      }
      *stride = (size_t)step;
    } else {
      CHECK_FIXNUM(argv[0]);
      length = FIX2INT(argv[0]);
      if((length < 0 && -length > (int) size) || (length > 0 && length > (int) size)) {
        rb_raise(rb_eRangeError,
            "length %d is out of range for Vector of length %d",
     length, (int) size);
      } else if(length < 0) {
        begin = length;
        *n = (size_t)(-length);
      } else {
        // begin was initialized to 0
        *n = (size_t)length;
      }
    }
    break;
  case 2:
    if(rb_obj_is_kind_of(argv[0], rb_cRange)) {
      get_range_beg_en_n_for_size(argv[0], &begin, &end, n, &step, size);
      if(begin < 0 || (size_t)begin >= size) {
        rb_raise(rb_eRangeError,
            "begin value %d is out of range for Vector of length %d",
     (int) begin, (int) size);
      }
      if(end < 0 || (size_t)end >= size) {
        rb_raise(rb_eRangeError,
            "end value %d is out of range for Vector of length %d",
     (int) end, (int) size);
      }
      CHECK_FIXNUM(argv[1]);
      step = FIX2INT(argv[1]);
      if(step == 0 && begin != end) {
        rb_raise(rb_eArgError, "stride must be non-zero");
      } else if((step < 0 && begin <= end) || (step > 0 && end < begin)) {
        step = -step;
      }
      if(step < 0) {
        *n = (*n-1)/(-step) + 1;
      } else if(step > 0) {
        *n = (*n-1)/step + 1;
      }
      *stride = (size_t)step;
    } else {
      CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
      begin = FIX2INT(argv[0]);
      length = FIX2INT(argv[1]);
      if(length < 0) {
        length = -length;
        *stride = -1;
      }
      *n = (size_t)length;
    }
    break;
  case 3:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
    begin = FIX2INT(argv[0]);
    step = FIX2INT(argv[1]);
    length = FIX2INT(argv[2]);
    if(length < 0) {
      step = -step;
      length = -length;
    }
    *stride = (size_t)step;
    *n = (size_t)length;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0-3)", argc);
    break;
  }
  if(begin < 0) {
    begin += size;
  }
  *offset = (size_t)begin;
}

#endif

void FUNCTION(set_ptr_data,by_range)(BASE *ptr, size_t n, VALUE range)
{
  size_t n2, i;
  BASE beg, en, val;
  int step;
  FUNCTION(get_range,beg_en_n)(range, &beg, &en, &n2, &step);
  val = beg;
  for (i = 0; i < n; i++) {
    if (i < n2) ptr[i] = val;
    else ptr[i] = (BASE) 0;
    val += step;
  }
}

void FUNCTION(cvector,set_from_rarray)(GSL_TYPE(gsl_vector) *v, VALUE ary)
{
  size_t i;
  if (CLASS_OF(ary) == rb_cRange) ary = rb_gsl_range2ary(ary);
  Check_Type(ary, T_ARRAY);
  if (RARRAY_LEN(ary) == 0) return;
  for (i = 0; i < v->size; i++) FUNCTION(gsl_vector,set)(v, i, NUMCONV(rb_ary_entry(ary, i)));
}

GSL_TYPE(gsl_vector)* FUNCTION(make_cvector,from_rarray)(VALUE ary)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  if (CLASS_OF(ary) == rb_cRange) ary = rb_gsl_range2ary(ary);
  Check_Type(ary, T_ARRAY);
  v = FUNCTION(gsl_vector,alloc)(RARRAY_LEN(ary));
  if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
  FUNCTION(cvector,set_from_rarray)(v, ary);
  return v;
}

VALUE FUNCTION(rb_gsl_vector,new)(int argc, VALUE *argv, VALUE klass)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vtmp = NULL;
  BASE xnative;
  size_t n, i;
  BASE beg, en;
  int step;
#ifdef HAVE_NARRAY_H
  VALUE ary2;
#endif
  switch (argc) {
  case 1:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0])) {
      n = NA_TOTAL(argv[0]);
      v = FUNCTION(gsl_vector,alloc)(n);
      if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
#ifdef BASE_DOUBLE
      ary2 = na_change_type(argv[0], NA_DFLOAT);
#else
      ary2 = na_change_type(argv[0], NA_LINT);
#endif
      memcpy(v->data, NA_PTR_TYPE(ary2,BASE*), n*sizeof(BASE));
      return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_vector,free), v);
    }
#endif
    switch (TYPE(argv[0])) {
    case T_FIXNUM:  
      /*! if an integer n is given, create an empty vector of length n */
      CHECK_FIXNUM(argv[0]);
      n = FIX2INT(argv[0]);
      v = FUNCTION(gsl_vector,calloc)(n);
      if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
      break;
    case T_BIGNUM:
      rb_raise(rb_eRangeError, "vector length is limited within the range of Fixnum.");
      break;
    case T_FLOAT:
      v = FUNCTION(gsl_vector,alloc)(1);
      switch(TYPE(argv[0])) {
        case T_FIXNUM: case T_BIGNUM: case T_FLOAT:
          xnative = NUMCONV2(argv[0]);
          break;
        default:
          xnative = (BASE)0;
      }
      FUNCTION(gsl_vector,set)(v, 0, xnative);
      break;
#ifdef BASE_DOUBLE
    case T_ARRAY: 
      v = make_cvector_from_rarrays(argv[0]);
      break;
#endif
    default:
      if (CLASS_OF(argv[0]) == rb_cRange) {
  FUNCTION(get_range,beg_en_n)(argv[0], &beg, &en, &n, &step);
  v = FUNCTION(gsl_vector,alloc)(n);
  FUNCTION(set_ptr_data,by_range)(v->data, v->size, argv[0]);
  return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_vector,free), v);
      } else if (VEC_P(argv[0])) {
  /*! Create a new vector with the same elements of the vector given */ 
  Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), vtmp);
  v = FUNCTION(gsl_vector,alloc)(vtmp->size);
  for (i = 0; i < vtmp->size; i++) 
    FUNCTION(gsl_vector,set)(v, i, FUNCTION(gsl_vector,get)(vtmp, i));
  return Data_Wrap_Struct(VEC_ROW_COL(argv[0]), 0, FUNCTION(gsl_vector,free), v);
      } else {
  rb_raise(rb_eTypeError, 
     "wrong argument type %s", rb_class2name(CLASS_OF(argv[0])));
      }
      break;
    }
    break;
  default:
    v = FUNCTION(gsl_vector,alloc)(argc);
    if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
    for (i = 0; (int) i < argc; i++) {
      switch(TYPE(argv[i])) {
        case T_FIXNUM: case T_BIGNUM: case T_FLOAT:
          xnative = NUMCONV2(argv[i]);
          break;
        default:
          xnative = (BASE)0;
      }
      FUNCTION(gsl_vector,set)(v, i, xnative);
    }
    break;
  }
  return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_vector,free), v);
}

static VALUE FUNCTION(rb_gsl_vector,calloc)(VALUE klass, VALUE nn)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  CHECK_FIXNUM(nn);
  v = FUNCTION(gsl_vector,calloc)(FIX2INT(nn));
  if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_calloc failed");
  return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_vector,free), v);
}

static VALUE FUNCTION(rb_gsl_vector,subvector)(int argc, VALUE *argv, VALUE obj);
static VALUE FUNCTION(rb_gsl_vector,get)(int argc, VALUE *argv, VALUE obj)
{
  VALUE retval = Qnil;
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  //  QUALIFIED_VIEW(gsl_vector,view) *vv;
  gsl_index *p;
  int i;  /*! not size_t, since a negative index is allowed */
  size_t j, k;
  // If argc is not 1 or argv[0] is a Range
  if( argc != 1 || rb_obj_is_kind_of(argv[0], rb_cRange)) {
    // Treat as call to subvector
    retval = FUNCTION(rb_gsl_vector,subvector)(argc, argv, obj);
  } else {
    Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);

    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      i = FIX2INT(argv[0]);
      if (i < 0) 
        retval = C_TO_VALUE2(FUNCTION(gsl_vector,get)(v, (size_t) (v->size + i)));
      else
        retval = C_TO_VALUE2(FUNCTION(gsl_vector,get)(v, (size_t) (i)));
      break;
    case T_ARRAY:
      vnew = FUNCTION(gsl_vector,alloc)(RARRAY_LEN(argv[0]));
      for (j = 0; j < vnew->size; j++) {
        i = NUMCONV(rb_ary_entry(argv[0], j));
        if (i < 0) k = v->size + i;
        //  else k = j;
        else k = i;
        FUNCTION(gsl_vector,set)(vnew, j, FUNCTION(gsl_vector,get)(v, k));
      }
      retval = Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), vnew);
      break;
    default:
      if (PERMUTATION_P(argv[0])) {
        Data_Get_Struct(argv[0], gsl_index, p);
        vnew = FUNCTION(gsl_vector,alloc)(p->size);
        for (j = 0; j < p->size; j++) {
          k = p->data[j];
          //if (k < 0) k = p->size + j;
          FUNCTION(gsl_vector,set)(vnew, j, FUNCTION(gsl_vector,get)(v, k));
        }
        retval = Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), vnew);
      } else {
        // TODO Support Vector::Int (and even Vector?)
        rb_raise(rb_eTypeError, "wrong argument type %s (Array, Range, GSL::Permutation, or Fixnum expected)", rb_class2name(CLASS_OF(argv[0])));
      }
      break;
    }
  }
  return retval;
}

static VALUE FUNCTION(rb_gsl_vector,size)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return INT2FIX(v->size);
}

static VALUE FUNCTION(rb_gsl_vector,stride)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return INT2FIX(v->stride);
}

static VALUE FUNCTION(rb_gsl_vector,set_stride)(VALUE obj, VALUE ss)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  v->stride = (size_t) FIX2INT(ss);
  return obj;
}

static VALUE FUNCTION(rb_gsl_vector,owner)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return INT2FIX(v->owner);
}

void FUNCTION(rb_gsl_vector,set_subvector)(int argc, VALUE *argv, GSL_TYPE(gsl_vector) *v, VALUE other)
{
  GSL_TYPE(gsl_vector) *vother;
  QUALIFIED_VIEW(gsl_vector,view) vv;
  int step;
  size_t i, offset, stride, n, nother;
  BASE beg, end;

  // assignment to v.subvector(...)
  parse_subvector_args(argc, argv, v->size, &offset, &stride, &n);
  vv = FUNCTION(gsl_vector,subvector_with_stride)(v, offset, stride, n);
  if(rb_obj_is_kind_of(other, GSL_TYPE(cgsl_vector))) {
    Data_Get_Struct(other, GSL_TYPE(gsl_vector), vother);
    if(n != vother->size) {
      rb_raise(rb_eRangeError, "lengths do not match (%d != %d)",(int) n, (int) vother->size);
    }
    // TODO Change to gsl_vector_memmove if/when GSL has such a function
    // because gsl_vector_memcpy does not handle overlapping regions (e.g.
    // Views) well.
    FUNCTION(gsl_vector,memcpy)(&vv.vector, vother);
  } else if(rb_obj_is_kind_of(other, rb_cArray)) {
    if((int) n != RARRAY_LEN(other)) {
      rb_raise(rb_eRangeError, "lengths do not match (%d != %d)", (int) n, (int) RARRAY_LEN(other));
    }
    for(i = 0; i < n; i++) {
      FUNCTION(gsl_vector,set)(&vv.vector, i, NUMCONV2(rb_ary_entry(other, i)));
    }
  } else if(rb_obj_is_kind_of(other, rb_cRange)) {
    FUNCTION(get_range,beg_en_n)(other, &beg, &end, &nother, &step);
    if(n != nother) {
      rb_raise(rb_eRangeError, "lengths do not match (%d != %d)", (int) n, (int) nother);
    }
    for(i = 0; i < n; i++) {
      FUNCTION(gsl_vector,set)(&vv.vector, i, beg);
      beg += step;
    }
  } else {
    FUNCTION(gsl_vector,set_all)(&vv.vector, NUMCONV2(other));
  }
}

static VALUE FUNCTION(rb_gsl_vector,set)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v;
  VALUE other;
  int ii;

  if(argc < 1 || argc > 4) {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1-4)", argc);
  }

  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  other = argv[argc-1];

  if(argc == 1) {
    // // If assigning from another vector
    if(VECTOR_P(other) || VECTOR_INT_P(other)) {
      // treat as assignment to v.subvector(...)
      FUNCTION(rb_gsl_vector,set_subvector)(argc-1, argv, v, other);
    } else {
      FUNCTION(gsl_vector,set_all)(v, NUMCONV2(other));
    }
  } else if(argc == 2 && TYPE(argv[0]) == T_FIXNUM) {
    // v[i] = x
    ii = FIX2INT(argv[0]);
    if(ii < 0) ii += v->size;
    FUNCTION(gsl_vector,set)(v, (size_t)ii, NUMCONV2(other));
  } else {
    // assignment to v.subvector(...)
    FUNCTION(rb_gsl_vector,set_subvector)(argc-1, argv, v, other);
  }

  return obj;
}

static VALUE FUNCTION(rb_gsl_vector,set_all)(VALUE obj, VALUE xx)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  BASE xnative = NUMCONV2(xx);
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(gsl_vector,set_all)(v, xnative);
  return obj;
}

VALUE FUNCTION(rb_gsl_vector,do_something)(VALUE obj, void (*func)(GSL_TYPE(gsl_vector)*))
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  (*func)(v);
  return obj;
}

static VALUE FUNCTION(rb_gsl_vector,set_zero)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,do_something)(obj, FUNCTION(gsl_vector,set_zero));
}

static VALUE FUNCTION(rb_gsl_vector,set_basis)(VALUE obj, VALUE ii)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  CHECK_FIXNUM(ii);
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(gsl_vector,set_basis)(v, (size_t) FIX2INT(ii));
  return obj;
}

static VALUE FUNCTION(rb_gsl_vector,each)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  for (i = 0; i < v->size; i++) rb_yield(C_TO_VALUE2(FUNCTION(gsl_vector,get)(v, i)));
  return Qnil;
}

static VALUE FUNCTION(rb_gsl_vector,reverse_each)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  for (i = v->size-1;; i--) {
    rb_yield(C_TO_VALUE2(FUNCTION(gsl_vector,get)(v, i)));
    if (i == 0) break;
  }
  return Qnil;
}

static VALUE FUNCTION(rb_gsl_vector,each_index)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  for (i = 0; i < v->size; i++) rb_yield(INT2FIX(i));
  return Qnil;
}

static VALUE FUNCTION(rb_gsl_vector,reverse_each_index)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  for (i = v->size-1;; i--) {
    rb_yield(INT2FIX(i));
    if (i == 0) break;
  }
  return Qnil;
}

static VALUE FUNCTION(rb_gsl_vector,to_a)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t i;
  VALUE ary;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  ary = rb_ary_new2(v->size);
  for (i = 0; i < v->size; i++)
    rb_ary_store(ary, i, C_TO_VALUE2(FUNCTION(gsl_vector,get)(v, i)));
  return ary;
}

static VALUE FUNCTION(rb_gsl_vector,reverse_bang)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(gsl_vector,reverse)(v);
  return obj;
}

static VALUE FUNCTION(rb_gsl_vector,reverse)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  if (vnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_int_alloc failed");
  FUNCTION(gsl_vector,memcpy)(vnew, v);
  FUNCTION(gsl_vector,reverse)(vnew);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), vnew);
}

static VALUE FUNCTION(rb_gsl_vector,max)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return C_TO_VALUE2(FUNCTION(gsl_vector,max)(v));
}

static VALUE FUNCTION(rb_gsl_vector,min)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return C_TO_VALUE2(FUNCTION(gsl_vector,min)(v));
}

static VALUE FUNCTION(rb_gsl_vector,minmax)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  BASE min, max;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(gsl_vector,minmax)(v, &min, &max);
  return rb_ary_new3(2, C_TO_VALUE2(min), C_TO_VALUE2(max));
}

static VALUE FUNCTION(rb_gsl_vector,maxmin)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  BASE min, max;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(gsl_vector,minmax)(v, &min, &max);
  return rb_ary_new3(2, C_TO_VALUE2(max), C_TO_VALUE2(min));
}

static VALUE FUNCTION(rb_gsl_vector,max_index)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return INT2FIX(FUNCTION(gsl_vector,max_index)(v));
}

static VALUE FUNCTION(rb_gsl_vector,min_index)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return INT2FIX(FUNCTION(gsl_vector,min_index)(v));
}

static VALUE FUNCTION(rb_gsl_vector,minmax_index)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t imin, imax;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(gsl_vector,minmax_index)(v, &imin, &imax);
  return rb_ary_new3(2, INT2FIX(imin), INT2FIX(imax));
}

static VALUE FUNCTION(rb_gsl_vector,maxmin_index)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t imin, imax;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(gsl_vector,minmax_index)(v, &imin, &imax);
  return rb_ary_new3(2, INT2FIX(imax), INT2FIX(imin));
}

static VALUE FUNCTION(rb_gsl_vector,isnull)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return INT2FIX(FUNCTION(gsl_vector,isnull)(v));
}

static VALUE FUNCTION(rb_gsl_vector,isnull2)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (FUNCTION(gsl_vector,isnull)(v)) return Qtrue;
  else return Qfalse;
}

static VALUE FUNCTION(rb_gsl_vector,trans)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(make_vector,clone)(v);
#ifdef BASE_DOUBLE
  if (VECTOR_COL_P(obj))
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
  else return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, vnew);
#elif defined(BASE_INT)
  if (VECTOR_INT_ROW_P(obj))
    return Data_Wrap_Struct(cgsl_vector_int_col, 0, gsl_vector_int_free, vnew);
  else if (VECTOR_INT_COL_P(obj))
    return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, vnew);
  else rb_raise(rb_eTypeError, 
    "wrong type %s (Vector::Int or Vector::Int::Col expected)",
    rb_class2name(CLASS_OF(obj)));
#endif
  return Qnil;
}

static VALUE FUNCTION(rb_gsl_vector,trans_bang)(VALUE obj)
{
#ifdef BASE_DOUBLE
  if (CLASS_OF(obj) == cgsl_vector) RBGSL_SET_CLASS(obj, cgsl_vector_col);
  else if (CLASS_OF(obj) == cgsl_vector_col) RBGSL_SET_CLASS(obj, cgsl_vector);
  else {
    rb_raise(rb_eRuntimeError, "method trans! for %s is not permitted.",
       rb_class2name(CLASS_OF(obj)));
  }  
#elif defined(BASE_INT)
  if (CLASS_OF(obj) == cgsl_vector_int) RBGSL_SET_CLASS(obj, cgsl_vector_int_col);
  else if (CLASS_OF(obj) == cgsl_vector_int_col) RBGSL_SET_CLASS(obj, cgsl_vector_int);
  else {
    rb_raise(rb_eRuntimeError, "method trans! for %s is not permitted.",
       rb_class2name(CLASS_OF(obj)));
  }
#endif
  return obj;
}

static VALUE FUNCTION(rb_gsl_vector,uplus)(VALUE obj)
{
  return obj;
}

EXTERN VALUE cgsl_poly;

VALUE FUNCTION(rb_gsl_vector,uminus)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  for (i = 0; i < v->size; i++) {
    FUNCTION(gsl_vector,set)(vnew, i, -FUNCTION(gsl_vector,get)(v, i));
  }
  if (CLASS_OF(obj) == GSL_TYPE(cgsl_poly))
    return Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), vnew);
  else 
    return Data_Wrap_Struct(VEC_ROW_COL(obj), 0, FUNCTION(gsl_vector,free), vnew);
}

static VALUE FUNCTION(rb_gsl_vector,sum)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  BASE sum = 0;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  for (i = 0; i < v->size; i++) sum += FUNCTION(gsl_vector,get)(v, i);
  return C_TO_VALUE2(sum);
}

/* Vector#sumsq is defined in blas1.c */
#ifdef BASE_INT
static VALUE FUNCTION(rb_gsl_vector,sumsq)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  BASE sum = 0, x;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  for (i = 0; i < v->size; i++) {
    x = FUNCTION(gsl_vector,get)(v, i);
    sum += x*x;
  }
  return C_TO_VALUE2(sum);
}
#endif

static VALUE FUNCTION(rb_gsl_vector,prod)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  BASE x = 1;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  for (i = 0; i < v->size; i++) x *= FUNCTION(gsl_vector,get)(v, i);
  return C_TO_VALUE(x);
}

static VALUE FUNCTION(rb_gsl_vector,connect)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  BASE *ptr = NULL;
  size_t i, total = 0;
  if (VEC_P(obj)) {
    Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
    total += v->size;
  }
  for (i = 0; (int) i < argc; i++) {
    CHECK_VEC(argv[i]);
    Data_Get_Struct(argv[i], GSL_TYPE(gsl_vector), v);
    total += v->size;
  }
  vnew = FUNCTION(gsl_vector,alloc)(total);
  ptr = vnew->data;
  if (VEC_P(obj)) {
    Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
    memcpy(ptr, v->data, sizeof(BASE)*v->size);
    ptr += v->size;
  }
  for (i = 0; (int) i < argc; i++) {
    Data_Get_Struct(argv[i], GSL_TYPE(gsl_vector), v);
    memcpy(ptr, v->data, sizeof(BASE)*v->size);
    ptr += v->size;
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), vnew);
}

GSL_TYPE(gsl_vector)* FUNCTION(mygsl_vector,up)(GSL_TYPE(gsl_vector) *p)
{
  GSL_TYPE(gsl_vector) *pnew;
  pnew = FUNCTION(gsl_vector,alloc)(p->size + 1);
  FUNCTION(gsl_vector,set)(pnew, 0, 0);
  memcpy(pnew->data+1, p->data, sizeof(BASE)*p->size);
  return pnew;
}

void FUNCTION(mygsl_vector,up2)(GSL_TYPE(gsl_vector) *pnew, GSL_TYPE(gsl_vector) *p)
{
  FUNCTION(gsl_vector,set_all)(pnew, 0);
  memcpy(pnew->data+1, p->data, sizeof(BASE)*p->size);
}

GSL_TYPE(gsl_vector)* FUNCTION(mygsl_vector,down)(GSL_TYPE(gsl_vector) *p)
{
  GSL_TYPE(gsl_vector) *pnew;
  if (p->size <= 1) {
    rb_raise(rb_eRangeError, "Length <= 1, cannot be shortened.");
  }
  pnew = FUNCTION(gsl_vector,alloc)(p->size - 1);
  memcpy(pnew->data, p->data + 1, sizeof(BASE)*(p->size-1));
  return pnew;
}

static VALUE FUNCTION(rb_gsl_vector,sgn)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  BASE x;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  for (i = 0; i < v->size; i++) {
    x = FUNCTION(gsl_vector,get)(v, i);
    FUNCTION(gsl_vector,set)(vnew, i, (BASE)(x>0 ? 1 : (x<0 ? -1 : 0)));
  }
  return Data_Wrap_Struct(VEC_ROW_COL(obj), 0, FUNCTION(gsl_vector,free), vnew);  
}

static VALUE FUNCTION(rb_gsl_vector,abs)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  for (i = 0; i < v->size; i++) {
    FUNCTION(gsl_vector,set)(vnew, i, (BASE) fabs(FUNCTION(gsl_vector,get)(v, i)));
  }
  return Data_Wrap_Struct(VEC_ROW_COL(obj), 0, FUNCTION(gsl_vector,free), vnew);  
}

static VALUE FUNCTION(rb_gsl_vector,square)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  for (i = 0; i < v->size; i++) {
    FUNCTION(gsl_vector,set)(vnew, i, gsl_pow_2(FUNCTION(gsl_vector,get)(v, i)));
  }
  return Data_Wrap_Struct(VEC_ROW_COL(obj), 0, FUNCTION(gsl_vector,free), vnew);  
}

static VALUE FUNCTION(rb_gsl_vector,sqrt)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  for (i = 0; i < v->size; i++) {
    FUNCTION(gsl_vector,set)(vnew, i, sqrt(FUNCTION(gsl_vector,get)(v, i)));
  }
  return Data_Wrap_Struct(VEC_ROW_COL(obj), 0, FUNCTION(gsl_vector,free), vnew);  
}

static VALUE FUNCTION(rb_gsl_vector,memcpy)(VALUE obj, VALUE dest, VALUE src)
{
  GSL_TYPE(gsl_vector) *vdest = NULL, *vsrc = NULL;
  Data_Get_Struct(dest, GSL_TYPE(gsl_vector), vdest);
  Data_Get_Struct(src, GSL_TYPE(gsl_vector), vsrc);
  FUNCTION(gsl_vector,memcpy)(vdest, vsrc);
  return dest;
}

static VALUE FUNCTION(rb_gsl_vector,clone)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  if (vnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
  FUNCTION(gsl_vector,memcpy)(vnew, v);
  if (!VEC_VIEW_P(obj))
    return Data_Wrap_Struct(CLASS_OF(obj), 0, FUNCTION(gsl_vector,free), vnew);
  else 
    return Data_Wrap_Struct(VEC_ROW_COL(obj), 0, FUNCTION(gsl_vector,free), vnew);
}

/* singleton */
static VALUE FUNCTION(rb_gsl_vector,swap)(VALUE obj, VALUE vv, VALUE ww)
{
  GSL_TYPE(gsl_vector) *v = NULL, *w = NULL;
  Data_Get_Struct(vv, GSL_TYPE(gsl_vector), v);
  Data_Get_Struct(ww, GSL_TYPE(gsl_vector), w);
  FUNCTION(gsl_vector,swap)(v, w);
  return obj;
}

static VALUE FUNCTION(rb_gsl_vector,swap_elements)(VALUE obj, VALUE i, VALUE j)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  CHECK_FIXNUM(i);  CHECK_FIXNUM(j);
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(gsl_vector,swap_elements)(v, FIX2INT(i), FIX2INT(j));
  return obj;
}

static VALUE FUNCTION(rb_gsl_vector,fwrite)(VALUE obj, VALUE io)
{
  GSL_TYPE(gsl_vector) *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), h);
  f = rb_gsl_open_writefile(io, &flag);
  status = FUNCTION(gsl_vector,fwrite)(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_gsl_vector,fread)(VALUE obj, VALUE io)
{
  GSL_TYPE(gsl_vector) *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), h);
  f = rb_gsl_open_readfile(io, &flag);
  status = FUNCTION(gsl_vector,fread)(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_gsl_vector,fprintf)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  if (argc != 1 && argc != 2) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), h);
  fp = rb_gsl_open_writefile(argv[0], &flag);
  if (argc == 2) {
    if (TYPE(argv[1]) == T_STRING)
      status = FUNCTION(gsl_vector,fprintf)(fp, h, STR2CSTR(argv[1]));
    else
      rb_raise(rb_eTypeError, "argv 2 String expected");
  } else {
    status = FUNCTION(gsl_vector,fprintf)(fp, h, "%g");
  }
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_gsl_vector,printf)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *h = NULL;
  int status;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), h);
  if (argc == 1) {
    if (TYPE(argv[0]) != T_STRING) 
      rb_raise(rb_eTypeError, "String expected");
    else
      status = FUNCTION(gsl_vector,fprintf)(stdout, h, STR2CSTR(argv[0]));
  } else {
    status = FUNCTION(gsl_vector,fprintf)(stdout, h, "%g");
  }
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_gsl_vector,fscanf)(VALUE obj, VALUE io)
{
  GSL_TYPE(gsl_vector) *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), h);
  fp = rb_gsl_open_readfile(io, &flag);
  status = FUNCTION(gsl_vector,fscanf)(fp, h);
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

/* 2.Aug.2004 */
VALUE FUNCTION(rb_gsl_vector,inner_product)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *v2 = NULL;
  BASE prod = 0;
#ifndef BASE_DOUBLE
  size_t i;
#endif
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    CHECK_VEC(argv[0]);
    CHECK_VEC(argv[1]);
    Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), v);
    Data_Get_Struct(argv[1], GSL_TYPE(gsl_vector), v2);
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    CHECK_VEC(argv[0]);
    Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
    Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), v2);
    break;
  }
  if (v->size != v2->size) rb_raise(rb_eRangeError, "vector lengths are different.");
#ifdef BASE_DOUBLE
  gsl_blas_ddot(v, v2, &prod);
#else
  for (i = 0; i < v->size; i++) {
    prod += FUNCTION(gsl_vector,get)(v, i)*FUNCTION(gsl_vector,get)(v2, i);
  }
#endif
  return C_TO_VALUE2(prod);
}

int FUNCTION(rbgsl_vector,equal)(const GSL_TYPE(gsl_vector) *v1, const GSL_TYPE(gsl_vector) *v2, double eps)
{
  size_t i;
  BASE x, y;
  if (v1->size != v2->size) return 0;
  for (i = 0; i < v2->size; i++) {
    x = FUNCTION(gsl_vector,get)(v1, i);
    y = FUNCTION(gsl_vector,get)(v2, i);
    if (fabs(x - y) > eps) return 0;
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

static VALUE FUNCTION(rb_gsl_vector,equal)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v1, *v2;
  VALUE other;
  size_t i;
  double eps = 1e-10;
  double x;
  switch (argc) {
  case 2:
    other = argv[0];
    eps = NUM2DBL(argv[1]);
    break;
  case 1:
    other = argv[0];
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
#ifdef HAVE_TENSOR_TENSOR_H
  if (TEN_P(other)) {
    return FUNCTION(rb_gsl_tensor,equal)(argc, argv, obj);
  }
#endif
  switch (TYPE(other)) {
  case T_FIXNUM:
  case T_FLOAT:
    x = NUM2DBL(other);
    Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v1);
    for (i = 0; i < v1->size; i++) 
      if (fabs(x-FUNCTION(gsl_vector,get)(v1, i)) > eps) return Qfalse;
    return Qtrue;
    break;
  default:
    CHECK_VEC(other);
    Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v1);
    Data_Get_Struct(other, GSL_TYPE(gsl_vector), v2);
    if (FUNCTION(rbgsl_vector,equal)(v1, v2, eps)) return Qtrue;
    else return Qfalse;
    break;
  }
  return Qnil;
}

#ifdef HAVE_TENSOR_TENSOR_H
#ifdef TEN_P
#undef TEN_P
#endif
#endif

static VALUE FUNCTION(rb_gsl_vector,to_poly)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  GSL_TYPE(gsl_poly) *p = NULL;
  if (CLASS_OF(obj) == GSL_TYPE(cgsl_poly)) return obj;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  p = FUNCTION(make_vector,clone)(v);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), p);
}

static VALUE FUNCTION(rb_gsl_vector,graph)(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_GNU_GRAPH
  GSL_TYPE(gsl_vector) *x = NULL, *y = NULL;
  FILE *fp = NULL;
  size_t i;
  char command[1024];
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), y);
  switch (argc) {
  case 0:
    strcpy(command, "graph -T X -g 3");
    break;
  case 1:
    if (TYPE(argv[0]) == T_STRING) {
      make_graphcommand(command, argv[0]);
    } else if (VEC_P(argv[0])) {
      strcpy(command, "graph -T X -g 3");
      Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), x);
    } else {
    }
    break;
  case 2:
    if (TYPE(argv[1]) == T_STRING) {
      make_graphcommand(command, argv[1]);
      if (VEC_P(argv[0])) {
  Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), x);
      } else {
  rb_raise(rb_eTypeError, "argv[0] wrong type %s (String or Vector expected)", 
     rb_class2name(CLASS_OF(argv[0])));
      }
    } else {
      rb_raise(rb_eTypeError, "argv[1] wrong type %s (String or Vector expected)", 
         rb_class2name(CLASS_OF(argv[1])));
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    break;
  }
  if (y == NULL) rb_raise(rb_eRuntimeError, "ydata not given");
  fp = popen(command, "w");
  for (i = 0; i < y->size; i++) {
    if (x == NULL) 
      fprintf(fp, "%d %e\n", (int) i, (double) FUNCTION(gsl_vector,get)(y, i));
    else
      fprintf(fp, "%e %e\n", (double) FUNCTION(gsl_vector,get)(x, i), (double) FUNCTION(gsl_vector,get)(y, i));
  }
  fflush(fp);
  pclose(fp);
  fp = NULL;
  return Qtrue;
#else
  rb_raise(rb_eNoMethodError, "not implemented");
  return Qfalse;
#endif
}

static VALUE FUNCTION(rb_gsl_vector,graph_step)(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_GNU_GRAPH
  GSL_TYPE(gsl_vector) *x = NULL, *y = NULL;
  FILE *fp = NULL;
  size_t i;
  char command[1024];
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), y);
  switch (argc) {
  case 0:
    strcpy(command, "graph -T X -g 3");
    break;
  case 1:
    if (TYPE(argv[0]) == T_STRING) {
      make_graphcommand(command, argv[0]);
    } else if (VECTOR_P(argv[0])) {
      strcpy(command, "graph -T X -g 3");
      Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), x);
    } else {
    }
    break;
  case 2:
    if (TYPE(argv[1]) == T_STRING) {
      make_graphcommand(command, argv[1]);
      if (VEC_P(argv[0])) {
  Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), x);
      } else {
  rb_raise(rb_eTypeError, "argv[0] wrong type %s (String or Vector expected)", 
     rb_class2name(CLASS_OF(argv[0])));
      }
    } else {
      rb_raise(rb_eTypeError, "argv[1] wrong type %s (String or Vector expected)", 
         rb_class2name(CLASS_OF(argv[1])));
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    break;
  }
  if (y == NULL) rb_raise(rb_eRuntimeError, "ydata not given");
  fp = popen(command, "w");
  for (i = 0; i < y->size; i++) {
    if (x == NULL) {
      fprintf(fp, "%d %e\n%d %e\n", (int) i, (double) FUNCTION(gsl_vector,get)(y, i),
        (int) (i+1), (double) FUNCTION(gsl_vector,get)(y, i));
    } else {
      if (i != y->size-1) 
  fprintf(fp, "%e %e\n%e %e\n", (double) FUNCTION(gsl_vector,get)(x, i), 
  (double) FUNCTION(gsl_vector,get)(y, i),
  (double) FUNCTION(gsl_vector,get)(x, i+1), 
  (double) FUNCTION(gsl_vector,get)(y, i));
      else
  fprintf(fp, "%e %e\n%e %e", 
  (double) FUNCTION(gsl_vector,get)(x, i), 
  (double) FUNCTION(gsl_vector,get)(y, i),
  2.0*FUNCTION(gsl_vector,get)(x, i)-FUNCTION(gsl_vector,get)(x, i-1), 
  (double) FUNCTION(gsl_vector,get)(y, i));
    }
  }
  fflush(fp);
  pclose(fp);
  fp = NULL;
  return Qtrue;
#else
  rb_raise(rb_eNoMethodError, "not implemented");
  return Qfalse;
#endif
}

static VALUE FUNCTION(rb_gsl_vector,plot)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *x = NULL, *y = NULL;
  FILE *fp = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), y);
  fp = popen("gnuplot -persist", "w");
  switch (argc) {
  case 0:
    fprintf(fp, "plot '-'\n");
    break;
  case 1:
    if (TYPE(argv[0]) == T_STRING) {
      fprintf(fp, "plot '-' %s\n", STR2CSTR(argv[0]));
    } else if (VEC_P(argv[0])) {
      fprintf(fp, "plot '-'\n");
      Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), x);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (String or Vector expected)",
         rb_class2name(CLASS_OF(argv[0])));
    }
    break;
  case 2:
    if (TYPE(argv[1]) == T_STRING)
      fprintf(fp, "plot '-' %s\n", STR2CSTR(argv[1]));
    if (VEC_P(argv[0]))
      Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), x);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    break;
  }
  if (y == NULL) rb_raise(rb_eRuntimeError, "ydata not given");
  for (i = 0; i < y->size; i++) {
    if (x == NULL) 
      fprintf(fp, "%d %e\n", (int) i, (double) FUNCTION(gsl_vector,get)(y, i));
    else
      fprintf(fp, "%e %e\n", (double) FUNCTION(gsl_vector,get)(x, i), 
  (double) FUNCTION(gsl_vector,get)(y, i));
  }
  fprintf(fp, "e\n");
  fflush(fp);
  pclose(fp);
  fp = NULL;
  return Qtrue;
}

void FUNCTION(gsl_vector,print)(const GSL_TYPE(gsl_vector) *v, VALUE klass)
{
  size_t i;
  printf("[ ");
  if (klass == cgsl_vector_col || klass == cgsl_vector_col_view
  || klass == cgsl_vector_col_view_ro
  || klass == cgsl_vector_int_col || klass == cgsl_vector_int_col_view
  || klass == cgsl_vector_int_col_view_ro) {
    printf(PRINTF_FORMAT, FUNCTION(gsl_vector,get)(v, 0));
    for (i = 1; i < v->size; i++) {  
      printf(PRINTF_FORMAT, FUNCTION(gsl_vector,get)(v, i));
      if (i != v->size-1) printf("\n");
    }
  } else {
    for (i = 0; i < v->size; i++) printf(PRINTF_FORMAT, FUNCTION(gsl_vector,get)(v, i));
  }
  printf("]\n");
}

VALUE FUNCTION(rb_gsl_vector,print)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(gsl_vector,print)(v, CLASS_OF(obj));
  return Qnil;
}

#ifdef BASE_DOUBLE
#define SHOW_ELM 6
#else
#define SHOW_ELM 15
#endif

VALUE FUNCTION(rb_gsl_vector,to_s)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  char buf[32], format[32], format2[32];
  size_t i;
  VALUE str;
  BASE x;
  int dig = 8;
#ifdef BASE_INT
  BASE min;
  BASE max;
  dig = 1;
#endif
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (v->size == 0) return rb_str_new2("[ ]");
  str = rb_str_new2("[ ");
  if (VEC_COL_P(obj)) {
#ifdef BASE_INT
    min = FUNCTION(gsl_vector,min)(v);
    max = gsl_vector_int_max(v);
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
    for (i = 0; i < v->size; i++) {
      if (i != 0) {
  strcpy(buf, "  ");
  rb_str_cat(str, buf, strlen(buf));
      }
      x = FUNCTION(gsl_vector,get)(v, i);
      if (x < 0) sprintf(buf, format, x); 
      else sprintf(buf, format2, x); 
      if (i != v->size-1) strcat(buf, "\n");
      rb_str_cat(str, buf, strlen(buf));
      if (i >= 20 && i != v->size-1) {
        strcpy(buf, "  ...");
        rb_str_cat(str, buf, strlen(buf));
        break;
      }
    }
  } else {
    sprintf(buf,  PRINTF_FORMAT, FUNCTION(gsl_vector,get)(v, 0));
    rb_str_cat(str, buf, strlen(buf));
    for (i = 1; i < v->size; i++) {
      sprintf(buf,  PRINTF_FORMAT, FUNCTION(gsl_vector,get)(v, i));
      rb_str_cat(str, buf, strlen(buf));
      if ((int) i >= (55/dig) && i != v->size-1) {
        strcpy(buf, "... ");
        rb_str_cat(str, buf, strlen(buf));
        break;
      }
    }
  }
  sprintf(buf, "]");
  rb_str_cat(str, buf, strlen(buf));
  return str;
}
#undef SHOW_ELM

static VALUE FUNCTION(rb_gsl_vector,inspect)(VALUE obj)
{
  VALUE str;
  char buf[64];
  sprintf(buf, "%s\n", rb_class2name(CLASS_OF(obj)));
  str = rb_str_new2(buf);
  return rb_str_concat(str, FUNCTION(rb_gsl_vector,to_s)(obj));
}

static VALUE FUNCTION(rb_gsl_vector,subvector)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  QUALIFIED_VIEW(gsl_vector,view) *vv = NULL;
  size_t offset, stride, n;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  parse_subvector_args(argc, argv, v->size, &offset, &stride, &n);
  vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
  *vv = FUNCTION(gsl_vector,subvector_with_stride)(v, offset, stride, n);
  if (VEC_COL_P(obj))
    return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,col_view), 0, free, vv);
  else
    return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,view), 0, free, vv);
}

static VALUE FUNCTION(rb_gsl_vector,subvector_with_stride)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  QUALIFIED_VIEW(gsl_vector,view) *vv = NULL;
  int offset = 0, step, length;
  size_t stride = 1, n;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  switch (argc) {
  case 1:
    CHECK_FIXNUM(argv[0]);
    step = FIX2INT(argv[0]);
    if(step == 0) {
      rb_raise(rb_eArgError, "stride must be non-zero");
    }
    stride = (size_t)step;
    //n = v->size/stride;
    n = (v->size-1)/stride + 1;
    break;
  case 2:
    CHECK_FIXNUM(argv[0]);    CHECK_FIXNUM(argv[1]);
    offset = FIX2INT(argv[0]);
    step = FIX2INT(argv[1]);
    if(offset < 0) {
      offset += v->size;
      if(offset < 0) {
        rb_raise(rb_eRangeError, "offset %d out of range", offset - (int)v->size);
      }
    } else if(offset >= (int) v->size) {
      rb_raise(rb_eRangeError, "offset %d out of range", offset);
    }
    if(step == 0) {
      rb_raise(rb_eArgError, "stride must be non-zero");
    }
    stride = (size_t)step;
    //n = (v->size-(size_t)offset)/stride;
    n = (v->size-(size_t)offset-1)/stride + 1;
    break;
  case 3:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
    offset = FIX2INT(argv[0]);
    step = FIX2INT(argv[1]);
    length = FIX2INT(argv[2]);
    if(offset < 0) {
      offset += v->size;
      if(offset < 0) {
        rb_raise(rb_eRangeError, "offset %d out of range", offset - (int)v->size);
      }
    }
    if(step == 0) {
      rb_raise(rb_eArgError, "stride must be non-zero");
    }
    if(length < 0) {
      rb_raise(rb_eArgError, "length must be non-negative");
    }
    stride = (size_t)step;
    n = (size_t)length;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 - 3)", argc);
    break;
  }
  vv = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
  *vv = FUNCTION(gsl_vector,subvector_with_stride)(v, (size_t)offset, stride, n);
  if (VEC_COL_P(obj))
    return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,col_view), 0, free, vv);
  else
    return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,view), 0, free, vv);
}

static VALUE FUNCTION(rb_gsl_vector,matrix_view)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  QUALIFIED_VIEW(gsl_matrix,view) *mv = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  switch (argc) {
  case 2:
    mv = ALLOC(QUALIFIED_VIEW(gsl_matrix,view));
    *mv = FUNCTION(gsl_matrix,view_vector)(v, FIX2INT(argv[0]), FIX2INT(argv[1]));
    break;
  case 3:
    mv = ALLOC(QUALIFIED_VIEW(gsl_matrix,view));
    *mv = FUNCTION(gsl_matrix,view_vector_with_tda)(v, FIX2INT(argv[0]), FIX2INT(argv[1]), 
         FIX2INT(argv[2]));
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_matrix,view), 0, free, mv);
}

static VALUE FUNCTION(rb_gsl_vector,matrix_view_with_tda)(VALUE obj, VALUE nn1, VALUE nn2,
            VALUE tda)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  QUALIFIED_VIEW(gsl_matrix,view) *mv = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  mv = ALLOC(QUALIFIED_VIEW(gsl_matrix,view));
  *mv = FUNCTION(gsl_matrix,view_vector_with_tda)(v, FIX2INT(nn1), FIX2INT(nn2), FIX2INT(tda));
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_matrix,view), 0, free, mv);
}

void FUNCTION(mygsl_vector,shift)(GSL_TYPE(gsl_vector) *p, size_t n)
{
  size_t i;
  for (i = n;; i--) {
    FUNCTION(gsl_vector,set)(p, i+1, FUNCTION(gsl_vector,get)(p, i));
    if (i == 0) break;
  }
  FUNCTION(gsl_vector,set)(p, 0, 0);
}

void FUNCTION(mygsl_vector,shift_scale2)(GSL_TYPE(gsl_vector) *p, size_t n)
{
  size_t i;
  for (i = n;; i--) {
    FUNCTION(gsl_vector,set)(p, i+1, 2*FUNCTION(gsl_vector,get)(p, i));
    if (i == 0) break;
  }
  FUNCTION(gsl_vector,set)(p, 0, 0);
}

GSL_TYPE(gsl_vector)* FUNCTION(make_vector,clone)(const GSL_TYPE(gsl_vector) *v)
{
  GSL_TYPE(gsl_vector) *vnew = NULL;
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  if (v->stride == 1) memcpy(vnew->data, v->data, sizeof(BASE)*v->size);
  else FUNCTION(gsl_vector,memcpy)(vnew, v);
  return vnew;
}

VALUE FUNCTION(rb_gsl_vector,scale)(VALUE obj, VALUE x)
{
  GSL_TYPE(gsl_vector) *v, *vnew;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(make_vector,clone)(v);
  FUNCTION(gsl_vector,scale)(vnew, NUMCONV(x));
  //  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), vnew);
    return Data_Wrap_Struct(VEC_ROW_COL(obj), 0, FUNCTION(gsl_vector,free), vnew);
}

VALUE FUNCTION(rb_gsl_vector,scale_bang)(VALUE obj, VALUE x)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(gsl_vector,scale)(v, NUMCONV(x));
  return obj;
}

VALUE FUNCTION(rb_gsl_vector,add_constant)(VALUE obj, VALUE x)
{
  GSL_TYPE(gsl_vector) *v, *vnew;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(make_vector,clone)(v);
  FUNCTION(gsl_vector,add_constant)(vnew, NUMCONV(x));
  //  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), vnew);
  return Data_Wrap_Struct(VEC_ROW_COL(obj), 0, FUNCTION(gsl_vector,free), vnew);
}

VALUE FUNCTION(rb_gsl_vector,add_constant_bang)(VALUE obj, VALUE x)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(gsl_vector,add_constant)(v, NUMCONV(x));
  return obj;
}

QUALIFIED_VIEW(gsl_vector,view)* FUNCTION(rb_gsl_make_vector,view)(BASE *data, size_t size, size_t stride)
{
  QUALIFIED_VIEW(gsl_vector,view) *v = NULL;
  v = ALLOC(QUALIFIED_VIEW(gsl_vector,view));
  v->vector.size = size;
  v->vector.stride = stride;
  v->vector.owner = 0;
  v->vector.data = data;
  return v;
}

#ifdef HAVE_TENSOR_TENSOR_H
#include "include/rb_gsl_tensor.h"
static VALUE FUNCTION(rb_gsl_vector,to_tensor)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  GSL_TYPE(rbgsl_tensor) *t;
  unsigned int rank;
  size_t dim;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  switch (argc) {
  case 0:
    rank = 1;
    dim = v->size;
    break;
  case 2:
    rank = FIX2UINT(argv[0]);
    dim = FIX2UINT(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 2)", argc);
    break;
  }
  t = FUNCTION(rbgsl_tensor,alloc)(rank, dim);
  memcpy(t->tensor->data, v->data, sizeof(BASE)*v->size);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), t);
}
#endif

#ifdef BASE_DOUBLE
#define PRINTF_FORMAT2 "%g "
#else
#define PRINTF_FORMAT2 "%d "
#endif
static VALUE FUNCTION(rb_gsl_vector,to_gplot)(int argc, VALUE *argv, VALUE obj)
{
  char buf[1024] = "";
  size_t i, j, len = 0, nv, istart;
  VALUE str, tmp;
  GSL_TYPE(gsl_vector) *v, **vp;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "no vectors given");
    if (TYPE(argv[0]) == T_ARRAY) nv = RARRAY_LEN(argv[0]);
    else nv = argc;
    vp = (GSL_TYPE(gsl_vector)**) ALLOC_N(GSL_TYPE(gsl_vector)*, nv);
    istart = 0;
    break;
  default:
    CHECK_VEC(obj);
    Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
    if (argc >= 1 && TYPE(argv[0]) == T_ARRAY) nv = 1 + RARRAY_LEN(argv[0]);
    else nv = argc + 1;
    vp = (GSL_TYPE(gsl_vector)**) ALLOC_N(GSL_TYPE(gsl_vector)*, nv);
    vp[0] = v; len = v->size;
    istart = 1;
    break;
  }
  for (i = 0; (int) i < argc; i++) {
    if (TYPE(argv[0]) == T_ARRAY) tmp = rb_ary_entry(argv[0], i);
    else tmp = argv[i];
    CHECK_VEC(tmp);
    Data_Get_Struct(tmp, GSL_TYPE(gsl_vector), v);
    if (len == 0) len = v->size;
    if (len != v->size) 
      rb_raise(rb_eRuntimeError, "vectors must have equal lengths");
    vp[i+istart] = v;
  }
  str = rb_str_new2(buf);
  for (j = 0; j < len; j++) {
    for (i = 0; i < nv; i++) {
      sprintf(buf, PRINTF_FORMAT2, FUNCTION(gsl_vector,get)(vp[i], j));
      rb_str_buf_cat(str, buf, strlen(buf));
    }
    rb_str_buf_cat2(str, "\n");
  }
  rb_str_buf_cat2(str, "\n");
  free((GSL_TYPE(gsl_vector)**)vp);
  return str;
}
#undef PRINTF_FORMAT2

static VALUE FUNCTION(rb_gsl_vector,to_m_diagonal)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  GSL_TYPE(gsl_matrix) *m = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  m = FUNCTION(gsl_matrix,calloc)(v->size, v->size);
  for (i = 0; i < v->size; i++)
    FUNCTION(gsl_matrix,set)(m, i, i, FUNCTION(gsl_vector,get)(v, i));
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

static VALUE FUNCTION(rb_gsl_vector,collect)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  for (i = 0; i < v->size; i++) {
    FUNCTION(gsl_vector,set)(vnew, i, NUMCONV(rb_yield(C_TO_VALUE(FUNCTION(gsl_vector,get)(v, i)))));
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), vnew);
}

/* 2004/May/03 */
static VALUE FUNCTION(rb_gsl_vector,collect_bang)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  for (i = 0; i < v->size; i++) {
    FUNCTION(gsl_vector,set)(v, i, NUMCONV(rb_yield(C_TO_VALUE(FUNCTION(gsl_vector,get)(v, i)))));
  }
  return obj;
}

/* Modified 2006/Sep/26 */
GSL_TYPE(gsl_vector)* FUNCTION(mygsl_vector,mul_matrix)(GSL_TYPE(gsl_vector) *v,
              GSL_TYPE(gsl_matrix) *m)
{
  GSL_TYPE(gsl_vector) *vnew;
  size_t i, j;
  BASE sum;
  if (v->size != m->size1) rb_raise(rb_eRuntimeError, "vector/matrix sizes are different.");
  vnew = FUNCTION(gsl_vector,alloc)(m->size2);
  for (i = 0; i < m->size2; i++) {
    sum = 0;
    for (j = 0; j < m->size1; j++) {
      sum += FUNCTION(gsl_vector,get)(v, j)*FUNCTION(gsl_matrix,get)(m, j, i);
    }
    FUNCTION(gsl_vector,set)(vnew, i, sum);
  }
  return vnew;
}

void FUNCTION(mygsl_vector,to_m_circulant)(GSL_TYPE(gsl_matrix) *m, GSL_TYPE(gsl_vector) *v)
{
  size_t i, j;
  for (i = v->size-1;; i--) {
    for (j = 0; j < v->size; j++) {
      if (j <= i) FUNCTION(gsl_matrix,set)(m, i, j, FUNCTION(gsl_vector,get)(v, v->size-1-i+j));
      else FUNCTION(gsl_matrix,set)(m, i, j, FUNCTION(gsl_vector,get)(v, j-i-1));
    }
    if (i == 0) break;
  }
}

static VALUE FUNCTION(rb_gsl_vector,to_m_circulant)(VALUE obj)
{
  GSL_TYPE(gsl_matrix) *m;
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  m = FUNCTION(gsl_matrix,alloc)(v->size, v->size);
  FUNCTION(mygsl_vector,to_m_circulant)(m, v);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

static void FUNCTION(mygsl_vector,indgen)(GSL_TYPE(gsl_vector) *v,
            BASE start, BASE step)
{
  size_t k = 0;
  BASE i;
  i = start;
  for (k = 0; k < v->size; k++) {
    FUNCTION(gsl_vector,set)(v, k, i);
    i += step;
  }
}

static VALUE FUNCTION(rb_gsl_vector,indgen_singleton)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t n;
  BASE start = 0, step = 1;
  switch (argc) {
  case 3:
    step = NUMCONV2(argv[2]);
    /* no break */
  case 2:
    start = NUMCONV2(argv[1]);
    /* no break */
  case 1:
    n = NUM2INT(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1-3)", argc);
    break;
  }
  v = FUNCTION(gsl_vector,alloc)(n);
  FUNCTION(mygsl_vector,indgen)(v, start, step);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), v);
}

static VALUE FUNCTION(rb_gsl_vector,indgen)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v, *vnew;
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
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  FUNCTION(mygsl_vector,indgen)(vnew, start, step);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), vnew);
}

static VALUE FUNCTION(rb_gsl_vector,indgen_bang)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
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
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  FUNCTION(mygsl_vector,indgen)(v, start, step);
  return obj;
}

static VALUE FUNCTION(rb_gsl_vector,to_m)(VALUE obj, VALUE ii, VALUE jj)
{
  GSL_TYPE(gsl_matrix) *m;
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t i, j, n;
  CHECK_FIXNUM(ii); CHECK_FIXNUM(jj);
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  i = (size_t) FIX2INT(ii); j = (size_t) FIX2INT(jj);
  n = i*j;
  m = FUNCTION(gsl_matrix,alloc)(i, j);
  memcpy(m->data, v->data, sizeof(BASE)*v->size);
  for (i = n; i < v->size; i++) m->data[i] = (BASE) 0;
  return Data_Wrap_Struct(GSL_TYPE(cgsl_matrix), 0, FUNCTION(gsl_matrix,free), m);
}

static VALUE FUNCTION(rb_gsl_vector,block)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_block), 0, NULL, v->block);
}

/*****/
static VALUE GSL_TYPE(rb_gsl_sort_vector)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  GSL_TYPE(gsl_sort_vector)(v);
  return obj;
}

static VALUE GSL_TYPE(rb_gsl_sort_vector2)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  FUNCTION(gsl_vector,memcpy)(vnew, v);
  GSL_TYPE(gsl_sort_vector)(vnew);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), vnew);
}

static VALUE FUNCTION(rb_gsl_sort_vector,index)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_index *p = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  p = gsl_permutation_alloc(v->size);
  FUNCTION(gsl_sort_vector,index)(p, v);
  return Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, p);
}

static VALUE FUNCTION(rb_gsl_sort_vector,smallest)(VALUE obj, VALUE kk)
{
  GSL_TYPE(gsl_vector) *v = NULL, *v2 = NULL;
  size_t k;
  CHECK_FIXNUM(kk);
  k = FIX2INT(kk);
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  v2 = FUNCTION(gsl_vector,alloc)(k);
  FUNCTION(gsl_sort_vector,smallest)(v2->data, k, v);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), v2);
}

static VALUE FUNCTION(rb_gsl_sort_vector,largest)(VALUE obj, VALUE kk)
{
  GSL_TYPE(gsl_vector) *v = NULL, *v2 = NULL;
  size_t k;
  CHECK_FIXNUM(kk);
  k = FIX2INT(kk);
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  v2 = FUNCTION(gsl_vector,alloc)(k);
  FUNCTION(gsl_sort_vector,largest)(v2->data, k, v);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), v2);
}

static VALUE FUNCTION(rb_gsl_sort_vector,smallest_index)(VALUE obj, VALUE kk)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_index *p = NULL;
  size_t k;
  CHECK_FIXNUM(kk);
  k = FIX2INT(kk);
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  p = gsl_permutation_alloc(k);
  FUNCTION(gsl_sort_vector,smallest_index)(p->data, k, v);
  return Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, p);
}

static VALUE FUNCTION(rb_gsl_sort_vector,largest_index)(VALUE obj, VALUE kk)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_index *p = NULL;
  size_t k;
  CHECK_FIXNUM(kk);
  k = FIX2INT(kk);
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  p = gsl_permutation_alloc(k);
  FUNCTION(gsl_sort_vector,largest_index)(p->data, k, v);
  return Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, p);
}

static VALUE FUNCTION(rb_gsl_vector,histogram)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_histogram *h;
  gsl_vector *ranges;
  double min, max;
  size_t i, n;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  switch (argc) {
  case 1:
    if (rb_obj_is_kind_of(argv[0], rb_cRange)) 
      argv[0] = rb_gsl_range2ary(argv[0]);
    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      n = NUM2INT(argv[0]);
      min = FUNCTION(gsl_vector,min)(v) - 4*GSL_DBL_EPSILON;
      max = FUNCTION(gsl_vector,max)(v) + 4*GSL_DBL_EPSILON;
      h = gsl_histogram_alloc(n);
      gsl_histogram_set_ranges_uniform(h, min, max);
      break;
    case T_ARRAY:
      n = RARRAY_LEN(argv[0]) - 1;
      h = gsl_histogram_alloc(n);
      for (i = 0; i <= n; i++) h->range[i] = NUM2DBL(rb_ary_entry(argv[0], i));
      break;
    default:
      CHECK_VECTOR(argv[0]);
      Data_Get_Struct(argv[0], gsl_vector, ranges);
      n = ranges->size - 1;
      h = gsl_histogram_alloc(n);
      gsl_histogram_set_ranges(h, ranges->data, ranges->size);
      break;
    }
    break;
  case 2:
    n = NUM2INT(argv[0]);
    switch (TYPE(argv[1])) {
    case T_ARRAY:
      min = NUM2DBL(rb_ary_entry(argv[1], 0));
      max = NUM2DBL(rb_ary_entry(argv[1], 1));
      break;
    default:
      rb_raise(rb_eTypeError, "wrong argument type %s (Array expected)",
         rb_class2name(CLASS_OF(argv[1])));
      break;
    }
    h = gsl_histogram_alloc(n);
    gsl_histogram_set_ranges_uniform(h, min, max);
    break;
  case 3:
    n = NUM2INT(argv[0]);
    min = NUM2DBL(argv[1]); max = NUM2DBL(argv[2]);
    h = gsl_histogram_alloc(n);
    gsl_histogram_set_ranges_uniform(h, min, max);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments %d", argc);
    break;
  }
  for (i = 0; i < v->size; i++)
    gsl_histogram_increment(h, FUNCTION(gsl_vector,get)(v, i));
  return Data_Wrap_Struct(cgsl_histogram, 0, gsl_histogram_free, h);
}

static VALUE FUNCTION(rb_gsl_vector,last)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return C_TO_VALUE(FUNCTION(gsl_vector,get)(v, v->size-1));
}

static VALUE FUNCTION(rb_gsl_vector,first)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return C_TO_VALUE(FUNCTION(gsl_vector,get)(v, 0));
}

static VALUE FUNCTION(rb_gsl_vector,concat)(VALUE obj, VALUE other)
{
  GSL_TYPE(gsl_vector) *v = NULL, *v2 = NULL, *vnew = NULL;
  QUALIFIED_VIEW(gsl_vector,view) vv;
  VALUE x;
  BASE beg, end;
  int step;
  size_t i, size2;

  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);

  switch(TYPE(other)) {
    case T_FIXNUM:
    case T_BIGNUM:
    case T_FLOAT:
      vnew = FUNCTION(gsl_vector,alloc)(v->size + 1);
      vv = FUNCTION(gsl_vector,subvector)(vnew, 0, v->size);
      FUNCTION(gsl_vector,memcpy)(&vv.vector, v);
      FUNCTION(gsl_vector,set)(vnew, v->size, NUMCONV2(other));
      break;

    case T_ARRAY:
      size2 = RARRAY_LEN(other);
      vnew = FUNCTION(gsl_vector,alloc)(v->size + size2);
      vv = FUNCTION(gsl_vector,subvector)(vnew, 0, v->size);
      FUNCTION(gsl_vector,memcpy)(&vv.vector, v);
      for (i = 0; i < size2; i++) {
        x = rb_ary_entry(other, i);
        FUNCTION(gsl_vector,set)(vnew, v->size + i, NUMCONV2(x));
      }
      break;

    default:
      if(rb_obj_is_kind_of(other, rb_cRange)) {
        FUNCTION(get_range,beg_en_n)(other, &beg, &end, &size2, &step);
        vnew = FUNCTION(gsl_vector,alloc)(v->size + size2);
        vv = FUNCTION(gsl_vector,subvector)(vnew, 0, v->size);
        FUNCTION(gsl_vector,memcpy)(&vv.vector, v);
        for (i = 0; i < size2; i++) {
          FUNCTION(gsl_vector,set)(vnew, v->size + i, beg);
          beg += step;
        }
      } else if (rb_obj_is_kind_of(other, GSL_TYPE(cgsl_vector))) {
        Data_Get_Struct(other, GSL_TYPE(gsl_vector), v2);
        size2 = v2->size;
        vnew = FUNCTION(gsl_vector,alloc)(v->size + size2);
        vv = FUNCTION(gsl_vector,subvector)(vnew, 0, v->size);
        FUNCTION(gsl_vector,memcpy)(&vv.vector, v);
        vv = FUNCTION(gsl_vector,subvector)(vnew, v->size, size2);
        FUNCTION(gsl_vector,memcpy)(&vv.vector, v2);
      } else {
        rb_raise(rb_eTypeError, "wrong argument type %s (Array, Numeric, Range, or %s expected)",
            rb_class2name(CLASS_OF(other)), rb_class2name(GSL_TYPE(cgsl_vector)));
      }
      break;
  }

  return Data_Wrap_Struct(VEC_ROW_COL(obj), 0, FUNCTION(gsl_vector,free), vnew);  
}

void FUNCTION(mygsl_vector,diff)(GSL_TYPE(gsl_vector) *vdst,
         GSL_TYPE(gsl_vector) *vsrc, size_t n)
{
  BASE a, b;
  int coef, fac, nn, ff;
  size_t i, k;
  nn = gsl_sf_fact((unsigned int) n);
  if (GSL_IS_EVEN(n)) ff = 1;
  else ff = -1;
  for (i = 0; i < vsrc->size-n; i++) {
    fac = ff;
    a = 0;
    for (k = 0; k <= n; k++) {
      b = FUNCTION(gsl_vector,get)(vsrc, i+k);
      coef = nn/gsl_sf_fact(k)/gsl_sf_fact(n-k);
      a += fac*coef*b;
      fac *= -1;
    }
    FUNCTION(gsl_vector,set)(vdst, i, a);
  }
}

static VALUE FUNCTION(rb_gsl_vector,diff)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  size_t n;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
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
  if (v->size <= n) return obj;
  vnew = FUNCTION(gsl_vector,alloc)(v->size - n);
  FUNCTION(mygsl_vector,diff)(vnew, v, n);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), vnew);
}

static VALUE FUNCTION(rb_gsl_vector,test)(VALUE obj, int (*f)(const double))
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector_int *vi = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vi = gsl_vector_int_alloc(v->size);
  for (i = 0; i < v->size; i++) 
    gsl_vector_int_set(vi, i, (*f)(FUNCTION(gsl_vector,get)(v, i)));
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, vi);
}

static VALUE FUNCTION(rb_gsl_vector,test2)(VALUE obj, int (*f)(const double))
{
  GSL_TYPE(gsl_vector) *v = NULL;
  VALUE ary;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  ary = rb_ary_new2(v->size);
  for (i = 0; i < v->size; i++) {
    if ((*f)(FUNCTION(gsl_vector,get)(v, i)))
      rb_ary_store(ary, i, Qtrue);
    else
      rb_ary_store(ary, i, Qfalse);
  }
  return ary;
}

static VALUE FUNCTION(rb_gsl_vector,isnan)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,test)(obj, gsl_isnan);
}

static VALUE FUNCTION(rb_gsl_vector,isinf)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,test)(obj, gsl_isinf);
}

static VALUE FUNCTION(rb_gsl_vector,finite)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,test)(obj, gsl_finite);
}

static VALUE FUNCTION(rb_gsl_vector,isnan2)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,test2)(obj, gsl_isnan);
}

static VALUE FUNCTION(rb_gsl_vector,isinf2)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,test2)(obj, gsl_isinf);
}

static VALUE FUNCTION(rb_gsl_vector,finite2)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,test2)(obj, gsl_finite);
}

static VALUE FUNCTION(rb_gsl_vector,delete_at)(VALUE obj, VALUE ii)
{
  int i2;
  size_t i;
  GSL_TYPE(gsl_vector) *v;
  BASE x;
  if (rb_obj_is_kind_of(obj,QUALIFIED_VIEW(cgsl_vector,view)))
    rb_raise(rb_eRuntimeError, "prohibited for %s", rb_class2name(CLASS_OF(obj)));
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (v->stride != 1) rb_raise(rb_eRuntimeError, "vector must have stride 1");
  if (v->size == 0) return Qnil;
  CHECK_FIXNUM(ii);
  i2 = FIX2INT(ii);
  if (i2 < 0) {
    i2 += v->size;
  }
  if (i2 < 0 || i2 > (int) (v->size-1)) return Qnil;
  i = (size_t) i2;
  x = FUNCTION(gsl_vector,get)(v, i);
  memmove(v->data+i, v->data+i+1, sizeof(BASE)*(v->size-1-i));
  v->size -= 1;
  return C_TO_VALUE(x);
}

static VALUE FUNCTION(rb_gsl_vector,delete_if)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v;
  BASE x;
  VALUE val;
  size_t i, count = 0;
  if (!rb_block_given_p()) rb_raise(rb_eRuntimeError, "block is not given");
  if (rb_obj_is_kind_of(obj,QUALIFIED_VIEW(cgsl_vector,view)))
    rb_raise(rb_eRuntimeError, "prohibited for %s", rb_class2name(CLASS_OF(obj)));
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (v->stride != 1) rb_raise(rb_eRuntimeError, "vector must have stride 1");
  for (i = 0; i < v->size; i++) {
    x = FUNCTION(gsl_vector,get)(v, i);
    val = rb_yield(C_TO_VALUE(x));
    if(RTEST(val)) {
      count++;
    } else if(count > 0) {
      FUNCTION(gsl_vector,set)(v, i-count, x);
    }
  }
  v->size -= count;
  return obj;
}

static VALUE FUNCTION(rb_gsl_vector,delete)(VALUE obj, VALUE yy)
{
  GSL_TYPE(gsl_vector) *v;
  BASE x, y;
  size_t i, count = 0;
  y = NUMCONV(yy);
  if (rb_obj_is_kind_of(obj,QUALIFIED_VIEW(cgsl_vector,view)))
    rb_raise(rb_eRuntimeError, "prohibited for %s", rb_class2name(CLASS_OF(obj)));
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (v->stride != 1) rb_raise(rb_eRuntimeError, "vector must have stride 1");
  if (v->size == 0) return obj;
  for (i = 0; i < v->size; i++) {
    x = FUNCTION(gsl_vector,get)(v, i);
    if (x == y) {
      count++;
    } else if(count > 0) {
      FUNCTION(gsl_vector,set)(v, i-count, x);
    }
  }
  v->size -= count;
  return count ? (VALUE) y : Qnil;
}

/* singleton method */
#ifdef BASE_INT
#define FORMAT_TMP "%d"
#else
#define FORMAT_TMP "%lf"
#endif
static VALUE FUNCTION(rb_gsl_vector,filescan)(VALUE klass, VALUE file)
{
  FILE *fp = NULL;
  int nn, k;
  char buf[1024], filename[1024];
  size_t n, lines, i, j, ii = 0, jj;
  long pos;
  GSL_TYPE(gsl_vector) **x;
  BASE val;
  VALUE ary;
  Check_Type(file, T_STRING);
  strcpy(filename, STR2CSTR(file));
  sprintf(buf, "sed '/^#/d' %s | wc", filename);
  if ((fp = popen(buf, "r")) == NULL) 
    rb_raise(rb_eIOError, "popen failed.");
  if (fgets(buf, 1024, fp) == NULL)
    rb_sys_fail(0);
  pclose(fp);
  sscanf(buf, "%d", &nn);
  lines = (size_t) nn;  /*  vector length */
  if ((fp = fopen(filename, "r")) == NULL) 
    rb_raise(rb_eIOError, "cannot open file %s.", filename);
  while (1) {
    if (fgets(buf, 1024, fp) == NULL)    /* read the first line to count number of columns */
      rb_sys_fail(0);
    if (buf[0] == '#') continue;
    else break;
  }
  n = count_columns(buf);  /* number of vectors to be created */
  x = (GSL_TYPE(gsl_vector)**) xmalloc(sizeof(GSL_TYPE(gsl_vector)*)*n);
  ary = rb_ary_new2(n);
  for (j = 0; j < n; j++) {
    x[j] = FUNCTION(gsl_vector,alloc)(lines);
    rb_ary_store(ary, j, Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), x[j]));
  }
  rewind(fp);
  for (i = 0, ii = 0; ii < lines; i++) {
    pos = ftell(fp);
    if (fgets(buf, 1024, fp) == NULL)
      rb_sys_fail(0);
    if (buf[0] == '#') continue;
    fseek(fp, pos, SEEK_SET);
    for (j = 0, jj = 0; jj < n; j++) {
      k = fscanf(fp, FORMAT_TMP, &val);
      if (k != 1) continue;
      FUNCTION(gsl_vector,set)(x[jj++], ii, val);
    }
    ii += 1;
  }
  fclose(fp);
  free(x);
  return ary;
}
#undef FORMAT_TMP

static int FUNCTION(gsl_vector,eq)(const GSL_TYPE(gsl_vector) *a, 
           const GSL_TYPE(gsl_vector) *b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b->data[i*b->stride];
    c->data[i] = (x > y || x < y) ? 0 : 1;
  }
  return 0;
}

static int FUNCTION(gsl_vector,ne)(const GSL_TYPE(gsl_vector) *a, 
           const GSL_TYPE(gsl_vector) *b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b->data[i*b->stride];
    c->data[i] = (x > y || x < y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,gt)(const GSL_TYPE(gsl_vector) *a, 
           const GSL_TYPE(gsl_vector) *b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b->data[i*b->stride];
    c->data[i] =  (x > y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,ge)(const GSL_TYPE(gsl_vector) *a, 
           const GSL_TYPE(gsl_vector) *b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b->data[i*b->stride];
    c->data[i] =  (x >= y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,lt)(const GSL_TYPE(gsl_vector) *a, 
           const GSL_TYPE(gsl_vector) *b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b->data[i*b->stride];
    c->data[i] =  (x < y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,le)(const GSL_TYPE(gsl_vector) *a, 
           const GSL_TYPE(gsl_vector) *b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b->data[i*b->stride];
    c->data[i] =  (x <= y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,and)(const GSL_TYPE(gsl_vector) *a, 
           const GSL_TYPE(gsl_vector) *b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b->data[i*b->stride];
    c->data[i] = (x != 0 && y != 0) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,or)(const GSL_TYPE(gsl_vector) *a, 
           const GSL_TYPE(gsl_vector) *b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b->data[i*b->stride];
    c->data[i] = (x != 0 || y != 0) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,xor)(const GSL_TYPE(gsl_vector) *a, 
            const GSL_TYPE(gsl_vector) *b,
            gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b->data[i*b->stride];
    c->data[i] = ((x != 0) == (y != 0)) ? 0 : 1;
  }
  return 0;
}

static int FUNCTION(gsl_vector,eq2)(const GSL_TYPE(gsl_vector) *a, 
            BASE b, gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b;
    c->data[i] = (x > y || x < y) ? 0 : 1;
  }
  return 0;
}

static int FUNCTION(gsl_vector,ne2)(const GSL_TYPE(gsl_vector) *a, 
           BASE b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b;
    c->data[i] = (x > y || x < y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,gt2)(const GSL_TYPE(gsl_vector) *a, 
           BASE b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b;
    c->data[i] =  (x > y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,ge2)(const GSL_TYPE(gsl_vector) *a, 
           BASE b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b;
    c->data[i] =  (x >= y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,lt2)(const GSL_TYPE(gsl_vector) *a, 
           BASE b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b;
    c->data[i] =  (x < y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,le2)(const GSL_TYPE(gsl_vector) *a, 
            BASE b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b;
    c->data[i] =  (x <= y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,and2)(const GSL_TYPE(gsl_vector) *a, 
            BASE b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b;
    c->data[i] = (x != 0 && y != 0) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,or2)(const GSL_TYPE(gsl_vector) *a, 
            BASE b,
           gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b;
    c->data[i] = (x != 0 || y != 0) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_vector,xor2)(const GSL_TYPE(gsl_vector) *a, 
             BASE b,
            gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i*a->stride];
    y = b;
    c->data[i] = ((x != 0) == (y != 0)) ? 0 : 1;
  }
  return 0;
}


static VALUE FUNCTION(rb_gsl_vector,compare)(VALUE aa, VALUE bb,
               int (*cmp)(const GSL_TYPE(gsl_vector)*,
              const GSL_TYPE(gsl_vector)*,
              gsl_block_uchar*),
               int (*cmp2)(const GSL_TYPE(gsl_vector)*,
               BASE,
               gsl_block_uchar*))
{
  GSL_TYPE(gsl_vector) *a, *b;
  /*  gsl_vector_int *c;*/
  gsl_block_uchar *c;
  //int status;
  Data_Get_Struct(aa, GSL_TYPE(gsl_vector), a);
  c = gsl_block_uchar_alloc(a->size);
  if (VEC_P(bb)) {
    Data_Get_Struct(bb, GSL_TYPE(gsl_vector), b);
    if (a->size != b->size) 
      rb_raise(rb_eRuntimeError, "Vector size mismatch, %d and %d", (int) a->size, 
         (int) b->size);
    /*status =*/ (*cmp)(a, b, c);
  } else {
    /*status =*/ (*cmp2)(a, NUMCONV(bb), c);
  }
  return Data_Wrap_Struct(cgsl_block_uchar, 0, gsl_block_uchar_free, c);
}

static VALUE FUNCTION(rb_gsl_vector,eq)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_vector,compare)(aa, bb, FUNCTION(gsl_vector,eq),
           FUNCTION(gsl_vector,eq2));
}

static VALUE FUNCTION(rb_gsl_vector,ne)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_vector,compare)(aa, bb, FUNCTION(gsl_vector,ne),
           FUNCTION(gsl_vector,ne2));
}

static VALUE FUNCTION(rb_gsl_vector,gt)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_vector,compare)(aa, bb, FUNCTION(gsl_vector,gt),
           FUNCTION(gsl_vector,gt2));
}

static VALUE FUNCTION(rb_gsl_vector,ge)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_vector,compare)(aa, bb, FUNCTION(gsl_vector,ge),
           FUNCTION(gsl_vector,ge2));
}

static VALUE FUNCTION(rb_gsl_vector,lt)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_vector,compare)(aa, bb, FUNCTION(gsl_vector,lt),
           FUNCTION(gsl_vector,lt2));
}

static VALUE FUNCTION(rb_gsl_vector,le)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_vector,compare)(aa, bb, FUNCTION(gsl_vector,le),
           FUNCTION(gsl_vector,le2));
}

static VALUE FUNCTION(rb_gsl_vector,and)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_vector,compare)(aa, bb, FUNCTION(gsl_vector,and),
           FUNCTION(gsl_vector,and2));
}

static VALUE FUNCTION(rb_gsl_vector,or)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_vector,compare)(aa, bb, FUNCTION(gsl_vector,or),
           FUNCTION(gsl_vector,or2));
}

static VALUE FUNCTION(rb_gsl_vector,xor)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_vector,compare)(aa, bb, FUNCTION(gsl_vector,xor),
           FUNCTION(gsl_vector,xor2));
}

static VALUE FUNCTION(rb_gsl_vector,not)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v;
  gsl_block_uchar *vv;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vv = gsl_block_uchar_alloc(v->size);
  for (i = 0; i < v->size; i++) vv->data[i] = (v->data[i*v->stride] != 0) ? 0 : 1;
  return Data_Wrap_Struct(cgsl_block_uchar, 0, gsl_block_uchar_free, vv);
}

static VALUE FUNCTION(rb_gsl_vector,any)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t i;

  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (rb_block_given_p()) {
    for (i = 0; i < v->size; i++) {
      if (rb_yield(C_TO_VALUE(FUNCTION(gsl_vector,get)(v, i)))) return INT2FIX(1);
    }
    return INT2FIX(0);
  } else {
    if (FUNCTION(gsl_vector,isnull)(v)) return INT2FIX(0);
    return INT2FIX(1);
  }
}

static VALUE FUNCTION(rb_gsl_vector,any2)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  size_t i;

  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (rb_block_given_p()) {
    for (i = 0; i < v->size; i++) 
      if (rb_yield(C_TO_VALUE(FUNCTION(gsl_vector,get)(v, i)))) return Qtrue;
    return Qfalse;
  } else {
    for (i = 0; i < v->size; i++) 
      if (v->data[i*v->stride]) return Qtrue;
    return Qfalse;
  }
}

static VALUE FUNCTION(rb_gsl_vector,none)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v;
  size_t i;

  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (rb_block_given_p()) {
    for (i = 0; i < v->size; i++) 
      if (rb_yield(C_TO_VALUE(FUNCTION(gsl_vector,get)(v, i)))) return Qfalse;
    return Qtrue;
  } else {
    for (i = 0; i < v->size; i++) 
      if (v->data[i*v->stride]) return Qfalse;
    return Qtrue;
  }
}

static VALUE FUNCTION(rb_gsl_vector,all)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v;
  size_t i;

  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (rb_block_given_p()) {
    for (i = 0; i < v->size; i++) 
      if (!rb_yield(C_TO_VALUE(FUNCTION(gsl_vector,get)(v, i)))) return Qfalse;
    return Qtrue;
  } else {
    for (i = 0; i < v->size; i++) 
      if (!v->data[i*v->stride]) return Qfalse;
    return Qtrue;
  }
}

static VALUE FUNCTION(rb_gsl_vector,where)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v;
  gsl_index *vv;
  gsl_block_uchar *btmp = NULL;
  size_t i, j, n = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  /* count true elements */
  if (rb_block_given_p()) {
    btmp = gsl_block_uchar_alloc(v->size);
    for (i = 0; i < v->size; i++) {
      if (rb_yield(C_TO_VALUE(FUNCTION(gsl_vector,get)(v, i)))) {
  btmp->data[i] = 1;
  n++;
      } else {
  btmp->data[i] = 0;
      }
    }  /* for */
  } else { /* block is not given */
    for (i = 0; i < v->size; i++) { if (FUNCTION(gsl_vector,get)(v, i)) n++; }
  }
  if (n == 0) {
    if (btmp) gsl_block_uchar_free(btmp);
    return Qnil;
  }
  vv = gsl_permutation_alloc(n);
  for (i = 0, j = 0; i < v->size; i++) {
    if ((!btmp && FUNCTION(gsl_vector,get)(v, i)) || (btmp && btmp->data[i])) {
      vv->data[j++] = i;
    }
  }
  if (btmp) gsl_block_uchar_free(btmp);
  return Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, vv);
}

static VALUE FUNCTION(rb_gsl_vector,where2)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v;
  gsl_index *v1, *v2;
  gsl_block_uchar *btmp = NULL;
  VALUE vv1, vv2;
  size_t i, j, k, n = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (rb_block_given_p()) {
    btmp = gsl_block_uchar_alloc(v->size);
    for (i = 0; i < v->size; i++) {
      if (rb_yield(C_TO_VALUE(FUNCTION(gsl_vector,get)(v, i)))) {
  btmp->data[i] = 1;
  n++;
      } else {
  btmp->data[i] = 0;
      }
    } /* for */
  } else {  /* block is not given */
    for (i = 0; i < v->size; i++) { if (FUNCTION(gsl_vector,get)(v, i)) n++; }
  }
  /* true and false logic.  need to handle both */
  if (n == 0) {
    v2 = gsl_permutation_calloc(v->size);  /* calloc() initializes v2 */
    vv1 = Qnil;
    vv2 = Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, v2);
  } else if (v->size-n == 0) { 
    v1 = gsl_permutation_calloc(n);           /* calloc() initializes v1 */
    vv1 = Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, v1);     
    vv2 = Qnil;
  } else { 
    /* same case as 'where' */
    v1 = gsl_permutation_alloc(n);
    v2 = gsl_permutation_alloc(v->size-n);
    for (i = 0, j = 0, k = 0; i < v->size; i++) {
      if ((!btmp && FUNCTION(gsl_vector,get)(v, i)) || (btmp && btmp->data[i])) v1->data[j++] = i;
      else v2->data[k++] = i;
    }
    vv1 = Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, v1);
    vv2 = Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, v2);
  }
  if (btmp) gsl_block_uchar_free(btmp);
  return rb_ary_new3(2, vv1, vv2);
}

static VALUE FUNCTION(rb_gsl_vector,op_inplace)(VALUE vv1, VALUE vv2, 
            int (*f)(GSL_TYPE(gsl_vector)*, const GSL_TYPE(gsl_vector)*))
{
  GSL_TYPE(gsl_vector) *v1, *v2;
  Data_Get_Struct(vv1, GSL_TYPE(gsl_vector), v1);
  Data_Get_Struct(vv2, GSL_TYPE(gsl_vector), v2);
  (*f)(v1, v2);
  return vv1;
}

static VALUE FUNCTION(rb_gsl_vector,add_inplace)(VALUE vv1, VALUE vv2)
{
  GSL_TYPE(gsl_vector) *v1;
  double x;
  if(VEC_P(vv2)) {
    return FUNCTION(rb_gsl_vector,op_inplace)(vv1, vv2, FUNCTION(gsl_vector,add));
  } else {
    x = NUM2DBL(vv2);
    Data_Get_Struct(vv1, GSL_TYPE(gsl_vector), v1);
    FUNCTION(gsl_vector,add_constant)(v1, x);
    return vv1;
  }
}

static VALUE FUNCTION(rb_gsl_vector,sub_inplace)(VALUE vv1, VALUE vv2)
{
  GSL_TYPE(gsl_vector) *v1;
  double x;
  if(VEC_P(vv2)) {
    return FUNCTION(rb_gsl_vector,op_inplace)(vv1, vv2, FUNCTION(gsl_vector,sub));
  } else {
    x = NUM2DBL(vv2);
    Data_Get_Struct(vv1, GSL_TYPE(gsl_vector), v1);
    FUNCTION(gsl_vector,add_constant)(v1, -x);
    return vv1;
  }
}

static VALUE FUNCTION(rb_gsl_vector,mul_inplace)(VALUE vv1, VALUE vv2)
{
  GSL_TYPE(gsl_vector) *v1;
  double x;
  if(VEC_P(vv2)) {
    return FUNCTION(rb_gsl_vector,op_inplace)(vv1, vv2, FUNCTION(gsl_vector,mul));
  } else {
    x = NUM2DBL(vv2);
    Data_Get_Struct(vv1, GSL_TYPE(gsl_vector), v1);
    FUNCTION(gsl_vector,scale)(v1, x);
    return vv1;
  }
}

static VALUE FUNCTION(rb_gsl_vector,div_inplace)(VALUE vv1, VALUE vv2)
{
  GSL_TYPE(gsl_vector) *v1;
  double x;
  if(VEC_P(vv2)) {
    return FUNCTION(rb_gsl_vector,op_inplace)(vv1, vv2, FUNCTION(gsl_vector,div));
  } else {
    x = NUM2DBL(vv2);
    Data_Get_Struct(vv1, GSL_TYPE(gsl_vector), v1);
    FUNCTION(gsl_vector,scale)(v1, 1.0/x);
    return vv1;
  }
}

static VALUE FUNCTION(rb_gsl_vector,zip)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v0, **vp, *vnew;
  VALUE ary;
  size_t i, j;
  int argc2;
  VALUE *argv2;
  if (VEC_P(obj)) {
    Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v0);
    argc2 = argc;
    argv2 = argv;
  } else {
    if (argc < 1) rb_raise(rb_eArgError, "Too few arguments.");
    Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), v0);    
    argc2 = argc - 1;
    argv2 = argv + 1;
  }
  for (i = 0; (int) i < argc2; i++) {
    CHECK_VEC(argv2[i]);
  }
  vp = (GSL_TYPE(gsl_vector)**) malloc(sizeof(GSL_TYPE(gsl_vector)**));
  for (i = 0; (int) i < argc2; i++) {
    Data_Get_Struct(argv2[i], GSL_TYPE(gsl_vector), vp[i]);
  }
  ary = rb_ary_new2(v0->size);
  for (i = 0; i < v0->size; i++) {
    vnew = FUNCTION(gsl_vector,alloc)(argc2 + 1);
    FUNCTION(gsl_vector,set)(vnew, 0, FUNCTION(gsl_vector,get)(v0, i));
    for (j = 0; (int) j < argc2; j++) {
      if (i < vp[j]->size) {
  FUNCTION(gsl_vector,set)(vnew, j+1, FUNCTION(gsl_vector,get)(vp[j], i));
      } else {
  FUNCTION(gsl_vector,set)(vnew, j+1, 0.0);
      }
    }
    rb_ary_store(ary, i, Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), vnew));
  }
  
  free((GSL_TYPE(gsl_vector)**) vp);
  return ary;
}

static VALUE FUNCTION(rb_gsl_vector,join)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v;
  VALUE str, sep;
  char *p, buf[16];
  size_t i;
  switch (argc) {
  case 0:
    sep = rb_str_new2(" ");
    break;
  case 1:
    sep = argv[0];
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0 or 1)", argc);
  }
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  //  p = (char *) malloc((10+RSTRING(sep))*v->size + 1);
  p = (char *) malloc((10+RSTRING_LEN(sep))*v->size + 1);
  str = rb_str_new2(p);
  for (i = 0; i < v->size; i++) {
#ifdef BASE_DOUBLE
    sprintf(buf, "%4.3e", FUNCTION(gsl_vector,get)(v, i));
#else
    sprintf(buf, "%d", FUNCTION(gsl_vector,get)(v, i));
#endif
    rb_str_concat(str, rb_str_new2(buf));
    if (i != v->size-1) rb_str_concat(str, sep);
  }
  return str;
}

static VALUE FUNCTION(rb_gsl_vector,cumsum)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v, *vnew;
  BASE sum = 0;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  for (i = 0; i < v->size; i++) {
    sum += FUNCTION(gsl_vector,get)(v, i);
    FUNCTION(gsl_vector,set)(vnew, i, sum);
  }
  return Data_Wrap_Struct(VEC_ROW_COL(obj), 0, FUNCTION(gsl_vector,free), vnew);
}

static VALUE FUNCTION(rb_gsl_vector,cumprod)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v, *vnew;
  BASE prod = 1;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_vector,alloc)(v->size);
  for (i = 0; i < v->size; i++) {
    prod *= FUNCTION(gsl_vector,get)(v, i);
    FUNCTION(gsl_vector,set)(vnew, i, prod);
  }
  return Data_Wrap_Struct(VEC_ROW_COL(obj), 0, FUNCTION(gsl_vector,free), vnew);
}

#ifdef GSL_1_9_LATER
static VALUE FUNCTION(rb_gsl_vector,property)(VALUE obj,
  int (*f)(const GSL_TYPE(gsl_vector) *)) {
  GSL_TYPE(gsl_vector) *v;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return INT2FIX((*f)(v));  
}

static VALUE FUNCTION(rb_gsl_vector,property2)(VALUE obj,
  int (*f)(const GSL_TYPE(gsl_vector) *)) {
  GSL_TYPE(gsl_vector) *v;  
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if ((*f)(v)) return Qtrue;
  else return Qfalse;
}

static VALUE FUNCTION(rb_gsl_vector,ispos)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,property)(obj, FUNCTION(gsl_vector, ispos));
}
static VALUE FUNCTION(rb_gsl_vector,ispos2)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,property2)(obj, FUNCTION(gsl_vector, ispos));  
}
static VALUE FUNCTION(rb_gsl_vector,isneg)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,property)(obj, FUNCTION(gsl_vector, isneg));    
}
static VALUE FUNCTION(rb_gsl_vector,isneg2)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,property2)(obj, FUNCTION(gsl_vector, isneg));      
}
#endif

#ifdef GSL_1_10_LATER
static VALUE FUNCTION(rb_gsl_vector,isnonneg)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,property)(obj, FUNCTION(gsl_vector, isnonneg));
}
static VALUE FUNCTION(rb_gsl_vector,isnonneg2)(VALUE obj)
{
  return FUNCTION(rb_gsl_vector,property2)(obj, FUNCTION(gsl_vector, isnonneg));
}  
#endif

void FUNCTION(Init_gsl_vector,init)(VALUE module)
{
  /*  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "new", 
      FUNCTION(rb_gsl_vector,new), -1);*/
  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "[]", 
           FUNCTION(rb_gsl_vector,new), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "alloc", 
           FUNCTION(rb_gsl_vector,new), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "calloc", 
           FUNCTION(rb_gsl_vector,calloc), 1);

/*****/
  rb_define_method(GSL_TYPE(cgsl_vector), "get", FUNCTION(rb_gsl_vector,get), -1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "[]", "get");
  rb_define_method(GSL_TYPE(cgsl_vector), "size", FUNCTION(rb_gsl_vector,size), 0);
  rb_define_alias(GSL_TYPE(cgsl_vector), "len", "size");
  rb_define_alias(GSL_TYPE(cgsl_vector), "length", "size");
  rb_define_method(GSL_TYPE(cgsl_vector), "stride", FUNCTION(rb_gsl_vector,stride), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "set_stride", FUNCTION(rb_gsl_vector,set_stride), 1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "stride=", "set_stride");
  rb_define_method(GSL_TYPE(cgsl_vector), "owner", FUNCTION(rb_gsl_vector,owner), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "set", FUNCTION(rb_gsl_vector,set), -1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "[]=", "set");
  rb_define_method(GSL_TYPE(cgsl_vector), "set_all", FUNCTION(rb_gsl_vector,set_all), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "set_zero", FUNCTION(rb_gsl_vector,set_zero), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "set_basis", FUNCTION(rb_gsl_vector,set_basis), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "each", FUNCTION(rb_gsl_vector,each), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "reverse_each", FUNCTION(rb_gsl_vector,reverse_each), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "each_index", FUNCTION(rb_gsl_vector,each_index), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "reverse_each_index", FUNCTION(rb_gsl_vector,reverse_each_index), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "to_a", FUNCTION(rb_gsl_vector,to_a), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "reverse", FUNCTION(rb_gsl_vector,reverse), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "reverse!", FUNCTION(rb_gsl_vector,reverse_bang), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "max", FUNCTION(rb_gsl_vector,max), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "min", FUNCTION(rb_gsl_vector,min), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "minmax", FUNCTION(rb_gsl_vector,minmax), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "maxmin", FUNCTION(rb_gsl_vector,maxmin), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "max_index", FUNCTION(rb_gsl_vector,max_index), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "min_index", FUNCTION(rb_gsl_vector,min_index), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "minmax_index", FUNCTION(rb_gsl_vector,minmax_index), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "maxmin_index", FUNCTION(rb_gsl_vector,maxmin_index), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "isnull", FUNCTION(rb_gsl_vector,isnull), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "isnull?", FUNCTION(rb_gsl_vector,isnull2), 0);
  /*  rb_define_alias(GSL_TYPE(cgsl_vector), "none?", "isnull?");*/
  /* none? method is define below, which can have a block. */

  rb_define_method(GSL_TYPE(cgsl_vector), "trans", FUNCTION(rb_gsl_vector,trans), 0);
  rb_define_alias(GSL_TYPE(cgsl_vector), "transpose", "trans");
  rb_define_alias(GSL_TYPE(cgsl_vector), "col", "trans");

  rb_define_method(GSL_TYPE(cgsl_vector), "trans!", FUNCTION(rb_gsl_vector,trans_bang), 0);
  rb_define_alias(GSL_TYPE(cgsl_vector), "transpose!", "trans!");
  rb_define_alias(GSL_TYPE(cgsl_vector), "col!", "trans!");
#ifdef BASE_DOUBLE
  rb_define_alias(cgsl_vector_col, "row", "trans");
  rb_define_alias(cgsl_vector_col, "row!", "trans!");
#elif defined(BASE_INT)
  rb_define_alias(cgsl_vector_int_col, "row", "trans");
  rb_define_alias(cgsl_vector_int_col, "row!", "trans!");
#endif

  rb_define_method(GSL_TYPE(cgsl_vector), "-@", FUNCTION(rb_gsl_vector,uminus), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "+@", FUNCTION(rb_gsl_vector,uplus), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "sum", FUNCTION(rb_gsl_vector,sum), 0);
#ifdef BASE_INT
  /* Vector#sumsq is defined in blas1.c */
  rb_define_method(GSL_TYPE(cgsl_vector), "sumsq", FUNCTION(rb_gsl_vector,sumsq), 0);
#endif
  rb_define_method(GSL_TYPE(cgsl_vector), "prod", FUNCTION(rb_gsl_vector,prod), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "cumsum", FUNCTION(rb_gsl_vector,cumsum), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "cumprod", FUNCTION(rb_gsl_vector,cumprod), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "connect", 
       FUNCTION(rb_gsl_vector,connect), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "connect", 
           FUNCTION(rb_gsl_vector,connect), -1);


  rb_define_method(GSL_TYPE(cgsl_vector), "sgn", FUNCTION(rb_gsl_vector,sgn), 0);
  rb_define_alias(GSL_TYPE(cgsl_vector), "signum", "sgn");
  rb_define_method(GSL_TYPE(cgsl_vector), "abs", FUNCTION(rb_gsl_vector,abs), 0);
  rb_define_alias(GSL_TYPE(cgsl_vector), "fabs", "abs");
  rb_define_method(GSL_TYPE(cgsl_vector), "square", 
       FUNCTION(rb_gsl_vector,square), 0);
  rb_define_alias(GSL_TYPE(cgsl_vector), "abs2", "square");
  rb_define_method(GSL_TYPE(cgsl_vector), "sqrt", FUNCTION(rb_gsl_vector,sqrt), 0);

  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "memcpy", 
           FUNCTION(rb_gsl_vector,memcpy), 2);
  rb_define_method(GSL_TYPE(cgsl_vector), "clone", 
       FUNCTION(rb_gsl_vector,clone), 0);
  rb_define_alias(GSL_TYPE(cgsl_vector), "duplicate", "clone");
  rb_define_alias(GSL_TYPE(cgsl_vector), "dup", "clone");
  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "swap", 
           FUNCTION(rb_gsl_vector,swap), 2);
  rb_define_method(GSL_TYPE(cgsl_vector), "swap_elements", 
       FUNCTION(rb_gsl_vector,swap_elements), 2);

  rb_define_method(GSL_TYPE(cgsl_vector), "fwrite", 
       FUNCTION(rb_gsl_vector,fwrite), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "fread", 
       FUNCTION(rb_gsl_vector,fread), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "fprintf", 
       FUNCTION(rb_gsl_vector,fprintf), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "printf", 
       FUNCTION(rb_gsl_vector,printf), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "fscanf", 
       FUNCTION(rb_gsl_vector,fscanf), 1);

  /* 2.Aug.2004 */
  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "inner_product", 
           FUNCTION(rb_gsl_vector,inner_product), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "dot", 
           FUNCTION(rb_gsl_vector,inner_product), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "inner_product",
       FUNCTION(rb_gsl_vector,inner_product), -1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "dot", "inner_product");

  rb_define_method(GSL_TYPE(cgsl_vector), "equal?", 
       FUNCTION(rb_gsl_vector,equal), -1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "==", "equal?");

  rb_define_method(GSL_TYPE(cgsl_vector), "to_poly", 
       FUNCTION(rb_gsl_vector,to_poly), 0);

/*****/
  rb_define_method(GSL_TYPE(cgsl_vector), "graph", 
       FUNCTION(rb_gsl_vector,graph), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "graph_step", 
       FUNCTION(rb_gsl_vector,graph_step), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "plot", 
       FUNCTION(rb_gsl_vector,plot), -1);

  rb_define_method(GSL_TYPE(cgsl_vector), "print",
       FUNCTION(rb_gsl_vector,print), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "inspect", 
       FUNCTION(rb_gsl_vector,inspect), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "to_s", 
       FUNCTION(rb_gsl_vector,to_s), 0);

/*****/
  rb_define_method(GSL_TYPE(cgsl_vector), "subvector",
       FUNCTION(rb_gsl_vector,subvector), -1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "view", "subvector");
  rb_define_method(GSL_TYPE(cgsl_vector), "subvector_with_stride", 
       FUNCTION(rb_gsl_vector,subvector_with_stride), -1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "view_with_stride", "subvector_with_stride");

  rb_define_method(GSL_TYPE(cgsl_vector), "matrix_view", 
       FUNCTION(rb_gsl_vector,matrix_view), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "matrix_view_with_tda", 
       FUNCTION(rb_gsl_vector,matrix_view_with_tda), 3);

#ifdef BASE_DOUBLE
  rb_undef_method(cgsl_vector_view_ro, "set");
  rb_undef_method(cgsl_vector_view_ro, "[]=");
  rb_undef_method(cgsl_vector_col_view_ro, "set");
  rb_undef_method(cgsl_vector_col_view_ro, "[]=");
#elif defined(BASE_INT)
  rb_undef_method(cgsl_vector_int_view_ro, "set");
  rb_undef_method(cgsl_vector_int_view_ro, "[]=");
  rb_undef_method(cgsl_vector_int_col_view_ro, "set");
  rb_undef_method(cgsl_vector_int_col_view_ro, "[]=");
#endif

  rb_define_method(GSL_TYPE(cgsl_vector), "scale", 
       FUNCTION(rb_gsl_vector,scale), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "scale!", 
       FUNCTION(rb_gsl_vector,scale_bang), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "add_constant", 
       FUNCTION(rb_gsl_vector,add_constant), 1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "add_const", "add_constant");
  rb_define_method(GSL_TYPE(cgsl_vector), "add_constant!", 
       FUNCTION(rb_gsl_vector,add_constant_bang), 1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "add_const!", "add_constant!");

#ifdef HAVE_TENSOR_TENSOR_H
  rb_define_method(GSL_TYPE(cgsl_vector), "to_tensor", 
       FUNCTION(rb_gsl_vector,to_tensor), -1);
#endif

  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "to_gplot", 
           FUNCTION(rb_gsl_vector,to_gplot), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "to_gsplot", 
           FUNCTION(rb_gsl_vector,to_gplot), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "to_gplot", 
       FUNCTION(rb_gsl_vector,to_gplot), -1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "to_gsplot", "to_gplot");

  /*****/
  rb_define_method(GSL_TYPE(cgsl_vector), "collect", 
       FUNCTION(rb_gsl_vector,collect), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "collect!", 
       FUNCTION(rb_gsl_vector,collect_bang), 0);
  rb_define_alias(GSL_TYPE(cgsl_vector), "map", "collect");
  rb_define_alias(GSL_TYPE(cgsl_vector), "map!", "collect!");

  rb_define_method(GSL_TYPE(cgsl_vector), "to_m", FUNCTION(rb_gsl_vector,to_m), 2);
  rb_define_alias(GSL_TYPE(cgsl_vector), "to_matrix", "to_m");
  rb_define_alias(GSL_TYPE(cgsl_vector), "reshape", "to_m");
  rb_define_method(GSL_TYPE(cgsl_vector), "to_m_diagonal", 
       FUNCTION(rb_gsl_vector,to_m_diagonal), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "block", 
       FUNCTION(rb_gsl_vector,block), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "to_m_circulant", 
       FUNCTION(rb_gsl_vector,to_m_circulant), 0);
  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "indgen", 
           FUNCTION(rb_gsl_vector,indgen_singleton), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "indgen", 
       FUNCTION(rb_gsl_vector,indgen), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "indgen!", 
       FUNCTION(rb_gsl_vector,indgen_bang), -1);
  /*****/
  rb_define_method(GSL_TYPE(cgsl_vector), "sort!", 
       GSL_TYPE(rb_gsl_sort_vector), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "sort", 
       GSL_TYPE(rb_gsl_sort_vector2), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "sort_index", 
       FUNCTION(rb_gsl_sort_vector,index), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "sort_smallest",
       FUNCTION(rb_gsl_sort_vector,smallest), 1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "smallest", "sort_smallest");
  rb_define_method(GSL_TYPE(cgsl_vector), "sort_largest",
       FUNCTION(rb_gsl_sort_vector,largest), 1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "largest", "sort_largest");
  rb_define_method(GSL_TYPE(cgsl_vector), "sort_smallest_index", 
       FUNCTION(rb_gsl_sort_vector,smallest_index), 1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "smallest_index", "sort_smallest_index");
  rb_define_method(GSL_TYPE(cgsl_vector), "sort_largest_index", 
       FUNCTION(rb_gsl_sort_vector,largest_index), 1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "largest_index", "sort_largest_index");

  /*****/
  rb_define_method(GSL_TYPE(cgsl_vector), "histogram",
       FUNCTION(rb_gsl_vector,histogram), -1);

  rb_define_method(GSL_TYPE(cgsl_vector), "last", 
       FUNCTION(rb_gsl_vector,last), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "first", 
       FUNCTION(rb_gsl_vector,first), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "concat", 
       FUNCTION(rb_gsl_vector,concat), 1);

  rb_define_method(GSL_TYPE(cgsl_vector), "diff", 
       FUNCTION(rb_gsl_vector,diff), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "isnan", 
       FUNCTION(rb_gsl_vector,isnan), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "isinf", 
       FUNCTION(rb_gsl_vector,isinf), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "finite", 
       FUNCTION(rb_gsl_vector,finite), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "isnan?", 
       FUNCTION(rb_gsl_vector,isnan2), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "isinf?", 
       FUNCTION(rb_gsl_vector,isinf2), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "finite?", 
       FUNCTION(rb_gsl_vector,finite2), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "delete_at", 
       FUNCTION(rb_gsl_vector,delete_at), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "delete_if", 
       FUNCTION(rb_gsl_vector,delete_if), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "delete", 
       FUNCTION(rb_gsl_vector,delete), 1);
  /***/
  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "filescan", FUNCTION(rb_gsl_vector,filescan), 1);

  /*****/
  rb_define_method(GSL_TYPE(cgsl_vector), "eq", FUNCTION(rb_gsl_vector,eq), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "ne", FUNCTION(rb_gsl_vector,ne), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "gt", FUNCTION(rb_gsl_vector,gt), 1);
  rb_define_alias(GSL_TYPE(cgsl_vector), ">", "gt");
  rb_define_method(GSL_TYPE(cgsl_vector), "ge", FUNCTION(rb_gsl_vector,ge), 1);
  rb_define_alias(GSL_TYPE(cgsl_vector), ">=", "ge");
  rb_define_method(GSL_TYPE(cgsl_vector), "lt", FUNCTION(rb_gsl_vector,lt), 1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "<", "lt");
  rb_define_method(GSL_TYPE(cgsl_vector), "le", FUNCTION(rb_gsl_vector,le), 1);
  rb_define_alias(GSL_TYPE(cgsl_vector), "<=", "le");

  rb_define_method(GSL_TYPE(cgsl_vector), "and", FUNCTION(rb_gsl_vector,and), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "or", FUNCTION(rb_gsl_vector,or), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "xor", FUNCTION(rb_gsl_vector,xor), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "not", FUNCTION(rb_gsl_vector,not), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "all?", FUNCTION(rb_gsl_vector,all), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "none?", FUNCTION(rb_gsl_vector,none), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "any",
       FUNCTION(rb_gsl_vector,any), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "any?",
       FUNCTION(rb_gsl_vector,any2), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "where", FUNCTION(rb_gsl_vector,where), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "where2", FUNCTION(rb_gsl_vector,where2), 0);

  rb_define_method(GSL_TYPE(cgsl_vector), "add!", FUNCTION(rb_gsl_vector,add_inplace), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "sub!", FUNCTION(rb_gsl_vector,sub_inplace), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "mul!", FUNCTION(rb_gsl_vector,mul_inplace), 1);
  rb_define_method(GSL_TYPE(cgsl_vector), "div!", FUNCTION(rb_gsl_vector,div_inplace), 1);

  rb_define_singleton_method(GSL_TYPE(cgsl_vector), "zip", FUNCTION(rb_gsl_vector,zip), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "zip", FUNCTION(rb_gsl_vector,zip), -1);
  rb_define_method(GSL_TYPE(cgsl_vector), "join", FUNCTION(rb_gsl_vector,join), -1);
#ifdef GSL_1_9_LATER
  rb_define_method(GSL_TYPE(cgsl_vector), "ispos", FUNCTION(rb_gsl_vector,ispos), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "ispos?", FUNCTION(rb_gsl_vector,ispos2), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "isneg", FUNCTION(rb_gsl_vector,isneg), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "isneg?", FUNCTION(rb_gsl_vector,isneg2), 0);
#endif

#ifdef GSL_1_10_LATER
  rb_define_method(GSL_TYPE(cgsl_vector), "isnonneg", FUNCTION(rb_gsl_vector,isnonneg), 0);
  rb_define_method(GSL_TYPE(cgsl_vector), "isnonneg?", FUNCTION(rb_gsl_vector,isnonneg2), 0);
#endif
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
