/*
  matrix_complex.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_array.h"
#include "include/rb_gsl_complex.h"

enum {
  GSL_MATRIX_COMPLEX_ADD,
  GSL_MATRIX_COMPLEX_SUB,
  GSL_MATRIX_COMPLEX_MUL,
  GSL_MATRIX_COMPLEX_DIV,
};

// From complex.c
gsl_complex rb_gsl_obj_to_gsl_complex(VALUE obj, gsl_complex *z);

// From ext/vector_source.c
void get_range_beg_en_n(VALUE range, double *beg, double *en, size_t *n, int *step);

// From ext/matrix_source.c
void parse_submatrix_args(int argc, VALUE *argv, size_t size1, size_t size2,
    size_t *i, size_t *j, size_t *n1, size_t *n2);

static VALUE rb_gsl_matrix_complex_arithmetics(int flag, VALUE obj, VALUE bb)
{
  gsl_matrix *m = NULL;
  gsl_matrix_complex *cm = NULL, *cmb = NULL, *cmnew = NULL;
  gsl_complex *c = NULL, z;
  gsl_vector *v = NULL;
  gsl_vector_complex *cv = NULL, *cvnew = NULL, *cvb = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, cm);
  switch (TYPE(bb)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    z = gsl_complex_rect(NUM2DBL(bb), 0.0);
    switch (flag) {
    case GSL_MATRIX_COMPLEX_ADD:
      cmnew = make_matrix_complex_clone(cm);
      gsl_matrix_complex_add_constant(cmnew, z);
      return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
      break;
    case GSL_MATRIX_COMPLEX_SUB:
      cmnew = make_matrix_complex_clone(cm);
      gsl_matrix_complex_add_constant(cmnew, gsl_complex_negative(z));
      return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
      break;
    case GSL_MATRIX_COMPLEX_MUL:
      cmnew = make_matrix_complex_clone(cm);
      gsl_matrix_complex_scale(cmnew, z);
      return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
      break;
    case GSL_MATRIX_COMPLEX_DIV:
      cmnew = make_matrix_complex_clone(cm);
      gsl_matrix_complex_scale(cmnew, gsl_complex_inverse(z));
      return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
      break;
    default:
      rb_raise(rb_eRuntimeError, "operation not defined");
      break;
    }
    break;
  default:
    if (rb_obj_is_kind_of(bb, cgsl_matrix_complex)) {
      Data_Get_Struct(bb, gsl_matrix_complex, cmb);
      switch (flag) {
      case GSL_MATRIX_COMPLEX_ADD:
  cmnew = make_matrix_complex_clone(cm);
  gsl_matrix_complex_add(cmnew, cmb);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
  break;
      case GSL_MATRIX_COMPLEX_SUB:
  cmnew = make_matrix_complex_clone(cm);
  gsl_matrix_complex_sub(cmnew,cmb);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
  break;
      case GSL_MATRIX_COMPLEX_MUL:
  cmnew = make_matrix_complex_clone(cm);
  gsl_matrix_complex_mul_elements(cmnew, cmb);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
  break;
      case GSL_MATRIX_COMPLEX_DIV:
  cmnew = make_matrix_complex_clone(cm);
  gsl_matrix_complex_div_elements(cmnew, cmb);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
  break;
      default:
  rb_raise(rb_eRuntimeError, "operation not defined");
  break;
      }
    } else if (rb_obj_is_kind_of(bb, cgsl_matrix)) {
      Data_Get_Struct(bb, gsl_matrix, m);
      cmb = matrix_to_complex(m);
      cmnew = make_matrix_complex_clone(cm);
      switch (flag) {
      case GSL_MATRIX_COMPLEX_ADD:
  gsl_matrix_complex_add(cmnew, cmb);
  break;
      case GSL_MATRIX_COMPLEX_SUB:
  gsl_matrix_complex_sub(cmnew,cmb);
  break;
      case GSL_MATRIX_COMPLEX_MUL:
  gsl_matrix_complex_mul_elements(cmnew, cmb);
  break;
      case GSL_MATRIX_COMPLEX_DIV:
  gsl_matrix_complex_div_elements(cmnew, cmb);
  break;
      default:
  rb_raise(rb_eRuntimeError, "operation not defined");
  break;
      }
      return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
    } else if (rb_obj_is_kind_of(bb, cgsl_complex)) {
      Data_Get_Struct(bb, gsl_complex, c);
      switch (flag) {
      case GSL_MATRIX_COMPLEX_ADD:
  cmnew = make_matrix_complex_clone(cm);
  gsl_matrix_complex_add_constant(cmnew, *c);
  break;
      case GSL_MATRIX_COMPLEX_SUB:
  cmnew = make_matrix_complex_clone(cm);
  gsl_matrix_complex_add_constant(cmnew, gsl_complex_negative(*c));
  break;
      case GSL_MATRIX_COMPLEX_MUL:
  cmnew = make_matrix_complex_clone(cm);
  gsl_matrix_complex_scale(cmnew, *c);
  break;
      case GSL_MATRIX_COMPLEX_DIV:
  cmnew = make_matrix_complex_clone(cm);
  gsl_matrix_complex_scale(cmnew, gsl_complex_inverse(*c));
  break;
      default:
  rb_raise(rb_eRuntimeError, "operation not defined");
  break;
      }
      return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
    } else if (rb_obj_is_kind_of(bb, cgsl_vector)) {
      Data_Get_Struct(bb, gsl_vector, v);
      switch (flag) {
      case GSL_MATRIX_COMPLEX_MUL:
  cvb = vector_to_complex(v);
  cvnew = gsl_vector_complex_alloc(v->size);
  if (cvnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
  gsl_matrix_complex_mul_vector(cvnew, cm, cvb);
  gsl_vector_complex_free(cvb);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cvnew);
  break;
      default:
  rb_raise(rb_eRuntimeError, 
     "operation is not defined %s and Matrix_Complex",
     rb_class2name(CLASS_OF(bb)));
  break;
      }
    } else if (rb_obj_is_kind_of(bb, cgsl_vector_complex)) {
      if (!VECTOR_COMPLEX_COL_P(bb))
  rb_raise(rb_eTypeError, 
     "Operation is not defined with %s (Vector::Complex::Col expected)", 
     rb_class2name(CLASS_OF(bb)));
      Data_Get_Struct(bb, gsl_vector_complex, cv);
      switch (flag) {
      case GSL_MATRIX_COMPLEX_MUL:
  cvnew = gsl_vector_complex_alloc(cv->size);
  if (cvnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
  gsl_matrix_complex_mul_vector(cvnew, cm, cv);
  return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, cvnew);;
      default:
  rb_raise(rb_eRuntimeError, 
     "operation is not defined %s and Matrix_Complex",
     rb_class2name(CLASS_OF(bb)));
  break;
      }
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s", rb_class2name(CLASS_OF(bb)));
    }
  }
  /* never reach here */
  return Qnil;
}

static VALUE rb_gsl_matrix_complex_new(VALUE klass, VALUE s1, VALUE s2)
{
  gsl_matrix_complex *m = NULL;
  CHECK_FIXNUM(s1); CHECK_FIXNUM(s2);
  m = gsl_matrix_complex_calloc(FIX2INT(s1),  FIX2INT(s2));
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
  return Data_Wrap_Struct(klass, 0, gsl_matrix_complex_free, m);
}

static VALUE rb_gsl_matrix_complex_eye(int argc, VALUE *argv, VALUE klass)
{
  size_t n, i;
  gsl_matrix_complex *m = NULL;
  gsl_complex z, *pz = &z;
  switch (argc) {
  case 1:
    CHECK_FIXNUM(argv[0]);
    n = FIX2INT(argv[0]);
    z = gsl_complex_rect(1.0, 0.0);
    break;
  case 2:
    CHECK_FIXNUM(argv[0]);
    n = FIX2INT(argv[0]);
    switch (TYPE(argv[1])) {
    case T_FIXNUM:
    case T_BIGNUM:
    case T_FLOAT:
      z = gsl_complex_rect(NUM2DBL(argv[1]), 0.0);
      break;
    case T_ARRAY:
      //      if (RARRAY(argv[1])->len < 2) rb_raise(rb_eArgError, "wrong argument");
      if (RARRAY_LEN(argv[1]) < 2) rb_raise(rb_eArgError, "wrong argument");
      z = gsl_complex_rect(NUM2DBL(rb_ary_entry(argv[1], 0)),
         NUM2DBL(rb_ary_entry(argv[1], 1)));
      break;
    default:
      if (rb_obj_is_kind_of(argv[1], cgsl_complex)) {
  Data_Get_Struct(argv[1], gsl_complex, pz);
        z = *pz;
      } else {
  rb_raise(rb_eTypeError, 
     "wrong argument type %s", rb_class2name(CLASS_OF(argv[1])));
      }
      break;
    }
    break;
  case 3:
    CHECK_FIXNUM(argv[0]);
    n = FIX2INT(argv[0]);
    Need_Float(argv[1]); Need_Float(argv[2]);
    z = gsl_complex_rect(NUM2DBL(argv[1]), NUM2DBL(argv[2]));
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1, 2, or 3)", argc);
    break;
  }
  m = gsl_matrix_complex_calloc(n, n);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
  for (i = 0; i < n; i++) gsl_matrix_complex_set(m, i, i, z);
  return  Data_Wrap_Struct(klass, 0, gsl_matrix_complex_free, m);
}

static VALUE rb_gsl_matrix_complex_identity(VALUE klass, VALUE nn)
{
  size_t n, i;
  gsl_matrix_complex *m = NULL;
  gsl_complex z;
  CHECK_FIXNUM(nn);
  n = FIX2INT(nn);
  m = gsl_matrix_complex_calloc(n, n);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_calloc failed");
  z = gsl_complex_rect(1.0, 0.0);
  for (i = 0; i < n; i++) gsl_matrix_complex_set(m, i, i, z);
  return  Data_Wrap_Struct(klass, 0, gsl_matrix_complex_free, m);
}

static VALUE rb_gsl_matrix_complex_set_zero(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  gsl_matrix_complex_set_zero(m);
  return obj;
}

static VALUE rb_gsl_matrix_complex_set_identity(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  gsl_matrix_complex_set_identity(m);
  return obj;
}

static VALUE rb_gsl_matrix_complex_set_all(VALUE obj, VALUE s)
{
  gsl_matrix_complex *m = NULL;
  gsl_complex c, *z = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  switch (TYPE(s)) {
  case T_FLOAT:
  case T_BIGNUM:
  case T_FIXNUM:
    c.dat[0] = NUM2DBL(s);
    c.dat[1] = 0.0;
    gsl_matrix_complex_set_all(m, c);
    break;
  case T_ARRAY:
    c.dat[0] = NUM2DBL(rb_ary_entry(s, 0));
    c.dat[1] = NUM2DBL(rb_ary_entry(s, 1));
    gsl_matrix_complex_set_all(m, c);
    break;
  default:
    if (rb_obj_is_kind_of(s, cgsl_complex)) {
      Data_Get_Struct(s, gsl_complex, z);
      gsl_matrix_complex_set_all(m, *z);
    } else {
      rb_raise(rb_eTypeError, 
         "wrong argument type %s", rb_class2name(CLASS_OF(s)));
    }
    break;
  }
  return obj;
}

void rb_gsl_vector_complex_set_subvector(int argc, VALUE *argv, gsl_vector_complex *v, VALUE other);
static VALUE rb_gsl_matrix_complex_set(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL, *mother;
  gsl_matrix_complex_view mv;
  gsl_vector_complex_view vv;
  gsl_complex tmp;
  VALUE other, row, row_set_argv[2];
  int ii, ij, step;
  size_t i, j, k, n1, n2, nother;
  double beg, end;

  if(argc < 1 || argc > 5) {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1-5)", argc);
  }

  Data_Get_Struct(obj, gsl_matrix_complex, m);
  other = argv[argc-1];

  if(argc == 1) {
    // m[] = x
    gsl_matrix_complex_set_all(m, rb_gsl_obj_to_gsl_complex(other, NULL));
  } else if(argc == 3 && TYPE(argv[0]) == T_FIXNUM && TYPE(argv[1]) == T_FIXNUM) {
    // m[i,j] = x
    ii = FIX2INT(argv[0]);
    ij = FIX2INT(argv[1]);
    if(ii < 0) ii += m->size1;
    if(ij < 0) ij += m->size2;
    gsl_matrix_complex_set(m, (size_t)ii, (size_t)ij, rb_gsl_obj_to_gsl_complex(other, NULL));
  } else if(TYPE(argv[0]) == T_ARRAY) {
    // m.set(row0...)
    row_set_argv[0] = INT2FIX(0);
    row_set_argv[1] = INT2FIX(m->size2);

    for(k = 0; (int) k < argc && k < m->size1; k++) {
      vv = gsl_matrix_complex_row(m, k);
      rb_gsl_vector_complex_set_subvector(2, row_set_argv, &vv.vector, argv[k]);
    }
  } else {
    // x -> assignment to m.submatrix(i...)
    parse_submatrix_args(argc-1, argv, m->size1, m->size2, &i, &j, &n1, &n2);
    if(n1 == 0) n1 = 1;
    if(n2 == 0) n2 = 1;
    mv = gsl_matrix_complex_submatrix(m, i, j, n1, n2);
    if(rb_obj_is_kind_of(other, cgsl_matrix_complex)) {
      // m[...] = m_other
      Data_Get_Struct(other, gsl_matrix_complex, mother);
      if(n1 * n2 != mother->size1 * mother->size2) {
        rb_raise(rb_eRangeError, "sizes do not match (%d x %d != %d x %d)",
     (int) n1, (int) n2, (int) mother->size1, (int) mother->size2);
      }
      // TODO Change to gsl_matrix_memmove if/when GSL has such a function
      // because gsl_matrix_memcpy does not handle overlapping regions (e.g.
      // Views) well.
      gsl_matrix_complex_memcpy(&mv.matrix, mother);
    } else if(rb_obj_is_kind_of(other, rb_cArray)) {
      row_set_argv[0] = INT2FIX(0);
      row_set_argv[1] = INT2FIX(n2);

      if(n1 == 1) {
        // m[...] = [col0, ...] # single row
        vv = gsl_matrix_complex_row(&mv.matrix, 0);
        rb_gsl_vector_complex_set_subvector(2, row_set_argv, &vv.vector, other);
      } else {
        // m[...] = [[row0], [row1], ...] # multiple rows
        if((int) n1 != RARRAY_LEN(other)) {
          rb_raise(rb_eRangeError, "row counts do not match (%d != %d)",
       (int) n1, (int) RARRAY_LEN(other));
        }
        for(k = 0; k < n1; k++) {
          vv = gsl_matrix_complex_row(&mv.matrix, k);
          row = rb_ary_entry(other, k);
          rb_gsl_vector_complex_set_subvector(2, row_set_argv, &vv.vector, row);
        }
      }
    } else if(rb_obj_is_kind_of(other, rb_cRange)) {
      // m[...] = beg..end
      get_range_beg_en_n(other, &beg, &end, &nother, &step);
      if(n1 * n2 != nother) {
        rb_raise(rb_eRangeError, "sizes do not match (%d x %d != %d)", (int) n1, (int) n2, (int) nother);
      }
      tmp = gsl_complex_rect(beg, 0.0);
      for(k = 0; k < nother; k++) {
        gsl_matrix_complex_set(&mv.matrix, k / n2, k % n2, tmp);
        GSL_SET_REAL(&tmp, GSL_REAL(tmp) + step);
      }
    } else {
      // m[...] = x
      gsl_matrix_complex_set_all(&mv.matrix, rb_gsl_obj_to_gsl_complex(other, NULL));
    }
  }

  return obj;
}

static VALUE rb_gsl_matrix_complex_set_row(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  size_t i, k;
  gsl_complex z, *pz = &z;
  if (argc < 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 2)",
       argc);
  CHECK_FIXNUM(argv[0]);
  Data_Get_Struct(obj, gsl_matrix_complex, A);
  i = FIX2INT(argv[0]);
  for (k = 1; (int) k < argc; k++) {
    if (k-1 >= A->size1) break;
    switch (TYPE(argv[k])) {
    case T_ARRAY:
      z = ary2complex(argv[k]);
      break;
    default:
      CHECK_COMPLEX(argv[k]);
      Data_Get_Struct(argv[k], gsl_complex, pz);
      z = *pz;
      break;
    }
    gsl_matrix_complex_set(A, i, k-1, z);
  }
  return obj;
}

static VALUE rb_gsl_matrix_complex_set_col(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  int j, k;
  gsl_complex z, *pz = &z;
  if (argc < 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 2)",
       argc);
  CHECK_FIXNUM(argv[0]);
  Data_Get_Struct(obj, gsl_matrix_complex, A);
  j = FIX2INT(argv[0]);
  for (k = 1; k < argc; k++) {
    if (k-1 >= (int) A->size2) break;
    switch (TYPE(argv[k])) {
    case T_ARRAY:
      z = ary2complex(argv[k]);
      break;
    default:
      CHECK_COMPLEX(argv[k]);
      Data_Get_Struct(argv[k], gsl_complex, pz);
      z = *pz;
      break;
    }
    gsl_matrix_complex_set(A, k-1, j, z);
  }
  return obj;
}

static VALUE rb_gsl_matrix_complex_submatrix(int argc, VALUE *argv, VALUE obj);
static VALUE rb_gsl_matrix_complex_get(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  gsl_complex *c = NULL;
  VALUE retval;
  int ii, ij;

  if(argc == 2 && TYPE(argv[0]) == T_FIXNUM && TYPE(argv[1]) == T_FIXNUM) {
    // m[i,j]
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    ii = FIX2INT(argv[0]);
    ij = FIX2INT(argv[1]);
    if(ii < 0) ii += m->size1;
    if(ij < 0) ij += m->size2;
    c = ALLOC(gsl_complex);
    *c = gsl_matrix_complex_get(m, (size_t)ii, (size_t)ij); 
    retval = Data_Wrap_Struct(cgsl_complex, 0, free, c);
  } else if(argc == 1 && TYPE(argv[0]) == T_FIXNUM) {
    // m[i]
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    ii = FIX2INT(argv[0]);
    if(ii < 0) ii += m->size1 * m->size2;
    c = ALLOC(gsl_complex);
    *c = gsl_matrix_complex_get(m, (size_t)(ii / m->size2), (size_t)(ii % m->size2));
    retval = Data_Wrap_Struct(cgsl_complex, 0, free, c);
  } else if(argc == 1 && TYPE(argv[0]) == T_ARRAY) {
    // m[[i,j]], to have parity with Real Matrix types
    if(RARRAY_LEN(argv[0]) == 2) {
      Data_Get_Struct(obj, gsl_matrix_complex, m);
      ii = FIX2INT(RARRAY_PTR(argv[0])[0]);
      ij = FIX2INT(RARRAY_PTR(argv[0])[1]);
      if(ii < 0) ii += m->size1;
      if(ij < 0) ij += m->size2;
      c = ALLOC(gsl_complex);
      *c = gsl_matrix_complex_get(m, (size_t)ii, (size_t)ij); 
      retval = Data_Wrap_Struct(cgsl_complex, 0, free, c);
    } else {
      rb_raise(rb_eArgError, "Array index must have length 2, not %d", (int) RARRAY_LEN(argv[0]));
    }
  } else {
    retval = rb_gsl_matrix_complex_submatrix(argc, argv, obj);
  }

  return retval;
}

static void rb_gsl_matrix_complex_collect_native(gsl_matrix_complex *src, gsl_matrix_complex *dst)
{
  VALUE vz;
  gsl_complex * zp;
  size_t i, j;
  for (i = 0; i < src->size1; i++) {
    for (j = 0; j < src->size2; j++) {
      vz = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, zp);
      *zp = gsl_matrix_complex_get(src, i, j);
      vz = rb_yield(vz);
      CHECK_COMPLEX(vz);
      Data_Get_Struct(vz, gsl_complex, zp);
      gsl_matrix_complex_set(dst, i, j, *zp);
    }
  }
}

static VALUE rb_gsl_matrix_complex_collect(VALUE obj)
{
  gsl_matrix_complex *m = NULL, *mnew;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
  rb_gsl_matrix_complex_collect_native(m, mnew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mnew);
}

static VALUE rb_gsl_matrix_complex_collect_bang(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  rb_gsl_matrix_complex_collect_native(m, m);
  return obj;
}

static VALUE rb_gsl_matrix_complex_to_a(VALUE obj)
{
  gsl_matrix_complex *m;
  gsl_complex *c;
  VALUE ma, ra;
  size_t i, j;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  ma = rb_ary_new2(m->size1);
  for(i=0; i < m->size1; i++) {
    ra = rb_ary_new2(m->size2);
    rb_ary_store(ma, i, ra);
    for(j=0; j < m->size2; j++) {
      c = ALLOC(gsl_complex);
      *c = gsl_matrix_complex_get(m, i, j); 
      rb_ary_store(ra, j, Data_Wrap_Struct(cgsl_complex, 0, free, c));
    }
  }
  return ma;
}

static VALUE rb_gsl_matrix_complex_ptr(VALUE obj, VALUE i, VALUE j)
{
  gsl_matrix_complex *m = NULL;
  gsl_complex *c = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  c = gsl_matrix_complex_ptr(m, FIX2INT(i), FIX2INT(j));
  return Data_Wrap_Struct(cgsl_complex, 0, NULL, c);
}

static VALUE rb_gsl_matrix_complex_printf(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *h = NULL;
  int status;
  Data_Get_Struct(obj, gsl_matrix_complex, h);
  if (argc == 1) {
    Check_Type(argv[0], T_STRING);
    status = gsl_matrix_complex_fprintf(stdout, h, STR2CSTR(argv[0]));
  } else {
    status = gsl_matrix_complex_fprintf(stdout, h, "%g");
  }
  return INT2FIX(status);
}

static VALUE rb_gsl_matrix_complex_print(VALUE obj)
{
  gsl_matrix_complex *h = NULL;
  gsl_complex *z = NULL;
  size_t i, j;
  Data_Get_Struct(obj, gsl_matrix_complex, h);
  printf("[ ");
  for (i = 0; i < h->size1; i++) {
    if (i != 0) printf("  ");
    for (j = 0; j < h->size2; j++) {
      z = gsl_matrix_complex_ptr(h, i, j);
      printf("[%4.3e %4.3e] ", GSL_REAL(*z), GSL_IMAG(*z));
    }
    if (i != h->size1 -1) printf("\n");
    else printf("]\n");
  }
  return obj;
}

static VALUE rb_gsl_matrix_complex_to_s(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  char buf[64];
  size_t i, j;
  VALUE str;
  gsl_complex z;
  int max_rows = 4;
  int max_cols = 4;

  switch(argc){
    case 2: max_cols = NUM2INT(argv[1]);
    case 1: max_rows = NUM2INT(argv[0]);
    case 0: break;
    default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0, 1, or 2)", argc);
  }

  Data_Get_Struct(obj, gsl_matrix_complex, m);
  if (m->size1 == 0 && m->size2 == 0) return rb_str_new2("[ ]");
  str = rb_str_new2("[ ");
  for (i = 0; i < m->size1; i++) {
    if (i != 0) {
      rb_str_cat(str, "\n  ", 3);
    }
    for (j = 0; j < m->size2; j++) {
      z = gsl_matrix_complex_get(m, i, j);
      sprintf(buf,
          "%s[ %4.3e %4.3e ]", (j==0) ? "" : " ", GSL_REAL(z), GSL_IMAG(z));
      rb_str_cat(str, buf, strlen(buf));
      // if too many cols
      if ((int) j >= max_cols-1 && j != m->size2-1) {
        rb_str_cat(str, " ...", 4);
        break;
      }
    }
    // if too many rows
    if ((int) i >= max_rows-1 && i != m->size1-1) {
      rb_str_cat(str, "\n  ...", 6);
      break;
    }
  }
  rb_str_cat(str, " ]", 2);
  return str;
}

static VALUE rb_gsl_matrix_complex_inspect(int argc, VALUE *argv, VALUE obj)
{
  VALUE str;
  char buf[128];
  gsl_matrix_complex *m;

  Data_Get_Struct(obj, gsl_matrix_complex, m);
  sprintf(buf, "#<%s[%lu,%lu]:%#lx>\n", rb_class2name(CLASS_OF(obj)), m->size1, m->size2, NUM2ULONG(rb_obj_id(obj)));
  str = rb_str_new2(buf);
  return rb_str_concat(str, rb_gsl_matrix_complex_to_s(argc, argv, obj));
}

static VALUE rb_gsl_matrix_complex_fwrite(VALUE obj, VALUE io)
{
  gsl_matrix_complex *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_matrix_complex, h);
  fp = rb_gsl_open_writefile(io, &flag);
  status = gsl_matrix_complex_fwrite(fp, h);
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE rb_gsl_matrix_complex_fread(VALUE obj, VALUE io)
{
  gsl_matrix_complex *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_matrix_complex, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = gsl_matrix_complex_fread(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_matrix_complex_fprintf(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  if (argc != 1 && argc != 2) {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
  Data_Get_Struct(obj, gsl_matrix_complex, h);
  fp = rb_gsl_open_writefile(argv[0], &flag);
  if (argc == 2) {
    Check_Type(argv[1], T_STRING);
    status = gsl_matrix_complex_fprintf(fp, h, STR2CSTR(argv[1]));
  } else {
    status = gsl_matrix_complex_fprintf(fp, h, "%g");
  }
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE rb_gsl_matrix_complex_fscanf(VALUE obj, VALUE io)
{
  gsl_matrix_complex *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_matrix_complex, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = gsl_matrix_complex_fscanf(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

gsl_matrix_complex_view* gsl_matrix_complex_view_alloc()
{
  gsl_matrix_complex_view *vv = NULL;
  vv = ALLOC(gsl_matrix_complex_view);
  if (vv == NULL) rb_raise(rb_eRuntimeError, "malloc failed");
  return vv;
}

void gsl_matrix_complex_view_free(gsl_matrix_view * vv)
{
  free((gsl_matrix_complex_view *) vv);
}

/* singleton */
static VALUE rb_gsl_matrix_complex_memcpy(VALUE obj, VALUE dst, VALUE src)
{
  gsl_matrix_complex *m, *dest;
  CHECK_MATRIX_COMPLEX(dst);
  CHECK_MATRIX_COMPLEX(src);
  Data_Get_Struct(dst, gsl_matrix_complex, dest);
  Data_Get_Struct(src, gsl_matrix_complex, m);
  gsl_matrix_complex_memcpy(dest, m);
  return dst;
}

static VALUE rb_gsl_matrix_complex_clone(VALUE obj)
{
  gsl_matrix_complex *m, *mnew = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
  if (mnew == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
  gsl_matrix_complex_memcpy(mnew, m);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mnew);
}

static VALUE rb_gsl_matrix_complex_swap_rows(VALUE obj, VALUE i, VALUE j)
{
  gsl_matrix_complex *m = NULL;
  CHECK_FIXNUM(i); CHECK_FIXNUM(j);
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  gsl_matrix_complex_swap_rows(m, FIX2INT(i), FIX2INT(j));
  return obj;
}

static VALUE rb_gsl_matrix_complex_swap_columns(VALUE obj, VALUE i, VALUE j)
{
  gsl_matrix_complex *m = NULL;
  CHECK_FIXNUM(i); CHECK_FIXNUM(j);
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  gsl_matrix_complex_swap_columns(m, FIX2INT(i), FIX2INT(j));
  return obj;
}

static VALUE rb_gsl_matrix_complex_swap_rowcol(VALUE obj, VALUE i, VALUE j)
{
  gsl_matrix_complex *m = NULL;
  CHECK_FIXNUM(i); CHECK_FIXNUM(j);
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  gsl_matrix_complex_swap_rowcol(m, FIX2INT(i), FIX2INT(j));
  return obj;
}

static VALUE rb_gsl_matrix_complex_transpose(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  gsl_matrix_complex_transpose(m);
  return obj;
}

static VALUE rb_gsl_matrix_complex_isnull(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  if (gsl_matrix_complex_isnull(m)) return Qtrue;
  else return Qfalse;
}

static VALUE rb_gsl_matrix_complex_add(VALUE obj, VALUE mmb)
{
  return rb_gsl_matrix_complex_arithmetics(GSL_MATRIX_COMPLEX_ADD, obj, mmb);
}

static VALUE rb_gsl_matrix_complex_sub(VALUE obj, VALUE mmb)
{
  return rb_gsl_matrix_complex_arithmetics(GSL_MATRIX_COMPLEX_SUB, obj, mmb);
}

static VALUE rb_gsl_matrix_complex_mul_elements(VALUE obj, VALUE mmb)
{
  return rb_gsl_matrix_complex_arithmetics(GSL_MATRIX_COMPLEX_MUL, obj, mmb);
}

static VALUE rb_gsl_matrix_complex_div_elements(VALUE obj, VALUE mmb)
{
  return rb_gsl_matrix_complex_arithmetics(GSL_MATRIX_COMPLEX_DIV, obj, mmb);
}

static VALUE rb_gsl_matrix_complex_scale(VALUE obj, VALUE s)
{
  return rb_gsl_matrix_complex_arithmetics(GSL_MATRIX_COMPLEX_MUL, obj, s);
}

static VALUE rb_gsl_matrix_complex_scale_bang(VALUE obj, VALUE s)
{
  gsl_matrix_complex *m;
  gsl_complex c, *z = &c;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  switch (TYPE(s)) {
  case T_FIXNUM:
  case T_FLOAT:
    GSL_SET_REAL(z, NUM2DBL(s));
    GSL_SET_IMAG(z, 0.0);
    break;
  default:
    CHECK_COMPLEX(s);
    Data_Get_Struct(s, gsl_complex, z);
    break;
  }
  gsl_matrix_complex_scale(m, *z);
  return obj;
}

static VALUE rb_gsl_matrix_complex_mul(VALUE obj, VALUE mb)
{
  gsl_matrix_complex *cm = NULL, *cmb = NULL, *cmnew = NULL;
  gsl_matrix *m = NULL;
  gsl_vector *v = NULL;
  gsl_vector_complex *vc = NULL, *vcnew = NULL;
  gsl_complex a, b;
  int flag = 0;
  if (COMPLEX_P(mb) || TYPE(mb) == T_FIXNUM || TYPE(mb) == T_FLOAT)
    return rb_gsl_matrix_complex_mul_elements(obj, mb);
  Data_Get_Struct(obj, gsl_matrix_complex, cm);
  if (VECTOR_P(mb)) {
    Data_Get_Struct(mb, gsl_vector, v);
    vc = vector_to_complex(v);
    vcnew = gsl_vector_complex_calloc(vc->size);
    a.dat[0] = 1.0; a.dat[1] = 0.0;
    b.dat[0] = 0.0; b.dat[1] = 0.0;
    gsl_blas_zgemv(CblasNoTrans, a, cm, vc, b, vcnew);
    gsl_vector_complex_free(vc);
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vcnew);
  }
  if (VECTOR_COMPLEX_P(mb)) {
    Data_Get_Struct(mb, gsl_vector_complex, vc);
    vcnew = gsl_vector_complex_calloc(vc->size);
    a.dat[0] = 1.0; a.dat[1] = 0.0;
    b.dat[0] = 0.0; b.dat[1] = 0.0;
    gsl_blas_zgemv(CblasNoTrans, a, cm, vc, b, vcnew);
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vcnew);
  }
  if (MATRIX_P(mb)) {
    Data_Get_Struct(mb, gsl_matrix, m);
    cmb = matrix_to_complex(m);
    flag = 1;
  } else {
    CHECK_MATRIX_COMPLEX(mb);
    Data_Get_Struct(mb, gsl_matrix_complex, cmb);
  }
  cmnew = gsl_matrix_complex_alloc(cm->size1, cm->size2);
  if (cmnew == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
  gsl_matrix_complex_mul(cmnew, cm, cmb);
  if (flag == 1) gsl_matrix_complex_free(cmb);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
}

static VALUE rb_gsl_matrix_complex_mul2(VALUE obj, VALUE mb)
{
  gsl_matrix_complex *cm = NULL, *cmb = NULL, *cmnew = NULL;
  gsl_matrix *m = NULL;
  int flag = 0;
  Data_Get_Struct(obj, gsl_matrix_complex, cm);
  if (MATRIX_P(mb)) {
    Data_Get_Struct(mb, gsl_matrix, m);
    cmb = matrix_to_complex(m);
    flag = 1;
  } else {
    CHECK_MATRIX_COMPLEX(mb);
    Data_Get_Struct(mb, gsl_matrix_complex, cmb);
  }
  cmnew = gsl_matrix_complex_alloc(cm->size1, cm->size2);
  if (cmnew == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
  gsl_matrix_complex_mul(cmnew, cm, cmb);
  gsl_matrix_complex_memcpy(cm, cmnew);
  if (flag == 1) gsl_matrix_complex_free(cmb);
  return obj;
}

static VALUE rb_gsl_matrix_complex_add_diagonal(VALUE obj, VALUE s)
{
  gsl_matrix_complex *m = NULL;
  gsl_complex c, *z = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  switch (TYPE(s)) {
  case T_FLOAT:
  case T_BIGNUM:
  case T_FIXNUM:
    c.dat[0] = NUM2DBL(s);
    c.dat[1] = 0.0;
    gsl_matrix_complex_add_diagonal(m, c);
    break;
  case T_ARRAY:
    c.dat[0] = NUM2DBL(rb_ary_entry(s, 0));
    c.dat[1] = NUM2DBL(rb_ary_entry(s, 1));
    gsl_matrix_complex_add_diagonal(m, c);
    break;
  default:
    if (rb_obj_is_kind_of(s, cgsl_complex)) {
      Data_Get_Struct(s, gsl_complex, z);
      gsl_matrix_complex_add_diagonal(m, *z);
    } else {
      rb_raise(rb_eTypeError, 
         "wrong argument type %s", rb_class2name(CLASS_OF(s)));
    }
    break;
  }
  return obj;
}

static VALUE rb_gsl_matrix_complex_submatrix(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  gsl_matrix_complex_view *mv = NULL;
  gsl_vector_complex_view *vv = NULL;
  size_t i, j, n1, n2;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  parse_submatrix_args(argc, argv, m->size1, m->size2, &i, &j, &n1, &n2);
  if(n1 == 0) {
    vv = ALLOC(gsl_vector_complex_view);
    *vv = gsl_matrix_complex_subrow(m, i, j, n2);
    return Data_Wrap_Struct(cgsl_vector_complex_view, 0, free, vv);
  }
  else if(n2 == 0) {
    vv = ALLOC(gsl_vector_complex_view);
    *vv = gsl_matrix_complex_subcolumn(m, j, i, n1);
    return Data_Wrap_Struct(cgsl_vector_complex_col_view, 0, free, vv);
  } else {
    mv = ALLOC(gsl_matrix_complex_view);
    *mv = gsl_matrix_complex_submatrix(m, i, j, n1, n2);
    return Data_Wrap_Struct(cgsl_matrix_complex_view, 0, free, mv);
  }
}

static VALUE rb_gsl_matrix_complex_row(VALUE obj, VALUE i)
{
  gsl_matrix_complex *m = NULL;
  gsl_vector_complex_view *vv = NULL;
  CHECK_FIXNUM(i);
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  vv = gsl_vector_complex_view_alloc();
  *vv = gsl_matrix_complex_row(m, FIX2INT(i));
  return Data_Wrap_Struct(cgsl_vector_complex_view, 0, gsl_vector_complex_view_free, vv);
}

static VALUE rb_gsl_matrix_complex_column(VALUE obj, VALUE i)
{
  gsl_matrix_complex *m = NULL;
  gsl_vector_complex_view *vv = NULL;
  CHECK_FIXNUM(i);
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  vv = gsl_vector_complex_view_alloc();
  *vv = gsl_matrix_complex_column(m, FIX2INT(i));
  return Data_Wrap_Struct(cgsl_vector_complex_col_view, 0, gsl_vector_complex_view_free, vv);
}

static VALUE rb_gsl_matrix_complex_diagonal(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  gsl_vector_complex_view *vv = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  vv = gsl_vector_complex_view_alloc();
  *vv = gsl_matrix_complex_diagonal(m);
  return Data_Wrap_Struct(cgsl_vector_complex_view, 0, gsl_vector_complex_view_free, vv);
}

static VALUE rb_gsl_matrix_complex_set_diagonal(VALUE obj, VALUE diag)
{
  gsl_matrix_complex *m = NULL;
  gsl_vector_complex *v;
  size_t i;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  if (VECTOR_COMPLEX_P(diag)) {
    Data_Get_Struct(diag, gsl_vector_complex, v);
    for (i = 0; i < m->size1; i++) gsl_matrix_complex_set(m, i, i, gsl_vector_complex_get(v, i));
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Vector_Complex or Array expected)",
       rb_class2name(CLASS_OF(diag)));
  }
  return obj;
}

static VALUE rb_gsl_matrix_complex_subdiagonal(VALUE obj, VALUE i)
{
  gsl_matrix_complex *m = NULL;
  gsl_vector_complex_view *vv = NULL;
  CHECK_FIXNUM(i);
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  vv = gsl_vector_complex_view_alloc();
  *vv = gsl_matrix_complex_subdiagonal(m, FIX2INT(i));
  return Data_Wrap_Struct(cgsl_vector_complex_view, 0, gsl_vector_complex_view_free, vv);
}

static VALUE rb_gsl_matrix_complex_superdiagonal(VALUE obj, VALUE i)
{
  gsl_matrix_complex *m = NULL;
  gsl_vector_complex_view *vv = NULL;
  CHECK_FIXNUM(i);
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  vv = gsl_vector_complex_view_alloc();
  *vv = gsl_matrix_complex_superdiagonal(m, FIX2INT(i));
  return Data_Wrap_Struct(cgsl_vector_complex_view, 0, gsl_vector_complex_view_free, vv);
}

static VALUE rb_gsl_matrix_complex_coerce(VALUE obj, VALUE other)
{
  gsl_matrix_complex *cm = NULL, *cmnew = NULL;
  gsl_matrix *m = NULL;
  gsl_complex z;
  VALUE vcm;
  Data_Get_Struct(obj, gsl_matrix_complex, cm);
  switch (TYPE(other)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    cmnew = gsl_matrix_complex_alloc(cm->size1, cm->size2);
    if (cmnew == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
    GSL_SET_REAL(&z, NUM2DBL(other));
    GSL_SET_IMAG(&z, 0.0);
    gsl_matrix_complex_set_all(cmnew, z);
    vcm = Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
    return rb_ary_new3(2, vcm, obj);
    break;
  default:
    if (rb_obj_is_kind_of(other, cgsl_matrix)) {
      Data_Get_Struct(other, gsl_matrix, m);
      cmnew = matrix_to_complex(m);
      vcm = Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
      return rb_ary_new3(2, vcm, obj);
    } else {
      rb_raise(rb_eTypeError, "cannot coerce %s to GSL::Matrix::Complex",
         rb_class2name(CLASS_OF(other)));
    }
    break;
  }
}

static VALUE rb_gsl_matrix_complex_real(VALUE obj)
{
  gsl_matrix_complex *cm = NULL;
  gsl_matrix *m = NULL;
  gsl_complex z;
  size_t i, j;
  Data_Get_Struct(obj, gsl_matrix_complex, cm);
  m = gsl_matrix_alloc(cm->size1, cm->size2);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
  for (i = 0; i < cm->size1; i++) {
    for (j = 0; j < cm->size2; j++) {
      z = gsl_matrix_complex_get(cm, i, j);
      gsl_matrix_set(m, i, j, GSL_REAL(z));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, m);
}

static VALUE rb_gsl_matrix_complex_imag(VALUE obj)
{
  gsl_matrix_complex *cm = NULL;
  gsl_matrix *m = NULL;
  gsl_complex z;
  size_t i, j;
  Data_Get_Struct(obj, gsl_matrix_complex, cm);
  m = gsl_matrix_alloc(cm->size1, cm->size2);
  if (m == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
  for (i = 0; i < cm->size1; i++) {
    for (j = 0; j < cm->size2; j++) {
      z = gsl_matrix_complex_get(cm, i, j);
      gsl_matrix_set(m, i, j, GSL_IMAG(z));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, m);
}

static void gsl_matrix_complex_conjugate(gsl_matrix_complex *cm)
{
  gsl_complex z;
  size_t i, j;
  for (i = 0; i < cm->size1; i++) {
    for (j = 0; j < cm->size2; j++) {
      z = gsl_matrix_complex_get(cm, i, j);
      gsl_matrix_complex_set(cm, i, j, gsl_complex_conjugate(z));
    }
  }
}

static void gsl_matrix_complex_conjugate2(gsl_matrix_complex *cmnew, gsl_matrix_complex *cm)
{
  gsl_complex z;
  size_t i, j;
  for (i = 0; i < cm->size1; i++) {
    for (j = 0; j < cm->size2; j++) {
      z = gsl_matrix_complex_get(cm, i, j);
      gsl_matrix_complex_set(cmnew, i, j, gsl_complex_conjugate(z));
    }
  }
}

static VALUE rb_gsl_matrix_complex_conjugate(VALUE obj)
{
  gsl_matrix_complex *cm = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, cm);
  gsl_matrix_complex_conjugate(cm);
  return obj;
}

static VALUE rb_gsl_matrix_complex_conjugate2(VALUE obj)
{
  gsl_matrix_complex *cm = NULL, *cmnew = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, cm);
  cmnew = gsl_matrix_complex_alloc(cm->size1, cm->size2);
  if (cmnew == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
  gsl_matrix_complex_conjugate2(cmnew, cm);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
}

static VALUE rb_gsl_matrix_complex_dagger(VALUE obj)
{
  gsl_matrix_complex *cm = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, cm);
  gsl_matrix_complex_conjugate(cm);
  gsl_matrix_complex_transpose(cm);
  return obj;
}

static VALUE rb_gsl_matrix_complex_dagger2(VALUE obj)
{
  gsl_matrix_complex *cm = NULL, *cmnew = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, cm);
  cmnew = gsl_matrix_complex_alloc(cm->size1, cm->size2);
  if (cmnew == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
  gsl_matrix_complex_conjugate2(cmnew, cm);
  gsl_matrix_complex_transpose(cmnew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
}

static VALUE rb_gsl_matrix_complex_trace(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  gsl_complex *trace = NULL;
  VALUE vtrace;
  size_t i;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  vtrace = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, trace);
  trace->dat[0] = 0.0; trace->dat[1] = 0.0;
  for (i = 0; i < m->size1; i++) *trace = gsl_complex_add(*trace, gsl_matrix_complex_get(m, i, i));
  return vtrace;
}

static VALUE rb_gsl_matrix_complex_each_row(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  gsl_vector_complex_view *vv;
  size_t i;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  for (i = 0; i < m->size1; i++) {
    vv = ALLOC(gsl_vector_complex_view);
    *vv = gsl_matrix_complex_row(m, i);
    rb_yield(Data_Wrap_Struct(cgsl_vector_complex_view, 0, free, vv));
  }
  return obj;
}

static VALUE rb_gsl_matrix_complex_each_col(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  gsl_vector_complex_view *vv;
  size_t i;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  for (i = 0; i < m->size2; i++) {
    vv = ALLOC(gsl_vector_complex_view);
    *vv = gsl_matrix_complex_column(m, i);
    rb_yield(Data_Wrap_Struct(cgsl_vector_complex_col_view, 0, free, vv));
  }
  return obj;
}

static VALUE rb_gsl_matrix_complex_size1(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  return INT2FIX(m->size1);
}
static VALUE rb_gsl_matrix_complex_size2(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  return INT2FIX(m->size2);
}
static VALUE rb_gsl_matrix_complex_shape(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  return rb_ary_new3(2, INT2FIX(m->size1), INT2FIX(m->size2));
}

static VALUE rb_gsl_matrix_complex_uplus(VALUE obj)
{
  return obj;
}

static VALUE rb_gsl_matrix_complex_uminus(VALUE obj)
{
  gsl_matrix_complex *m = NULL, *mnew;
  size_t i, j;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      gsl_matrix_complex_set(mnew, i, j, gsl_complex_negative(gsl_matrix_complex_get(m, i, j)));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mnew);
}

static VALUE rb_gsl_matrix_complex_XXX(VALUE obj, double (*f)(gsl_complex))
{
  gsl_matrix_complex *m;
  gsl_matrix *mnew;
  size_t i, j;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  mnew = gsl_matrix_alloc(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      gsl_matrix_set(mnew, i, j, (*f)(gsl_matrix_complex_get(m, i, j)));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
}

static VALUE rb_gsl_matrix_complex_arg(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX(obj, gsl_complex_arg);
}

static VALUE rb_gsl_matrix_complex_abs(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX(obj, gsl_complex_abs);
}

static VALUE rb_gsl_matrix_complex_abs2(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX(obj, gsl_complex_abs2);
}

static VALUE rb_gsl_matrix_complex_logabs(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX(obj, gsl_complex_logabs);
}

static VALUE rb_gsl_matrix_complex_XXX_complex(VALUE obj, 
                 gsl_complex (*f)(gsl_complex))
{
  gsl_matrix_complex *m, *mnew;
  size_t i, j;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      gsl_matrix_complex_set(mnew, i, j, (*f)(gsl_matrix_complex_get(m, i, j)));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mnew);
}

static VALUE rb_gsl_matrix_complex_sqrt(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_sqrt);
}

static VALUE rb_gsl_matrix_complex_exp(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_exp);
}

static VALUE rb_gsl_matrix_complex_log(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_log);
}

static VALUE rb_gsl_matrix_complex_log10(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_log10);
}

static VALUE rb_gsl_matrix_complex_sin(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_sin);
}

static VALUE rb_gsl_matrix_complex_cos(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_cos);
}

static VALUE rb_gsl_matrix_complex_tan(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_tan);
}

static VALUE rb_gsl_matrix_complex_sec(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_sec);
}

static VALUE rb_gsl_matrix_complex_csc(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_csc);
}

static VALUE rb_gsl_matrix_complex_cot(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_cot);
}

static VALUE rb_gsl_matrix_complex_arcsin(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arcsin);
}

static VALUE rb_gsl_matrix_complex_arccos(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arccos);
}

static VALUE rb_gsl_matrix_complex_arctan(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arctan);
}

static VALUE rb_gsl_matrix_complex_arcsec(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arcsec);
}

static VALUE rb_gsl_matrix_complex_arccsc(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arccsc);
}

static VALUE rb_gsl_matrix_complex_arccot(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arccot);
}

static VALUE rb_gsl_matrix_complex_sinh(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_sinh);
}

static VALUE rb_gsl_matrix_complex_cosh(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_cosh);
}

static VALUE rb_gsl_matrix_complex_tanh(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_tanh);
}

static VALUE rb_gsl_matrix_complex_sech(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_sech);
}

static VALUE rb_gsl_matrix_complex_csch(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_csch);
}
static VALUE rb_gsl_matrix_complex_coth(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_coth);
}

static VALUE rb_gsl_matrix_complex_arcsinh(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arcsinh);
}

static VALUE rb_gsl_matrix_complex_arccosh(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arccosh);
}

static VALUE rb_gsl_matrix_complex_arctanh(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arctanh);
}

static VALUE rb_gsl_matrix_complex_arcsech(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arcsech);
}

static VALUE rb_gsl_matrix_complex_arccsch(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arccsch);
}

static VALUE rb_gsl_matrix_complex_arccoth(VALUE obj)
{
  return rb_gsl_matrix_complex_XXX_complex(obj, gsl_complex_arccoth);
}

static VALUE rb_gsl_matrix_complex_indgen_bang(int argc, VALUE *argv[], VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  double start = 0, step = 1, x;
  size_t i, j;
  switch (argc) {
  case 0:
    break;
  case 1:
    start = NUM2DBL(argv[0]);
    break;
  case 2:
    start = NUM2DBL(argv[0]);
    step = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-2)", argc);
  }
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  x = start;
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      gsl_matrix_complex_set(m, i, j, gsl_complex_rect(x, 0));
      x += step;
    }
  }
  return obj;
}

static VALUE rb_gsl_matrix_complex_indgen(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL, *mnew;
  double start = 0, step = 1, x;
  size_t i, j;
  switch (argc) {
  case 0:
    break;
  case 1:
    start = NUM2DBL(argv[0]);
    break;
  case 2:
    start = NUM2DBL(argv[0]);
    step = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-2)", argc);
  }
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  mnew = gsl_matrix_complex_calloc(m->size1, m->size2);
  x = start;
  for (i = 0; i < mnew->size1; i++) {
    for (j = 0; j < mnew->size2; j++) {
      gsl_matrix_complex_set(mnew, i, j, gsl_complex_rect(x, 0));
      x += step;
    }
  }
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mnew);
}

static VALUE rb_gsl_matrix_complex_indgen_singleton(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *mnew;
  double start = 0, step = 1, x;
  size_t n1, n2, i, j;
  switch (argc) {
  case 2:
    n1 = (size_t) NUM2INT(argv[0]);
    n2 = (size_t) NUM2INT(argv[1]);
    break;
  case 3:
    n1 = (size_t) NUM2INT(argv[0]);
    n2 = (size_t) NUM2INT(argv[1]);
    start = NUM2DBL(argv[2]);
    break;
  case 4:
    n1 = (size_t) NUM2INT(argv[0]);
    n2 = (size_t) NUM2INT(argv[1]);
    start = NUM2DBL(argv[2]);
    step = NUM2DBL(argv[3]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-4)", argc);
  }
  mnew = gsl_matrix_complex_calloc(n1, n2);
  x = start;
  for (i = 0; i < mnew->size1; i++) {
    for (j = 0; j < mnew->size2; j++) {
      gsl_matrix_complex_set(mnew, i, j, gsl_complex_rect(x, 0));
      x += step;
    }
  }
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mnew);
}

// Starting with version 1.15, GSL provides a gsl_matrix_complex_equal
// function, but it only determines absolute equality (i.e. is has no epsilon
// argument).
static int gsl_matrix_complex_equal_eps(const gsl_matrix_complex *m1,
  const gsl_matrix_complex *m2, double eps)
{
  gsl_complex z1, z2;
  size_t i, j;
  if (m1->size1 != m2->size1) return 0;
  if (m1->size2 != m2->size2) return 0;  
  for (i = 0; i < m1->size1; i++) {
    for (j = 0; j < m1->size2; j++) {    
      z1 = gsl_matrix_complex_get(m1, i, j);
      z2 = gsl_matrix_complex_get(m2, i, j);
      if (!rbgsl_complex_equal(&z1, &z2, eps)) return 0;
    }
  }
  return 1;
}

static VALUE rb_gsl_matrix_complex_equal(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m1, *m2;
  double eps = 1e-8;
  int ret;
  switch (argc) {
  case 1:
    eps = 1e-8;
    break;
  case 2:
    eps = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 1 or 2)\n", argc);
  }
  Data_Get_Struct(obj, gsl_matrix_complex, m1);
  CHECK_MATRIX_COMPLEX(argv[0]);
  Data_Get_Struct(argv[0], gsl_matrix_complex, m2);
  ret = gsl_matrix_complex_equal_eps(m1, m2, eps);
  if (ret == 1) return Qtrue;
  else return Qfalse;
}

static VALUE rb_gsl_matrix_complex_not_equal(int argc, VALUE *argv, VALUE obj)
{
  VALUE ret;
  ret = rb_gsl_matrix_complex_equal(argc, argv, obj);
  if (ret == Qtrue) return Qfalse;
  else return Qtrue;
}

void Init_gsl_matrix_complex(VALUE module)
{
  //  rb_define_singleton_method(cgsl_matrix_complex, "new", rb_gsl_matrix_complex_new, 2);
  rb_define_singleton_method(cgsl_matrix_complex, "alloc", rb_gsl_matrix_complex_new, 2);
  rb_define_singleton_method(cgsl_matrix_complex, "[]", rb_gsl_matrix_complex_new, 2);
  rb_define_singleton_method(cgsl_matrix_complex, "calloc", rb_gsl_matrix_complex_new, 2);
  rb_define_singleton_method(cgsl_matrix_complex, "eye", rb_gsl_matrix_complex_eye, -1);
  rb_define_singleton_method(cgsl_matrix_complex, "diagonal", rb_gsl_matrix_complex_eye, -1);
  rb_define_singleton_method(cgsl_matrix_complex, "identity", rb_gsl_matrix_complex_identity, 1);
  rb_define_singleton_method(cgsl_matrix_complex, "unit", rb_gsl_matrix_complex_identity, 1);
  rb_define_singleton_method(cgsl_matrix_complex, "I", rb_gsl_matrix_complex_identity, 1);

  rb_define_method(cgsl_matrix_complex, "set", rb_gsl_matrix_complex_set, -1);
  rb_define_alias(cgsl_matrix_complex, "[]=", "set");

  rb_define_method(cgsl_matrix_complex, "set_row", rb_gsl_matrix_complex_set_row, -1);
  rb_define_method(cgsl_matrix_complex, "set_col", rb_gsl_matrix_complex_set_col, -1);

  rb_define_method(cgsl_matrix_complex, "get", rb_gsl_matrix_complex_get, -1);
  rb_define_alias(cgsl_matrix_complex, "[]", "get");
  rb_define_method(cgsl_matrix_complex, "ptr", rb_gsl_matrix_complex_ptr, 2);

  rb_define_method(cgsl_matrix_complex, "to_s", rb_gsl_matrix_complex_to_s, -1);
  rb_define_method(cgsl_matrix_complex, "fprintf", rb_gsl_matrix_complex_fprintf, -1);
  rb_define_method(cgsl_matrix_complex, "printf", rb_gsl_matrix_complex_printf, -1);
  rb_define_method(cgsl_matrix_complex, "print", rb_gsl_matrix_complex_print, 0);
  rb_define_method(cgsl_matrix_complex, "inspect", rb_gsl_matrix_complex_inspect, -1);
  rb_define_method(cgsl_matrix_complex, "fwrite", rb_gsl_matrix_complex_fwrite, 1);
  rb_define_method(cgsl_matrix_complex, "fread", rb_gsl_matrix_complex_fread, 1);
  rb_define_method(cgsl_matrix_complex, "fscanf", rb_gsl_matrix_complex_fscanf, 1);
  
  rb_define_singleton_method(cgsl_matrix_complex, "memcpy", rb_gsl_matrix_complex_memcpy, 2);
  rb_define_method(cgsl_matrix_complex, "clone", rb_gsl_matrix_complex_clone, 0);
  rb_define_alias(cgsl_matrix_complex, "duplicate", "clone");
  rb_define_alias(cgsl_matrix_complex, "dup", "clone");
  rb_define_method(cgsl_matrix_complex, "swap_rows", rb_gsl_matrix_complex_swap_rows, 2);
  rb_define_method(cgsl_matrix_complex, "swap_columns", rb_gsl_matrix_complex_swap_columns, 2);
  rb_define_method(cgsl_matrix_complex, "swap_rowcol", rb_gsl_matrix_complex_swap_rowcol, 2);
  
  rb_define_method(cgsl_matrix_complex, "transpose", rb_gsl_matrix_complex_transpose, 0);
  rb_define_method(cgsl_matrix_complex, "isnull", rb_gsl_matrix_complex_isnull, 0);
  
  rb_define_method(cgsl_matrix_complex, "add", rb_gsl_matrix_complex_add, 1);
  rb_define_alias(cgsl_matrix_complex, "add_constant", "add");
  rb_define_alias(cgsl_matrix_complex, "+", "add");
  rb_define_method(cgsl_matrix_complex, "sub", rb_gsl_matrix_complex_sub, 1);
  rb_define_alias(cgsl_matrix_complex, "-", "sub");
  rb_define_method(cgsl_matrix_complex, "mul_elements", rb_gsl_matrix_complex_mul_elements, 1);
  rb_define_method(cgsl_matrix_complex, "div_elements", rb_gsl_matrix_complex_div_elements, 1);
  rb_define_alias(cgsl_matrix_complex, "/", "div_elements");
  rb_define_method(cgsl_matrix_complex, "scale", rb_gsl_matrix_complex_scale, 1);
  rb_define_method(cgsl_matrix_complex, "scale!", rb_gsl_matrix_complex_scale_bang, 1);
  
  rb_define_method(cgsl_matrix_complex, "add_diagonal", rb_gsl_matrix_complex_add_diagonal, 1);
  
  rb_define_method(cgsl_matrix_complex, "set_zero", rb_gsl_matrix_complex_set_zero, 0);
  rb_define_method(cgsl_matrix_complex, "set_identity", rb_gsl_matrix_complex_set_identity, 0);
  rb_define_method(cgsl_matrix_complex, "set_all", rb_gsl_matrix_complex_set_all, 1);
  
  rb_define_method(cgsl_matrix_complex, "submatrix", rb_gsl_matrix_complex_submatrix, -1);
  rb_define_alias(cgsl_matrix_complex, "view", "submatrix");
  rb_define_method(cgsl_matrix_complex, "row", rb_gsl_matrix_complex_row, 1);
  /*  rb_define_alias(cgsl_matrix_complex, "[]", "row");*/
  rb_define_method(cgsl_matrix_complex, "column", rb_gsl_matrix_complex_column, 1);
  rb_define_alias(cgsl_matrix_complex, "col", "column");
  rb_define_method(cgsl_matrix_complex, "diagonal", rb_gsl_matrix_complex_diagonal, 0);
  rb_define_alias(cgsl_matrix_complex, "diag", "diagonal");
  rb_define_method(cgsl_matrix_complex, "set_diagonal", rb_gsl_matrix_complex_set_diagonal, 1);
  rb_define_method(cgsl_matrix_complex, "subdiagonal", rb_gsl_matrix_complex_subdiagonal, 1);
  rb_define_method(cgsl_matrix_complex, "superdiagonal", rb_gsl_matrix_complex_superdiagonal, 1);
  
  rb_define_method(cgsl_matrix_complex, "coerce", rb_gsl_matrix_complex_coerce, 1);
  
  rb_define_method(cgsl_matrix_complex, "mul", rb_gsl_matrix_complex_mul, 1);
  rb_define_alias(cgsl_matrix_complex, "*", "mul");
  rb_define_method(cgsl_matrix_complex, "mul!", rb_gsl_matrix_complex_mul2, 1);
  
  rb_define_method(cgsl_matrix_complex, "real", rb_gsl_matrix_complex_real, 0);
  rb_define_alias(cgsl_matrix_complex, "to_real", "real");
  rb_define_alias(cgsl_matrix_complex, "re", "real");
  rb_define_method(cgsl_matrix_complex, "imag", rb_gsl_matrix_complex_imag, 0);
  rb_define_alias(cgsl_matrix_complex, "im", "imag");

  /* 25.June.2004 */
  rb_define_method(cgsl_matrix_complex, "conjugate!", rb_gsl_matrix_complex_conjugate, 0);
  rb_define_alias(cgsl_matrix_complex, "conj!", "conjugate!");
  rb_define_method(cgsl_matrix_complex, "conjugate", rb_gsl_matrix_complex_conjugate2, 0);
  rb_define_alias(cgsl_matrix_complex, "conj", "conjugate");
  rb_define_method(cgsl_matrix_complex, "dagger!", rb_gsl_matrix_complex_dagger, 0);
  rb_define_method(cgsl_matrix_complex, "dagger", rb_gsl_matrix_complex_dagger2, 0);
  
  rb_define_method(cgsl_matrix_complex, "trace", rb_gsl_matrix_complex_trace, 0);
  rb_define_method(cgsl_matrix_complex, "each_row", rb_gsl_matrix_complex_each_row, 0);
  rb_define_method(cgsl_matrix_complex, "each_col", rb_gsl_matrix_complex_each_col, 0);
  rb_define_alias(cgsl_matrix_complex, "each_column", "each_col");
  rb_define_method(cgsl_matrix_complex, "collect", rb_gsl_matrix_complex_collect, 0);
  rb_define_method(cgsl_matrix_complex, "collect!", rb_gsl_matrix_complex_collect_bang, 0);
  rb_define_alias(cgsl_matrix_complex, "map", "collect");
  rb_define_alias(cgsl_matrix_complex, "map!", "collect!");
  
  rb_define_method(cgsl_matrix_complex, "to_a", rb_gsl_matrix_complex_to_a, 0);

  rb_define_method(cgsl_matrix_complex, "size1", rb_gsl_matrix_complex_size1, 0);
  rb_define_method(cgsl_matrix_complex, "size2", rb_gsl_matrix_complex_size2, 0);
  rb_define_method(cgsl_matrix_complex, "shape", rb_gsl_matrix_complex_shape, 0);
  rb_define_alias(cgsl_matrix_complex, "size", "shape");
  
  /*****/
  rb_define_method(cgsl_matrix_complex, "-@", rb_gsl_matrix_complex_uminus, 0);
  rb_define_method(cgsl_matrix_complex, "+@", rb_gsl_matrix_complex_uplus, 0);

  /****/
  rb_define_method(cgsl_matrix_complex, "arg", rb_gsl_matrix_complex_arg, 0);
  rb_define_alias(cgsl_matrix_complex, "angle", "arg");
  rb_define_alias(cgsl_matrix_complex, "phase", "arg");
  rb_define_method(cgsl_matrix_complex, "abs", rb_gsl_matrix_complex_abs, 0);
  rb_define_alias(cgsl_matrix_complex, "amp", "abs");
  rb_define_method(cgsl_matrix_complex, "abs2", rb_gsl_matrix_complex_abs2, 0);
  rb_define_method(cgsl_matrix_complex, "logabs", rb_gsl_matrix_complex_logabs, 0);

  rb_define_method(cgsl_matrix_complex, "sqrt", rb_gsl_matrix_complex_sqrt, 0);
  rb_define_method(cgsl_matrix_complex, "exp", rb_gsl_matrix_complex_exp, 0);
  rb_define_method(cgsl_matrix_complex, "log", rb_gsl_matrix_complex_log, 0);
  rb_define_method(cgsl_matrix_complex, "log10", rb_gsl_matrix_complex_log10, 0);

  rb_define_method(cgsl_matrix_complex, "sin", rb_gsl_matrix_complex_sin, 0);
  rb_define_method(cgsl_matrix_complex, "cos", rb_gsl_matrix_complex_cos, 0);
  rb_define_method(cgsl_matrix_complex, "tan", rb_gsl_matrix_complex_tan, 0);
  rb_define_method(cgsl_matrix_complex, "sec", rb_gsl_matrix_complex_sec, 0);
  rb_define_method(cgsl_matrix_complex, "csc", rb_gsl_matrix_complex_csc, 0);
  rb_define_method(cgsl_matrix_complex, "cot", rb_gsl_matrix_complex_cot, 0);

  rb_define_method(cgsl_matrix_complex, "arcsin", rb_gsl_matrix_complex_arcsin, 0);
  rb_define_method(cgsl_matrix_complex, "arccos", rb_gsl_matrix_complex_arccos, 0);
  rb_define_method(cgsl_matrix_complex, "arctan", rb_gsl_matrix_complex_arctan, 0);
  rb_define_method(cgsl_matrix_complex, "arcsec", rb_gsl_matrix_complex_arcsec, 0);
  rb_define_method(cgsl_matrix_complex, "arccsc", rb_gsl_matrix_complex_arccsc, 0);
  rb_define_method(cgsl_matrix_complex, "arccot", rb_gsl_matrix_complex_arccot, 0);

  rb_define_method(cgsl_matrix_complex, "sinh", rb_gsl_matrix_complex_sinh, 0);
  rb_define_method(cgsl_matrix_complex, "cosh", rb_gsl_matrix_complex_cosh, 0);
  rb_define_method(cgsl_matrix_complex, "tanh", rb_gsl_matrix_complex_tanh, 0);
  rb_define_method(cgsl_matrix_complex, "sech", rb_gsl_matrix_complex_sech, 0);
  rb_define_method(cgsl_matrix_complex, "csch", rb_gsl_matrix_complex_csch, 0);
  rb_define_method(cgsl_matrix_complex, "coth", rb_gsl_matrix_complex_coth, 0);

  rb_define_method(cgsl_matrix_complex, "arcsinh", rb_gsl_matrix_complex_arcsinh, 0);
  rb_define_method(cgsl_matrix_complex, "arccosh", rb_gsl_matrix_complex_arccosh, 0);
  rb_define_method(cgsl_matrix_complex, "arctanh", rb_gsl_matrix_complex_arctanh, 0);
  rb_define_method(cgsl_matrix_complex, "arcsech", rb_gsl_matrix_complex_arcsech, 0);
  rb_define_method(cgsl_matrix_complex, "arccsch", rb_gsl_matrix_complex_arccsch, 0);
  rb_define_method(cgsl_matrix_complex, "arccoth", rb_gsl_matrix_complex_arccoth, 0);

  rb_define_method(cgsl_matrix_complex, "indgen", rb_gsl_matrix_complex_indgen, -1);
  rb_define_method(cgsl_matrix_complex, "indgen!", rb_gsl_matrix_complex_indgen_bang, -1);
  rb_define_singleton_method(cgsl_matrix_complex, "indgen", rb_gsl_matrix_complex_indgen_singleton, -1);
  
  rb_define_method(cgsl_matrix_complex, "equal?", rb_gsl_matrix_complex_equal, -1);
  rb_define_alias(cgsl_matrix_complex, "==", "equal?");
  rb_define_method(cgsl_matrix_complex, "not_equal?", rb_gsl_matrix_complex_not_equal, -1);
  rb_define_alias(cgsl_matrix_complex, "!=", "not_equal?");    
}
