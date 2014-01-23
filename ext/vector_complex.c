/*
  vector_complex.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/
#include "rb_gsl_config.h"
#include "rb_gsl_array.h"
#include "rb_gsl_complex.h"

EXTERN VALUE cgsl_complex;
static VALUE rb_gsl_vector_complex_inner_product(int argc, VALUE *argv, VALUE obj);
static VALUE rb_gsl_vector_complex_product_to_m(int argc, VALUE *argv, VALUE obj);

// From vector_source.c
void get_range_beg_en_n(VALUE range, double *beg, double *en, size_t *n, int *step);
//
// From vector_source.c
void parse_subvector_args(int argc, VALUE *argv, size_t size,
    size_t *offset, size_t *stride, size_t *n);

// From complex.c
gsl_complex rb_gsl_obj_to_gsl_complex(VALUE obj, gsl_complex *z);

static VALUE rb_gsl_vector_complex_new(int argc, VALUE *argv, VALUE klass)
{
  VALUE tmp;
  gsl_vector_complex *v = NULL;
  gsl_vector *x, *y;
  gsl_complex z, *z2 = NULL;
  size_t n, i;
  switch (argc) {
  case 1:
    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      n = FIX2INT(argv[0]);
      v = gsl_vector_complex_calloc(n);
      if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
      break;
    case T_ARRAY: 
      n = RARRAY_LEN(argv[0]);
      v = gsl_vector_complex_alloc(n);
      if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
      for (i = 0; i < n; i++) {
	z2 = &z;
	tmp = rb_ary_entry(argv[0], i);
	if (TYPE(tmp) == T_ARRAY) {
	  GSL_SET_REAL(z2, NUM2DBL(rb_ary_entry(tmp, 0)));
	  GSL_SET_IMAG(z2, NUM2DBL(rb_ary_entry(tmp, 1)));
	} else if (COMPLEX_P(tmp)) {
	  Data_Get_Struct(tmp, gsl_complex, z2);
	} else {
	  rb_raise(rb_eTypeError, 
		   "wrong argument type %s (Array or Complex expected)", 
		   rb_class2name(CLASS_OF(tmp)));
	  
	}
	gsl_vector_complex_set(v, i, *z2);
      }
      break;
    default:
      rb_raise(rb_eTypeError, 
	       "wrong argument type %s", rb_class2name(CLASS_OF(argv[0])));
      break;
    }
    break;
  default:
    if (argc == 2 && (VECTOR_P(argv[0]) && VECTOR_P(argv[1]))) {
      Data_Get_Struct(argv[0], gsl_vector, x);
      Data_Get_Struct(argv[1], gsl_vector, y);
      n = GSL_MIN_INT(x->size, y->size);
      v = gsl_vector_complex_alloc(n);
      for (i = 0; i < n; i++) {
	z.dat[0] = gsl_vector_get(x, i);
	z.dat[1] = gsl_vector_get(y, i);
	gsl_vector_complex_set(v, i, z);
      }
      break;
    }
    n = argc;
    v = gsl_vector_complex_alloc(n);
    if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
    for (i = 0; i < n; i++) {
      if (TYPE(argv[i]) == T_ARRAY) {
	GSL_SET_REAL(&z, NUM2DBL(rb_ary_entry(argv[i], 0)));
	GSL_SET_IMAG(&z, NUM2DBL(rb_ary_entry(argv[i], 1)));
	z2 = &z;
      } else if (COMPLEX_P(argv[i])) {
	Data_Get_Struct(argv[i], gsl_complex, z2);
      } else {
	rb_raise(rb_eTypeError, 
		 "wrong argument type %s (Array or Complex expected)", 
		 rb_class2name(CLASS_OF(argv[i])));
      }
      gsl_vector_complex_set(v, i, *z2);
    }
    break;
  }
  return Data_Wrap_Struct(klass, 0, gsl_vector_complex_free, v);
}

static VALUE rb_gsl_vector_complex_row_new(int argc, VALUE *argv, VALUE klass)
{
  return rb_gsl_vector_complex_new(argc, argv, klass);
}

static VALUE rb_gsl_vector_complex_calloc(VALUE klass, VALUE nn)
{
  gsl_vector_complex *vc = NULL;
  CHECK_FIXNUM(nn);
  vc = gsl_vector_complex_calloc(FIX2INT(nn));
  if (vc == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
  return Data_Wrap_Struct(klass, 0, gsl_vector_complex_free, vc);
}

static VALUE rb_gsl_vector_complex_size(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  return INT2FIX(v->size);
}

static VALUE rb_gsl_vector_complex_stride(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  return INT2FIX(v->stride);
}

static VALUE rb_gsl_vector_complex_owner(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  return INT2FIX(v->owner);
}

static VALUE rb_gsl_vector_complex_ptr(VALUE obj, VALUE i)
{
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  return Data_Wrap_Struct(cgsl_complex, 0, NULL, gsl_vector_complex_ptr(v, FIX2INT(i)));
}

// TODO return Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj),..." where appropriate
static VALUE rb_gsl_vector_complex_subvector(int argc, VALUE *argv, VALUE obj);
static VALUE rb_gsl_vector_complex_get(int argc, VALUE *argv, VALUE obj)
{
  VALUE retval = Qnil;
  gsl_vector_complex *v = NULL, *vnew;
  gsl_complex *c = NULL;
  gsl_index *p;
  int i, k;
  size_t index, j;
  // If argc is not 1 or argv[0] is a Range
  if( argc != 1 || rb_obj_is_kind_of(argv[0], rb_cRange)) {
    // Treat as call to subvector
    retval = rb_gsl_vector_complex_subvector(argc, argv, obj);
  } else {
    Data_Get_Struct(obj, gsl_vector_complex, v);

    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      CHECK_FIXNUM(argv[0]);
      i = FIX2INT(argv[0]);
      if (i < 0) index = (size_t) (v->size + i);
      else index = (size_t) i;
      c = ALLOC(gsl_complex);
      *c = gsl_vector_complex_get(v, index);
      retval = Data_Wrap_Struct(cgsl_complex, 0, free, c);
      break;
    case T_ARRAY:
      vnew = gsl_vector_complex_alloc(RARRAY_LEN(argv[0]));
      for (j = 0; j < vnew->size; j++) {
	i = FIX2INT(rb_ary_entry(argv[0], j));
	if (i < 0) i = v->size + i;
	gsl_vector_complex_set(vnew, j, gsl_vector_complex_get(v, i));
      }
      retval = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
      break;
    default:
      if (PERMUTATION_P(argv[0])) {
        Data_Get_Struct(argv[0], gsl_index, p);
        vnew = gsl_vector_complex_alloc(p->size);
        for (j = 0; j < p->size; j++) {
          k = p->data[j];
          if (k < 0) k = p->size + j;
          gsl_vector_complex_set(vnew, j, gsl_vector_complex_get(v, k));
        }
        retval = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
      } else {
        // TODO Support Vector::Int (and even Vector?)
        rb_raise(rb_eTypeError, "wrong argument type %s (Array, Range, GSL::Permutation, or Fixnum expected)", rb_class2name(CLASS_OF(argv[0])));
      }
      break;
    }
  }
  return retval;
}

static VALUE rb_gsl_vector_complex_set_all(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_complex tmp;

  if (argc < 1) rb_raise(rb_eArgError, "wrong number of arguments");

  Data_Get_Struct(obj, gsl_vector_complex, v);

  switch (argc) {
  case 1:
    tmp = rb_gsl_obj_to_gsl_complex(argv[0], NULL);
    break;
  case 2:
    GSL_SET_COMPLEX(&tmp, NUM2DBL(argv[0]), NUM2DBL(argv[1]));
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }

  gsl_vector_complex_set_all(v, tmp);

  return obj;
}

void rb_gsl_vector_complex_set_subvector(int argc, VALUE *argv, gsl_vector_complex *v, VALUE other)
{
  gsl_vector_complex *vother;
  gsl_vector_complex_view vv;
  gsl_complex tmp;
  int step;
  size_t i, offset, stride, n, nother;
  double beg, end;

  // assignment to v.subvector(...)
  parse_subvector_args(argc, argv, v->size, &offset, &stride, &n);
  vv = gsl_vector_complex_subvector_with_stride(v, offset, stride, n);
  if(rb_obj_is_kind_of(other, cgsl_vector_complex)) {
    Data_Get_Struct(other, gsl_vector_complex, vother);
    if(n != vother->size) {
      rb_raise(rb_eRangeError, "lengths do not match (%d != %d)", (int) n, (int) vother->size);
    }
    // TODO Change to gsl_vector_complex_memmove if/when GSL has such a
    // function because gsl_vector_memcpy does not handle overlapping regions
    // (e.g.  Views) well.
    gsl_vector_complex_memcpy(&vv.vector, vother);
  } else if(rb_obj_is_kind_of(other, rb_cArray)) {
    // TODO Support other forms of Array contents as well
    if((int) n != RARRAY_LEN(other)) {
      rb_raise(rb_eRangeError, "lengths do not match (%d != %d)", (int) n, (int) RARRAY_LEN(other));
    }
    for(i = 0; i < n; i++) {
      tmp = rb_gsl_obj_to_gsl_complex(rb_ary_entry(other, i), NULL);
      gsl_vector_complex_set(&vv.vector, i, tmp);
    }
  } else if(rb_obj_is_kind_of(other, rb_cRange)) {
    get_range_beg_en_n(other, &beg, &end, &nother, &step);
    if(n != nother) {
      rb_raise(rb_eRangeError, "lengths do not match (%d != %d)", (int) n, (int) nother);
    }
    GSL_SET_IMAG(&tmp, 0.0);
    for(i = 0; i < n; i++) {
      GSL_SET_REAL(&tmp, beg);
      gsl_vector_complex_set(&vv.vector, i, tmp);
      beg += step;
    }
  } else {
    tmp = rb_gsl_obj_to_gsl_complex(argv[1], NULL);
    gsl_vector_complex_set_all(&vv.vector, tmp);
  }
}

static VALUE rb_gsl_vector_complex_set(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_complex tmp;
  VALUE other;
  int ii;

  if(argc < 1 || argc > 4) {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1-4)", argc);
  }

  Data_Get_Struct(obj, gsl_vector_complex, v);
  other = argv[argc-1];

  if(argc == 1) {
    // // If assigning from another vector
    if(VECTOR_P(other) || VECTOR_COMPLEX_P(other)) {
      // treat as assignment to v.subvector(...)
      rb_gsl_vector_complex_set_subvector(argc-1, argv, v, other);
    } else {
      // treat as set_all
      rb_gsl_vector_complex_set_all(argc, argv, obj);
    }
  } else if(argc == 2 && TYPE(argv[0]) == T_FIXNUM) {
    // v[i] = x
    ii = FIX2INT(argv[0]);
    if(ii < 0) ii += v->size;
    // Get/make GSL::Complex from second arg
    tmp = gsl_vector_complex_get(v, ii);
    rb_gsl_obj_to_gsl_complex(argv[1], &tmp);
    gsl_vector_complex_set(v, (size_t)ii, tmp);
  } else {
    // assignment to v.subvector(...)
    rb_gsl_vector_complex_set_subvector(argc-1, argv, v, other);
  }

  return obj;
}

static VALUE rb_gsl_vector_complex_each(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  VALUE vz;
  gsl_complex * zp;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  for (i = 0; i < v->size; i++) {
    vz = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, zp);
    *zp = gsl_vector_complex_get(v, i);
    rb_yield(vz);
  }
  return obj;
}

static VALUE rb_gsl_vector_complex_reverse_each(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  VALUE vz;
  gsl_complex * zp;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  for (i = v->size-1;; i--) {
    vz = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, zp);
    *zp = gsl_vector_complex_get(v, i);
    rb_yield(vz);
    if (i == 0) break;
  }
  return obj;
}

static VALUE rb_gsl_vector_complex_each_index(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  for (i = 0; i < v->size; i++) rb_yield(INT2FIX(i));
  return obj;
}

static VALUE rb_gsl_vector_complex_reverse_each_index(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  for (i = v->size-1;; i--) {
    rb_yield(INT2FIX(i));
    if (i == 0) break;
  }
  return obj;
}

static void rb_gsl_vector_complex_collect_native(gsl_vector_complex *src, gsl_vector_complex *dst)
{
  VALUE vz;
  gsl_complex * zp;
  size_t i;
  for (i = 0; i < src->size; i++) {
    vz = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, zp);
    *zp = gsl_vector_complex_get(src, i);
    vz = rb_yield(vz);
    CHECK_COMPLEX(vz);
    Data_Get_Struct(vz, gsl_complex, zp);
    gsl_vector_complex_set(dst, i, *zp);
  }
}

static VALUE rb_gsl_vector_complex_collect(VALUE obj)
{
  gsl_vector_complex *v = NULL, *vnew;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vnew = gsl_vector_complex_alloc(v->size);
  rb_gsl_vector_complex_collect_native(v, vnew);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
}

static VALUE rb_gsl_vector_complex_collect_bang(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  rb_gsl_vector_complex_collect_native(v, v);
  return obj;
}

static VALUE rb_gsl_vector_complex_set_zero(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  gsl_vector_complex_set_zero(v);
  return obj;
}

static VALUE rb_gsl_vector_complex_set_basis(VALUE obj, VALUE ii)
{
  gsl_vector_complex *v = NULL;
  CHECK_FIXNUM(ii);
  Data_Get_Struct(obj, gsl_vector_complex, v);
  gsl_vector_complex_set_basis(v, FIX2INT(ii));
  return obj;
}

static VALUE rb_gsl_vector_complex_to_s(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  char buf[64];
  size_t i;
  VALUE str;
  gsl_complex * z;

  Data_Get_Struct(obj, gsl_vector_complex, v);
  if (v->size == 0) return rb_str_new2("[ ]");
  str = rb_str_new2("[ ");
  if (VECTOR_COMPLEX_COL_P(obj)) {
    for (i = 0; i < v->size; i++) {
      if (i != 0) {
	rb_str_cat(str, "  ", 2);
      }
      z = GSL_COMPLEX_AT(v, i);
      sprintf(buf, "[%4.3e %4.3e]", GSL_REAL(*z), GSL_IMAG(*z));
      if (i != v->size-1) strcat(buf, "\n");
      rb_str_cat(str, buf, strlen(buf));
      if (i >= 10 && i != v->size-1) {
	rb_str_cat(str, "  ...", 5);
	break;
      }
    }
  } else {
    z = GSL_COMPLEX_AT(v, 0);
    sprintf(buf, "[%4.3e %4.3e]", GSL_REAL(*z), GSL_IMAG(*z));
    rb_str_cat(str, buf, strlen(buf));
    for (i = 1; i < v->size; i++) {
      z = GSL_COMPLEX_AT(v, i);
      sprintf(buf, " [%4.3e %4.3e]", GSL_REAL(*z), GSL_IMAG(*z));
      rb_str_cat(str, buf, strlen(buf));
      if (i >= 10 && i != v->size-1) {
	rb_str_cat(str, " ...", 4);
	break;
      }
    }
  }
  rb_str_cat(str, " ]", 2);
  return str;
}

static VALUE rb_gsl_vector_complex_inspect(VALUE obj)
{
  VALUE str;
  char buf[128];
  gsl_vector_complex *v;

  Data_Get_Struct(obj, gsl_vector_complex, v);
  sprintf(buf, "#<%s[%lu]:%#lx>\n", rb_class2name(CLASS_OF(obj)), v->size, NUM2ULONG(rb_obj_id(obj)));
  str = rb_str_new2(buf);
  return rb_str_concat(str, rb_gsl_vector_complex_to_s(obj));
}

/*static VALUE rb_gsl_vector_complex_fprintf(VALUE obj, VALUE io, VALUE format)*/
static VALUE rb_gsl_vector_complex_fprintf(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  if (argc != 1 && argc != 2) {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
  Data_Get_Struct(obj, gsl_vector_complex, h);
  fp = rb_gsl_open_writefile(argv[0], &flag);
  if (argc == 2) {
    Check_Type(argv[1], T_STRING);
    status = gsl_vector_complex_fprintf(fp, h, STR2CSTR(argv[1]));
  } else {
    status = gsl_vector_complex_fprintf(fp, h, "%4.3e");
  }
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

/*static VALUE rb_gsl_vector_complex_printf(VALUE obj, VALUE format)*/
static VALUE rb_gsl_vector_complex_printf(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *h = NULL;
  int status;
  Data_Get_Struct(obj, gsl_vector_complex, h);
  if (argc == 1) {
    Check_Type(argv[0], T_STRING);
    status = gsl_vector_complex_fprintf(stdout, h, STR2CSTR(argv[0]));
  } else {
    status = gsl_vector_complex_fprintf(stdout, h, "%4.3e");
  }
  return INT2FIX(status);
}
/*
static VALUE rb_gsl_vector_complex_print(VALUE obj);
static VALUE rb_gsl_vector_complex_inspect(VALUE obj)
{
  gsl_vector_complex *h = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, h);
  printf("%s\n", rb_class2name(CLASS_OF(obj)));
  rb_gsl_vector_complex_print(obj);
  return Qtrue;
}
*/
static VALUE rb_gsl_vector_complex_print(VALUE obj)
{
  gsl_vector_complex *h = NULL;
  gsl_complex *z = NULL;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, h);
  printf("[ ");
  if (VECTOR_COMPLEX_COL_P(obj)) {
    printf("\n");
    for (i = 0; i < h->size; i++) {
      z = GSL_COMPLEX_AT(h, i);
      printf("  [%4.3e %4.3e]\n", GSL_REAL(*z), GSL_IMAG(*z));
    }
  } else {
    for (i = 0; i < h->size; i++) {
      z = GSL_COMPLEX_AT(h, i);
      printf("[%4.3e %4.3e] ", GSL_REAL(*z), GSL_IMAG(*z));
    }
  }
  printf("]\n");
  return obj;
}

static VALUE rb_gsl_vector_complex_fwrite(VALUE obj, VALUE io)
{
  gsl_vector_complex *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_vector_complex, h);
  fp = rb_gsl_open_writefile(io, &flag);
  status = gsl_vector_complex_fwrite(fp, h);
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE rb_gsl_vector_complex_fread(VALUE obj, VALUE io)
{
  gsl_vector_complex *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_vector_complex, h);
  f = rb_gsl_open_writefile(io, &flag);
  status = gsl_vector_complex_fread(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_vector_complex_fscanf(VALUE obj, VALUE io)
{
  gsl_vector_complex *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_vector_complex, h);
  f = rb_gsl_open_writefile(io, &flag);
  status = gsl_vector_complex_fscanf(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_vector_complex_real(VALUE obj)
{
  gsl_vector_complex *c = NULL;
  gsl_vector_view *vv = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, c);
  vv = gsl_vector_view_alloc();
  *vv = gsl_vector_complex_real(c);
  if (VECTOR_COMPLEX_ROW_P(obj))
    return Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free, vv);
  else
    return Data_Wrap_Struct(cgsl_vector_col_view, 0, gsl_vector_view_free, vv);
}

static VALUE rb_gsl_vector_complex_imag(VALUE obj)
{
  gsl_vector_complex *c = NULL;
  gsl_vector_view *vv = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, c);
  vv = gsl_vector_view_alloc();
  *vv = gsl_vector_complex_imag(c);
  if (VECTOR_COMPLEX_ROW_P(obj))
    return Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free, vv);
  else
    return Data_Wrap_Struct(cgsl_vector_col_view, 0, gsl_vector_view_free, vv);
}

static VALUE rb_gsl_vector_complex_set_real(VALUE obj, VALUE val)
{
  gsl_vector_complex *v = NULL;
  gsl_vector_view vv;
  double d = NUM2DBL(val);
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vv = gsl_vector_complex_real(v);
  gsl_vector_set_all(&vv.vector, d);
  return obj;
}

static VALUE rb_gsl_vector_complex_set_imag(VALUE obj, VALUE val)
{
  gsl_vector_complex *v = NULL;
  gsl_vector_view vv;
  double d = NUM2DBL(val);
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vv = gsl_vector_complex_imag(v);
  gsl_vector_set_all(&vv.vector, d);
  return obj;
}

static VALUE rb_gsl_vector_complex_conj(VALUE obj)
{
  size_t i;
  gsl_vector_complex *vin = NULL;
  gsl_vector_complex *vout = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, vin);
  vout = gsl_vector_complex_alloc(vin->size);
  for(i=0; i<vin->size; i++) {
    gsl_vector_complex_set(vout, i,
        gsl_complex_conjugate(
          gsl_vector_complex_get(vin, i)));
  }
  return Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj), 0, gsl_vector_complex_free, vout);
}

static VALUE rb_gsl_vector_complex_conj_bang(VALUE obj)
{
  size_t i;
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  for(i=0; i<v->size; i++) {
    gsl_vector_complex_set(v, i,
        gsl_complex_conjugate(
          gsl_vector_complex_get(v, i)));
  }
  return obj;
}

static VALUE rb_gsl_vector_complex_to_a(VALUE obj)
{
  gsl_vector_complex *c = NULL;
  gsl_complex *z = NULL;
  size_t i, j;
  VALUE ary;
  Data_Get_Struct(obj, gsl_vector_complex, c);
  ary = rb_ary_new2(c->size*2);
  for (i = 0, j = 0; i < c->size; i++, j+=2) {
    z = GSL_COMPLEX_AT(c, i);
    rb_ary_store(ary, j, rb_float_new(GSL_REAL(*z)));
    rb_ary_store(ary, j+1, rb_float_new(GSL_IMAG(*z)));
  }
  return ary;
}

static VALUE rb_gsl_vector_complex_to_a2(VALUE obj)
{
  gsl_vector_complex *c = NULL;
  gsl_complex *znew = NULL, *z = NULL;
  size_t i;
  VALUE ary, vz;
  Data_Get_Struct(obj, gsl_vector_complex, c);
  ary = rb_ary_new2(c->size);
  for (i = 0; i < c->size; i++) {
    z = GSL_COMPLEX_AT(c, i);
    znew = make_complex(z->dat[0], z->dat[1]);
    vz = Data_Wrap_Struct(cgsl_complex, 0, free, znew);
    rb_ary_store(ary, i, vz);
  }
  return ary;
}

gsl_vector_complex_view* gsl_vector_complex_view_alloc()
{
  gsl_vector_complex_view *vv = NULL;
  vv = ALLOC(gsl_vector_complex_view);
  if (vv == NULL) rb_raise(rb_eRuntimeError, "malloc failed");
  return vv;
}

void gsl_vector_complex_view_free(gsl_vector_view * vv)
{
  free((gsl_vector_complex_view *) vv);
}

static VALUE rb_gsl_vector_complex_subvector(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_vector_complex_view *vv = NULL;
  size_t offset, stride, n;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  parse_subvector_args(argc, argv, v->size, &offset, &stride, &n);
  vv = gsl_vector_complex_view_alloc();
  *vv = gsl_vector_complex_subvector_with_stride(v, offset, stride, n);
  if (VECTOR_COMPLEX_ROW_P(obj))
    return Data_Wrap_Struct(cgsl_vector_complex_view, 0, gsl_vector_complex_view_free, vv);
  else
    return Data_Wrap_Struct(cgsl_vector_complex_col_view, 0, gsl_vector_complex_view_free, vv);
}

static VALUE rb_gsl_vector_complex_subvector_with_stride(VALUE obj, VALUE o, VALUE s, VALUE nn)
{
  gsl_vector_complex *v = NULL;
  gsl_vector_complex_view *vv = NULL;
  int offset;
  CHECK_FIXNUM(o); CHECK_FIXNUM(nn); CHECK_FIXNUM(s);
  offset = NUM2INT(o);
  Data_Get_Struct(obj, gsl_vector_complex, v);
  if(offset < 0) {
    offset += v->size;
  }
  vv = gsl_vector_complex_view_alloc();
  *vv = gsl_vector_complex_subvector_with_stride(v, (size_t)offset, FIX2INT(s), FIX2INT(nn));
  if (VECTOR_COMPLEX_ROW_P(obj))
    return Data_Wrap_Struct(cgsl_vector_complex_view, 0, gsl_vector_complex_view_free, vv);
  else
    return Data_Wrap_Struct(cgsl_vector_complex_col_view, 0, gsl_vector_complex_view_free, vv);
}

/* singleton */
static VALUE rb_gsl_vector_complex_memcpy(VALUE obj, VALUE dst, VALUE src)
{
  gsl_vector_complex *v = NULL, *dest = NULL;
  CHECK_VECTOR_COMPLEX(dst);
  CHECK_VECTOR_COMPLEX(src);
  Data_Get_Struct(dst, gsl_vector_complex, dest);
  Data_Get_Struct(src, gsl_vector_complex, v);
  gsl_vector_complex_memcpy(dest, v);
  return dst;
}

static VALUE rb_gsl_vector_complex_clone(VALUE obj)
{
  gsl_vector_complex *v = NULL, *vnew = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vnew = gsl_vector_complex_alloc(v->size);
  if (vnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
  gsl_vector_complex_memcpy(vnew, v);
  if (VECTOR_COMPLEX_ROW_P(obj))
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
  else
    return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, vnew);
}

static VALUE rb_gsl_vector_complex_reverse(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  gsl_vector_complex_reverse(v);
  return obj;
}

static VALUE rb_gsl_vector_complex_reverse2(VALUE obj)
{
  gsl_vector_complex *v = NULL, *vnew = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vnew = gsl_vector_complex_alloc(v->size);
  gsl_vector_complex_memcpy(vnew, v);
  gsl_vector_complex_reverse(vnew);
  if (VECTOR_COMPLEX_ROW_P(obj)) 
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
  else
    return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, vnew);
}

static VALUE rb_gsl_vector_complex_swap_elements(VALUE obj, VALUE i, VALUE j)
{
  gsl_vector_complex *v = NULL;
  CHECK_FIXNUM(i); CHECK_FIXNUM(j);
  Data_Get_Struct(obj, gsl_vector_complex, v);
  gsl_vector_complex_swap_elements(v, FIX2INT(i), FIX2INT(j));
  return obj;
}

static VALUE rb_gsl_vector_complex_fftshift_bang(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_complex tmp;
  size_t i, n;

  Data_Get_Struct(obj, gsl_vector_complex, v);
  n = v->size;
  if(n & 1) {
    // length is odd
    tmp = gsl_vector_complex_get(v,0);
    for(i = 0; i < n/2; i++) {
      gsl_vector_complex_set(v, i, gsl_vector_complex_get(v, i+n/2+1));
      gsl_vector_complex_set(v, i+n/2+1, gsl_vector_complex_get(v, i+1));
    }
    gsl_vector_complex_set(v, n/2, tmp);
  } else {
    // length is even
    for(i = 0; i < n/2; i++) {
      gsl_vector_complex_swap_elements(v, i, i+n/2);
    }
  }

  return obj;
}

static VALUE rb_gsl_vector_complex_fftshift(VALUE obj)
{
  gsl_vector_complex *v, *vnew;
  gsl_vector_complex_view vv, vvnew;
  size_t n;

  Data_Get_Struct(obj, gsl_vector_complex, v);
  n = v->size;
  vnew = gsl_vector_complex_alloc(n);
  // Copy low to high
  vv = gsl_vector_complex_subvector(v, 0, (n+1)/2);
  vvnew = gsl_vector_complex_subvector(vnew, n/2, (n+1)/2);
  gsl_vector_complex_memcpy(&vvnew.vector, &vv.vector);
  // Copy high to low
  vv = gsl_vector_complex_subvector(v, (n+1)/2, n/2);
  vvnew = gsl_vector_complex_subvector(vnew, 0, n/2);
  gsl_vector_complex_memcpy(&vvnew.vector, &vv.vector);

  return Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj), 0, gsl_vector_complex_free, vnew);  
}

static VALUE rb_gsl_vector_complex_ifftshift_bang(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_complex tmp;
  size_t i, n;

  Data_Get_Struct(obj, gsl_vector_complex, v);
  n = v->size;
  if(n & 1) {
    // length is odd
    tmp = gsl_vector_complex_get(v,n/2);
    for(i = n/2; i > 0; i--) {
      gsl_vector_complex_set(v, i, gsl_vector_complex_get(v, i+n/2));
      gsl_vector_complex_set(v, i+n/2, gsl_vector_complex_get(v, i-1));
    }
    gsl_vector_complex_set(v, 0, tmp);
  } else {
    // length is even
    for(i = 0; i < n/2; i++) {
      gsl_vector_complex_swap_elements(v, i, i+n/2);
    }
  }

  return obj;
}

static VALUE rb_gsl_vector_complex_ifftshift(VALUE obj)
{
  return rb_gsl_vector_complex_ifftshift_bang(rb_gsl_vector_complex_clone(obj));
  /*gsl_vector_complex *v, *vnew;
  gsl_vector_complex_view vv, vvnew;
  size_t n;

  Data_Get_Struct(obj, gsl_vector_complex, v);
  n = v->size;
  vnew = gsl_vector_complex_alloc(n);
  // Copy high to low
  vv = gsl_vector_complex_subvector(vnew, n/2, (n+1)/2);
  vvnew = gsl_vector_complex_subvector(v, 0, (n+1)/2);
  gsl_vector_complex_memcpy(&vvnew.vector, &vv.vector);
  // Copy low to high
  vv = gsl_vector_complex_subvector(vnew, 0, n/2);
  vvnew = gsl_vector_complex_subvector(v, (n+1)/2, n/2);
  gsl_vector_complex_memcpy(&vvnew.vector, &vv.vector);

  return Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj), 0, gsl_vector_complex_free, vnew);  */
}

static VALUE rb_gsl_vector_complex_isnull(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  if (gsl_vector_complex_isnull(v)) return Qtrue;
  else return Qfalse;
}

static VALUE rb_gsl_vector_complex_matrix_view(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_matrix_complex_view *mv = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  switch (argc) {
  case 2:
    mv = gsl_matrix_complex_view_alloc();
    *mv = gsl_matrix_complex_view_vector(v, FIX2INT(argv[0]), FIX2INT(argv[1]));
    break;
  case 3:
    mv = gsl_matrix_complex_view_alloc();
    *mv = gsl_matrix_complex_view_vector_with_tda(v, FIX2INT(argv[0]), FIX2INT(argv[1]), 
				 FIX2INT(argv[2]));
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  return Data_Wrap_Struct(cgsl_matrix_complex_view, 0, gsl_matrix_complex_view_free, mv);
}

static VALUE rb_gsl_vector_complex_matrix_view_with_tda(VALUE obj, VALUE nn1, VALUE nn2,
						VALUE tda)
{
  gsl_vector_complex *v = NULL;
  gsl_matrix_complex_view *mv = NULL;
  CHECK_FIXNUM(nn1); CHECK_FIXNUM(nn2); CHECK_FIXNUM(tda);
  Data_Get_Struct(obj, gsl_vector_complex, v);
  mv = gsl_matrix_complex_view_alloc();
  if (mv == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
  *mv = gsl_matrix_complex_view_vector_with_tda(v, FIX2INT(nn1), FIX2INT(nn2), FIX2INT(tda));
  return Data_Wrap_Struct(cgsl_matrix_complex_view, 0, gsl_matrix_complex_view_free, mv);
}

static VALUE rb_gsl_vector_complex_trans(VALUE obj)
{
  gsl_vector_complex *v = NULL, *vnew = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vnew = make_vector_complex_clone(v);
  if (VECTOR_COMPLEX_ROW_P(obj)) 
    return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, vnew);
  else if (VECTOR_COMPLEX_COL_P(obj)) 
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
  else {
    rb_raise(rb_eTypeError, "wrong type");
  }
}

static VALUE rb_gsl_vector_complex_trans2(VALUE obj)
{
  if (CLASS_OF(obj) == cgsl_vector_complex) 
    RBGSL_SET_CLASS(obj, cgsl_vector_complex_col);
  else if (CLASS_OF(obj) == cgsl_vector_complex_col) 
    RBGSL_SET_CLASS(obj, cgsl_vector_complex);
  else {
    rb_raise(rb_eRuntimeError, "method trans! for %s is forbidden",
	     rb_class2name(CLASS_OF(obj)));
  }
  return obj;
}

static VALUE rb_gsl_vector_complex_to_real(VALUE obj)
{
  gsl_vector_complex *cv = NULL;
  gsl_vector *v = NULL;
  gsl_complex z;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, cv);
  v = gsl_vector_alloc(cv->size);
  if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
  for (i = 0; i < cv->size; i++) {
    z = gsl_vector_complex_get(cv, i);
    gsl_vector_set(v, i, GSL_REAL(z));
  }
  if (VECTOR_COMPLEX_ROW_P(obj))
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
  else
    return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, v);
}

enum {
  GSL_VECTOR_COMPLEX_ADD,
  GSL_VECTOR_COMPLEX_SUB,
  GSL_VECTOR_COMPLEX_MUL,
  GSL_VECTOR_COMPLEX_DIV,
  GSL_VECTOR_COMPLEX_ADD_BANG,
  GSL_VECTOR_COMPLEX_SUB_BANG,
  GSL_VECTOR_COMPLEX_MUL_BANG,
  GSL_VECTOR_COMPLEX_DIV_BANG,
};

static VALUE rb_gsl_vector_complex_arithmetics(int flag, VALUE obj, VALUE bb);

static VALUE rb_gsl_vector_complex_arithmetics(int flag, VALUE obj, VALUE bb) 
{
  gsl_vector *b = NULL;
  gsl_vector_complex *cv = NULL, *cvnew = NULL, *cb = NULL;
  gsl_complex *c = NULL, z;
  Data_Get_Struct(obj, gsl_vector_complex, cv);
  switch (flag) {
  case GSL_VECTOR_COMPLEX_ADD:
  case GSL_VECTOR_COMPLEX_SUB:
  case GSL_VECTOR_COMPLEX_MUL:
  case GSL_VECTOR_COMPLEX_DIV:
    cvnew = make_vector_complex_clone(cv);
    obj = Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj), 0, gsl_vector_complex_free, cvnew);
    break;
  case GSL_VECTOR_COMPLEX_ADD_BANG:
  case GSL_VECTOR_COMPLEX_SUB_BANG:
  case GSL_VECTOR_COMPLEX_MUL_BANG:
  case GSL_VECTOR_COMPLEX_DIV_BANG:
    cvnew = cv;
    break;
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  switch (TYPE(bb)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    z = gsl_complex_rect(NUM2DBL(bb), 0.0);
    switch (flag) {
    case GSL_VECTOR_COMPLEX_ADD:
    case GSL_VECTOR_COMPLEX_ADD_BANG:
      gsl_vector_complex_add_constant(cvnew, z);
      break;
    case GSL_VECTOR_COMPLEX_SUB:
    case GSL_VECTOR_COMPLEX_SUB_BANG:
      gsl_vector_complex_add_constant(cvnew, gsl_complex_negative(z));
      break;
    case GSL_VECTOR_COMPLEX_MUL:
    case GSL_VECTOR_COMPLEX_MUL_BANG:
      gsl_vector_complex_scale(cvnew, z);
      break;
    case GSL_VECTOR_COMPLEX_DIV:
    case GSL_VECTOR_COMPLEX_DIV_BANG:
      gsl_vector_complex_scale(cvnew, gsl_complex_inverse(z));
      break;
    }
    break;
  default:
    if (VECTOR_P(bb)) {
      Data_Get_Struct(bb, gsl_vector, b);
      cb = vector_to_complex(b);
      switch (flag) {
      case GSL_VECTOR_COMPLEX_ADD:
      case GSL_VECTOR_COMPLEX_ADD_BANG:
	gsl_vector_complex_add(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_SUB:
      case GSL_VECTOR_COMPLEX_SUB_BANG:
	gsl_vector_complex_sub(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_MUL:
      case GSL_VECTOR_COMPLEX_MUL_BANG:
	gsl_vector_complex_mul(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_DIV:
      case GSL_VECTOR_COMPLEX_DIV_BANG:
	gsl_vector_complex_div(cvnew, cb);
	break;
      }
      gsl_vector_complex_free(cb);
    } else if (VECTOR_COMPLEX_P(bb)) {
      Data_Get_Struct(bb, gsl_vector_complex, cb);
      switch (flag) {
      case GSL_VECTOR_COMPLEX_ADD:
      case GSL_VECTOR_COMPLEX_ADD_BANG:
	gsl_vector_complex_add(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_SUB:
      case GSL_VECTOR_COMPLEX_SUB_BANG:
	gsl_vector_complex_sub(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_MUL:
      case GSL_VECTOR_COMPLEX_MUL_BANG:
	gsl_vector_complex_mul(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_DIV:
      case GSL_VECTOR_COMPLEX_DIV_BANG:
	gsl_vector_complex_div(cvnew, cb);
	break;
      }
    } else if (COMPLEX_P(bb)) {
      Data_Get_Struct(bb, gsl_complex, c);
      switch (flag) {
      case GSL_VECTOR_COMPLEX_ADD:
      case GSL_VECTOR_COMPLEX_ADD_BANG:
	gsl_vector_complex_add_constant(cvnew, *c);
	break;
      case GSL_VECTOR_COMPLEX_SUB:
      case GSL_VECTOR_COMPLEX_SUB_BANG:
	gsl_vector_complex_add_constant(cvnew, gsl_complex_negative(*c));
	break;
      case GSL_VECTOR_COMPLEX_MUL:
      case GSL_VECTOR_COMPLEX_MUL_BANG:
	gsl_vector_complex_scale(cvnew, *c);
	break;
      case GSL_VECTOR_COMPLEX_DIV:
      case GSL_VECTOR_COMPLEX_DIV_BANG:
	gsl_vector_complex_scale(cvnew, gsl_complex_inverse(*c));
	break;
      }
    } else {
      rb_raise(rb_eTypeError, "wrong type argument %s", rb_class2name(CLASS_OF(bb)));
    }
    break;
  }
  return obj;
}

static VALUE rb_gsl_vector_complex_add(VALUE obj, VALUE bb)
{
  return rb_gsl_vector_complex_arithmetics(GSL_VECTOR_COMPLEX_ADD, obj, bb);
}

static VALUE rb_gsl_vector_complex_sub(VALUE obj, VALUE bb)
{
  return rb_gsl_vector_complex_arithmetics(GSL_VECTOR_COMPLEX_SUB, obj, bb);
}

static VALUE rb_gsl_vector_complex_mul(VALUE obj, VALUE bb)
{
  VALUE argv[2];
  if (VECTOR_COMPLEX_ROW_P(obj) && VECTOR_COMPLEX_COL_P(bb)) {
    argv[0] = obj;
    argv[1] = bb;
    return rb_gsl_vector_complex_inner_product(2, argv, CLASS_OF(obj));
  }
  if (VECTOR_COMPLEX_COL_P(obj) && VECTOR_COMPLEX_ROW_P(bb)) {
    argv[0] = obj;
    argv[1] = bb;
    return rb_gsl_vector_complex_product_to_m(2, argv, CLASS_OF(obj));
  }
  return rb_gsl_vector_complex_arithmetics(GSL_VECTOR_COMPLEX_MUL, obj, bb);
}

static VALUE rb_gsl_vector_complex_div(VALUE obj, VALUE bb)
{
  return rb_gsl_vector_complex_arithmetics(GSL_VECTOR_COMPLEX_DIV, obj, bb);
}

static VALUE rb_gsl_vector_complex_add_bang(VALUE obj, VALUE bb)
{
  return rb_gsl_vector_complex_arithmetics(GSL_VECTOR_COMPLEX_ADD_BANG, obj, bb);
}

static VALUE rb_gsl_vector_complex_sub_bang(VALUE obj, VALUE bb)
{
  return rb_gsl_vector_complex_arithmetics(GSL_VECTOR_COMPLEX_SUB_BANG, obj, bb);
}

static VALUE rb_gsl_vector_complex_mul_bang(VALUE obj, VALUE bb)
{
  return rb_gsl_vector_complex_arithmetics(GSL_VECTOR_COMPLEX_MUL_BANG, obj, bb);
}

static VALUE rb_gsl_vector_complex_div_bang(VALUE obj, VALUE bb)
{
  return rb_gsl_vector_complex_arithmetics(GSL_VECTOR_COMPLEX_DIV_BANG, obj, bb);
}

static VALUE rb_gsl_vector_complex_coerce(VALUE obj, VALUE other)
{
  gsl_vector_complex *cv = NULL, *cb = NULL;
  gsl_complex z;
  VALUE vv;
  Data_Get_Struct(obj, gsl_vector_complex, cv);
  switch (TYPE(other)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    z = gsl_complex_rect(NUM2DBL(other), 0.0);
    cb = gsl_vector_complex_alloc(cv->size);
    if (cb == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
    gsl_vector_complex_set_all(cb, z);
    vv = Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj), 0, gsl_vector_complex_free, cb);
    return rb_ary_new3(2, vv, obj);
    break;
  default:
    rb_raise(rb_eTypeError, "GSL::Vector::Complex, operation not defined");
    break;
  }
}

/* 2.Aug.2004 */
static VALUE rb_gsl_vector_complex_inner_product(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL, *v2 = NULL;
  gsl_complex prod, a, b, *z = NULL;
  size_t i;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
			    argc);
    if (!VECTOR_COMPLEX_ROW_P(argv[0]))
      rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Vector::Complex expected)",
	       rb_class2name(CLASS_OF(argv[0])));
    if (!VECTOR_COMPLEX_COL_P(argv[1]))
      rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Vector::Complex::Col expected)",
	       rb_class2name(CLASS_OF(argv[1])));
    Data_Get_Struct(argv[0], gsl_vector_complex, v);
    Data_Get_Struct(argv[1], gsl_vector_complex, v2);
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    if (!VECTOR_COMPLEX_COL_P(argv[0]))
      rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Vector::Complex::Col expected)",
	       rb_class2name(CLASS_OF(argv[0])));
    Data_Get_Struct(obj, gsl_vector_complex, v);
    Data_Get_Struct(argv[0], gsl_vector_complex, v2);
    break;
  }
  if (v->size != v2->size) rb_raise(rb_eRangeError, "vector lengths are different.");
  prod = gsl_complex_rect(0.0, 0.0);
  for (i = 0; i < v->size; i++) {
    a = gsl_vector_complex_get(v, i);
    b = gsl_vector_complex_get(v2, i);
    prod = gsl_complex_add(prod, gsl_complex_mul(a, b));
  }
  z = ALLOC(gsl_complex);
  *z = prod;
  return Data_Wrap_Struct(cgsl_complex, 0, free, z);
}

static VALUE rb_gsl_vector_complex_product_to_m(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL, *v2 = NULL;
  gsl_matrix_complex *m = NULL;
  gsl_complex a, b;
  size_t i, j;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
			    argc);
    if (!VECTOR_COMPLEX_COL_P(argv[0]))
      rb_raise(rb_eTypeError, 
	       "wrong argument type %s (GSL::Vector::Complex::Col expected)",
	       rb_class2name(CLASS_OF(argv[0])));
    if (!VECTOR_COMPLEX_ROW_P(argv[1]))
      rb_raise(rb_eTypeError, 
	       "wrong argument type %s (GSL::Vector::Complex expected)",
	       rb_class2name(CLASS_OF(argv[1])));
    Data_Get_Struct(argv[0], gsl_vector_complex, v);
    Data_Get_Struct(argv[1], gsl_vector_complex, v2);
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    if (!VECTOR_COMPLEX_ROW_P(argv[0]))
      rb_raise(rb_eTypeError, 
	       "wrong argument type %s (GSL::Vector::Complex expected)",
	       rb_class2name(CLASS_OF(argv[0])));
    Data_Get_Struct(obj, gsl_vector_complex, v);
    Data_Get_Struct(argv[0], gsl_vector_complex, v2);
    break;
  }
  m = gsl_matrix_complex_alloc(v->size, v2->size);
  for (i = 0; i < v->size; i++) {
    for (j = 0; j < v2->size; j++) {
      a = gsl_vector_complex_get(v, i);
      b = gsl_vector_complex_get(v2, j);
      gsl_matrix_complex_set(m, i, j, gsl_complex_mul(a, b));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, m);
}

static VALUE rb_gsl_vector_complex_uplus(VALUE obj)
{
  return obj;
}

static VALUE rb_gsl_vector_complex_uminus(VALUE obj)
{
  gsl_vector_complex *v = NULL, *vnew;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vnew = gsl_vector_complex_alloc(v->size);
  for (i = 0; i < v->size; i++) {
    gsl_vector_complex_set(vnew, i, gsl_complex_negative(gsl_vector_complex_get(v, i)));
  }
  return Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj), 0, gsl_vector_complex_free, vnew);
}

/*****/
static VALUE rb_gsl_vector_complex_XXX(VALUE obj, double (*f)(gsl_complex))
{
  gsl_vector_complex *m;
  gsl_vector *v;
  gsl_complex c;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, m);
  v = gsl_vector_alloc(m->size);
  for (i = 0; i < m->size; i++) {
    c = gsl_vector_complex_get(m, i);
    gsl_vector_set(v, i, (*f)(c));
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

static VALUE rb_gsl_vector_complex_XXXz(VALUE obj, gsl_complex (*f)(gsl_complex))
{
  gsl_vector_complex *m, *v;
  gsl_complex c;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, m);
  v = gsl_vector_complex_alloc(m->size);
  for (i = 0; i < m->size; i++) {
    c = gsl_vector_complex_get(m, i);
    gsl_vector_complex_set(v, i, (*f)(c));
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, v);
}

/* In-place version of rb_gsl_vector_complex_XXXz */
static VALUE rb_gsl_vector_complex_XXXz_bang(VALUE obj, gsl_complex (*f)(gsl_complex))
{
  gsl_vector_complex *v;
  gsl_complex c;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  for (i = 0; i < v->size; i++) {
    c = gsl_vector_complex_get(v, i);
    gsl_vector_complex_set(v, i, (*f)(c));
  }
  return obj;
}

static VALUE rb_gsl_vector_complex_XXXz2(VALUE obj, VALUE a, 
					 gsl_complex (*f)(gsl_complex, gsl_complex))
{
  gsl_vector_complex *m, *v;
  gsl_complex c, *z;
  size_t i;
  CHECK_COMPLEX(a);
  Data_Get_Struct(obj, gsl_vector_complex, m);
  Data_Get_Struct(a, gsl_complex, z);
  v = gsl_vector_complex_alloc(m->size);
  for (i = 0; i < m->size; i++) {
    c = gsl_vector_complex_get(m, i);
    gsl_vector_complex_set(v, i, (*f)(c, *z));
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, v);
}

static VALUE rb_gsl_vector_complex_XXXz2_bang(VALUE obj, VALUE a, 
					 gsl_complex (*f)(gsl_complex, gsl_complex))
{
  gsl_vector_complex *v;
  gsl_complex c, *z;
  size_t i;
  CHECK_COMPLEX(a);
  Data_Get_Struct(obj, gsl_vector_complex, v);
  Data_Get_Struct(a, gsl_complex, z);
  for (i = 0; i < v->size; i++) {
    c = gsl_vector_complex_get(v, i);
    gsl_vector_complex_set(v, i, (*f)(c, *z));
  }
  return obj;
}

static VALUE rb_gsl_vector_complex_abs2(VALUE obj)
{
  return rb_gsl_vector_complex_XXX(obj, gsl_complex_abs2);
}

static VALUE rb_gsl_vector_complex_abs(VALUE obj)
{
  return rb_gsl_vector_complex_XXX(obj, gsl_complex_abs);
}

static VALUE rb_gsl_vector_complex_logabs(VALUE obj)
{
  return rb_gsl_vector_complex_XXX(obj, gsl_complex_logabs);
}

static VALUE rb_gsl_vector_complex_arg(VALUE obj)
{
  return rb_gsl_vector_complex_XXX(obj, gsl_complex_arg);
}

static VALUE rb_gsl_vector_complex_sqrt(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_sqrt);
}

static VALUE rb_gsl_vector_complex_sqrt_bang(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz_bang(obj, gsl_complex_sqrt);
}

static VALUE rb_gsl_vector_complex_exp(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_exp);
}

static VALUE rb_gsl_vector_complex_exp_bang(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz_bang(obj, gsl_complex_exp);
}

static VALUE rb_gsl_vector_complex_pow(VALUE obj, VALUE a)
{
  return rb_gsl_vector_complex_XXXz2(obj, a, gsl_complex_pow);
}

static VALUE rb_gsl_vector_complex_pow_bang(VALUE obj, VALUE a)
{
  return rb_gsl_vector_complex_XXXz2_bang(obj, a, gsl_complex_pow);
}

static VALUE rb_gsl_vector_complex_log(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_log);
}

static VALUE rb_gsl_vector_complex_log_bang(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz_bang(obj, gsl_complex_log);
}

static VALUE rb_gsl_vector_complex_log10(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_log10);
}

static VALUE rb_gsl_vector_complex_log10_bang(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz_bang(obj, gsl_complex_log10);
}

static VALUE rb_gsl_vector_complex_log_b(VALUE obj, VALUE a)
{
  return rb_gsl_vector_complex_XXXz2(obj, a, gsl_complex_log_b);
}

static VALUE rb_gsl_vector_complex_log_b_bang(VALUE obj, VALUE a)
{
  return rb_gsl_vector_complex_XXXz2_bang(obj, a, gsl_complex_log_b);
}

/* gsl_vector_complex_sum */
static gsl_complex rb_gsl_vector_complex_sum_gsl(gsl_vector_complex * v)
{
  size_t i;
  gsl_complex z = gsl_complex_rect(0.0,0.0);

  for(i=0; i<v->size; i++) {
    z = gsl_complex_add(z, gsl_vector_complex_get(v,i));
  }
  return z;
}

/* gsl_vector_complex_mean */
static gsl_complex rb_gsl_vector_complex_mean_gsl(gsl_vector_complex * v)
{
  gsl_complex z = rb_gsl_vector_complex_sum_gsl(v);
  return gsl_complex_div_real(z, (double)v->size);
}

/* gsl_vector_complex_tss_m */
static double rb_gsl_vector_complex_tss_m_gsl(gsl_vector_complex * v, gsl_complex mean)
{
  size_t i;
  double tss = 0.0;

  for(i=0; i<v->size; i++) {
    tss += gsl_complex_abs2(gsl_complex_sub(gsl_vector_complex_get(v,i), mean));
  }
  return tss;
}

/* gsl_vector_complex_variance_m */
static double rb_gsl_vector_complex_variance_m_gsl(gsl_vector_complex * v, gsl_complex mean)
{
  double tss = rb_gsl_vector_complex_tss_m_gsl(v, mean);
  return tss / (double)(v->size - 1);
}

/* gsl_vector_complex_variance_with_fixed_mean */
static double rb_gsl_vector_complex_variance_fm_gsl(gsl_vector_complex * v, gsl_complex mean)
{
  double tss = rb_gsl_vector_complex_tss_m_gsl(v, mean);
  return tss / (double)(v->size);
}

/* gsl_vector_complex_sd_m */
static double rb_gsl_vector_complex_sd_m_gsl(gsl_vector_complex * v, gsl_complex mean)
{
  double var = rb_gsl_vector_complex_variance_m_gsl(v, mean);
  return sqrt(var);
}

/* gsl_vector_complex_sd_with_fixed_mean */
static double rb_gsl_vector_complex_sd_fm_gsl(gsl_vector_complex * v, gsl_complex mean)
{
  double var = rb_gsl_vector_complex_variance_fm_gsl(v, mean);
  return sqrt(var);
}

/* gsl_vector_complex_tss */
static double rb_gsl_vector_complex_tss_gsl(gsl_vector_complex * v)
{
  gsl_complex mean = rb_gsl_vector_complex_mean_gsl(v);
  return rb_gsl_vector_complex_tss_m_gsl(v, mean);
}

/* gsl_vector_complex_variance */
static double rb_gsl_vector_complex_variance_gsl(gsl_vector_complex * v)
{
  double tss = rb_gsl_vector_complex_tss_gsl(v);
  return tss / (double)(v->size - 1);
}

/* gsl_vector_complex_sd */
static double rb_gsl_vector_complex_sd_gsl(gsl_vector_complex * v)
{
  double var = rb_gsl_vector_complex_variance_gsl(v);
  return sqrt(var);
}

/* Wrapper around stats funcs with prototype like
 * "gsl_complex func(gsl_vector_complex *v)"
 * (e.g. sum and mean)
 */
static VALUE rb_gsl_vector_complex_z_stats_v(VALUE obj,
    gsl_complex (*func)(gsl_vector_complex*))
{
  gsl_vector_complex * v;
  gsl_complex * zp;
  VALUE zv;

  CHECK_VECTOR_COMPLEX(obj);
  Data_Get_Struct(obj, gsl_vector_complex, v);

  zv = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, zp); 
  *zp = func(v);

  return zv;
}

/* Wrapper around stats funcs with prototype like
 * "double func(gsl_vector_complex *v)"
 * (e.g. tss, variance, sd)
 */
static VALUE rb_gsl_vector_complex_d_stats_v(VALUE obj,
    double (*func)(gsl_vector_complex*))
{
  gsl_vector_complex * v;
  double d;

  CHECK_VECTOR_COMPLEX(obj);
  Data_Get_Struct(obj, gsl_vector_complex, v);

  d = func(v);

  return rb_float_new(d);
}

/* Wrapper around stats funcs with prototype like
 * "double func(gsl_vector_complex *v, gsl_complex z)"
 * (e.g. tss_m, variance_m, sd_m, variance_fm, sd_fm)
 */
static VALUE rb_gsl_vector_complex_d_stats_v_z(VALUE obj, VALUE arg,
    double (*func)(gsl_vector_complex*, gsl_complex))
{
  gsl_vector_complex * v;
  gsl_complex z;
  gsl_complex * zp;
  double d;

  CHECK_VECTOR_COMPLEX(obj);
  Data_Get_Struct(obj, gsl_vector_complex, v);


  switch (TYPE(arg)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    z = gsl_complex_rect(NUM2DBL(arg), 0.0);
    zp = &z;
    break;
  default:
    CHECK_COMPLEX(arg);
    Data_Get_Struct(arg, gsl_complex, zp);
    break;
  }

  d = func(v,*zp);

  return rb_float_new(d);
}

static VALUE rb_gsl_vector_complex_sum(VALUE obj)
{
  return rb_gsl_vector_complex_z_stats_v(obj, rb_gsl_vector_complex_sum_gsl);
}

static VALUE rb_gsl_vector_complex_mean(VALUE obj)
{
  return rb_gsl_vector_complex_z_stats_v(obj, rb_gsl_vector_complex_mean_gsl);
}

static VALUE rb_gsl_vector_complex_tss_m(VALUE obj, VALUE arg)
{
  return rb_gsl_vector_complex_d_stats_v_z(obj, arg, rb_gsl_vector_complex_tss_m_gsl);
}

static VALUE rb_gsl_vector_complex_variance_m(VALUE obj, VALUE arg)
{
  return rb_gsl_vector_complex_d_stats_v_z(obj, arg, rb_gsl_vector_complex_variance_m_gsl);
}

static VALUE rb_gsl_vector_complex_variance_fm(VALUE obj, VALUE arg)
{
  return rb_gsl_vector_complex_d_stats_v_z(obj, arg, rb_gsl_vector_complex_variance_fm_gsl);
}

static VALUE rb_gsl_vector_complex_sd_m(VALUE obj, VALUE arg)
{
  return rb_gsl_vector_complex_d_stats_v_z(obj, arg, rb_gsl_vector_complex_sd_m_gsl);
}

static VALUE rb_gsl_vector_complex_sd_fm(VALUE obj, VALUE arg)
{
  return rb_gsl_vector_complex_d_stats_v_z(obj, arg, rb_gsl_vector_complex_sd_fm_gsl);
}

static VALUE rb_gsl_vector_complex_tss(VALUE obj)
{
  return rb_gsl_vector_complex_d_stats_v(obj, rb_gsl_vector_complex_tss_gsl);
}

static VALUE rb_gsl_vector_complex_variance(VALUE obj)
{
  return rb_gsl_vector_complex_d_stats_v(obj, rb_gsl_vector_complex_variance_gsl);
}

static VALUE rb_gsl_vector_complex_sd(VALUE obj)
{
  return rb_gsl_vector_complex_d_stats_v(obj, rb_gsl_vector_complex_sd_gsl);
}

/*****/

static VALUE rb_gsl_vector_complex_sin(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_sin);
}

static VALUE rb_gsl_vector_complex_cos(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_cos);
}

static VALUE rb_gsl_vector_complex_tan(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_tan);
}

static VALUE rb_gsl_vector_complex_sec(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_sec);
}

static VALUE rb_gsl_vector_complex_csc(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_csc);
}
static VALUE rb_gsl_vector_complex_cot(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_cot);
}

static VALUE rb_gsl_vector_complex_arcsin(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arcsin);
}

static VALUE rb_gsl_vector_complex_arccos(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arccos);
}

static VALUE rb_gsl_vector_complex_arctan(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arctan);
}

static VALUE rb_gsl_vector_complex_arcsec(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arcsec);
}

static VALUE rb_gsl_vector_complex_arccsc(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arccsc);
}

static VALUE rb_gsl_vector_complex_arccot(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arccot);
}

static VALUE rb_gsl_vector_complex_sinh(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_sinh);
}

static VALUE rb_gsl_vector_complex_cosh(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_cosh);
}

static VALUE rb_gsl_vector_complex_tanh(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_tanh);
}

static VALUE rb_gsl_vector_complex_sech(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_sech);
}

static VALUE rb_gsl_vector_complex_csch(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_csch);
}

static VALUE rb_gsl_vector_complex_coth(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_coth);
}

static VALUE rb_gsl_vector_complex_arcsinh(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arcsinh);
}

static VALUE rb_gsl_vector_complex_arccosh(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arccosh);
}

static VALUE rb_gsl_vector_complex_arctanh(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arctanh);
}

static VALUE rb_gsl_vector_complex_arcsech(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arcsech);
}

static VALUE rb_gsl_vector_complex_arccsch(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arccsch);
}

static VALUE rb_gsl_vector_complex_arccoth(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_arccoth);
}

static VALUE rb_gsl_vector_complex_concat(VALUE obj, VALUE other)
{
  gsl_vector_complex *v = NULL, *v2 = NULL, *vnew = NULL;
  gsl_vector_complex_view vv;
  gsl_complex tmp;
  VALUE x;
  double beg, end;
  int step;
  size_t i, size2;

  Data_Get_Struct(obj, gsl_vector_complex, v);

  switch(TYPE(other)) {
    case T_FIXNUM:
    case T_BIGNUM:
    case T_FLOAT:
      vnew = gsl_vector_complex_alloc(v->size + 1);
      vv = gsl_vector_complex_subvector(vnew, 0, v->size);
      gsl_vector_complex_memcpy(&vv.vector, v);
      gsl_vector_complex_set(vnew, v->size, rb_gsl_obj_to_gsl_complex(other, NULL));
      break;

    case T_ARRAY:
      size2 = RARRAY_LEN(other);
      vnew = gsl_vector_complex_alloc(v->size + size2);
      vv = gsl_vector_complex_subvector(vnew, 0, v->size);
      gsl_vector_complex_memcpy(&vv.vector, v);
      for (i = 0; i < size2; i++) {
        x = rb_ary_entry(other, i);
        gsl_vector_complex_set(vnew, v->size + i, rb_gsl_obj_to_gsl_complex(x, NULL));
      }
      break;

    default:
      if(rb_obj_is_kind_of(other, cgsl_complex)) {
        vnew = gsl_vector_complex_alloc(v->size + 1);
        vv = gsl_vector_complex_subvector(vnew, 0, v->size);
        gsl_vector_complex_memcpy(&vv.vector, v);
        gsl_vector_complex_set(vnew, v->size, rb_gsl_obj_to_gsl_complex(other, NULL));
      } else if(rb_obj_is_kind_of(other, rb_cRange)) {
        get_range_beg_en_n(other, &beg, &end, &size2, &step);
        vnew = gsl_vector_complex_alloc(v->size + size2);
        vv = gsl_vector_complex_subvector(vnew, 0, v->size);
        gsl_vector_complex_memcpy(&vv.vector, v);
        GSL_SET_COMPLEX(&tmp, beg, 0.0);
        for (i = 0; i < size2; i++) {
          gsl_vector_complex_set(vnew, v->size + i, tmp);
          GSL_SET_REAL(&tmp, GSL_REAL(tmp) + step);
        }
      } else if (rb_obj_is_kind_of(other, cgsl_vector_complex)) {
        Data_Get_Struct(other, gsl_vector_complex, v2);
        size2 = v2->size;
        vnew = gsl_vector_complex_alloc(v->size + size2);
        vv = gsl_vector_complex_subvector(vnew, 0, v->size);
        gsl_vector_complex_memcpy(&vv.vector, v);
        vv = gsl_vector_complex_subvector(vnew, v->size, size2);
        gsl_vector_complex_memcpy(&vv.vector, v2);
      } else {
        rb_raise(rb_eTypeError, "wrong argument type %s (Array, Numeric, Range, GSL::Complex, or %s expected)",
            rb_class2name(CLASS_OF(other)), rb_class2name(cgsl_vector_complex));
      }
      break;
  }

  return Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj), 0, gsl_vector_complex_free, vnew);  
}

static VALUE rb_gsl_vector_complex_block(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  return Data_Wrap_Struct(cgsl_block_complex, 0, NULL, v->block);
}

static VALUE rb_gsl_vector_complex_indgen_bang(int argc, VALUE *argv[], VALUE obj)
{
  gsl_vector_complex *v = NULL;
  double start = 0.0, step = 1.0, x;
  size_t i;
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
  Data_Get_Struct(obj, gsl_vector_complex, v);
  for (i = 0, x = start; i < v->size; i++, x += step) {
    gsl_vector_complex_set(v, i, gsl_complex_rect(x,0));
  }
  return obj;
}

static VALUE rb_gsl_vector_complex_indgen(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL, *vnew;
  double start = 0, step = 1, x;
  size_t i;
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
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vnew = gsl_vector_complex_calloc(v->size);
  for (i = 0, x = start; i < vnew->size; i++, x += step) {
    gsl_vector_complex_set(vnew, i, gsl_complex_rect(x,0));
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
}

static VALUE rb_gsl_vector_complex_indgen_singleton(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *vnew;
  double start = 0, step = 1, x;
  size_t n, i;
  switch (argc) {
  case 1:
    n = (size_t) FIX2INT(argv[0]);
    break;
  case 2:
    n = (size_t) FIX2INT(argv[0]);
    start = NUM2DBL(argv[1]);
    break;
  case 3:
    n = (size_t) FIX2INT(argv[0]);
    start = NUM2DBL(argv[1]);
    step = NUM2DBL(argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-3)",argc);
  }
  vnew = gsl_vector_complex_calloc(n);
  for (i = 0, x = start; i < vnew->size; i++, x += step) {
    gsl_vector_complex_set(vnew, i, gsl_complex_rect(x,0));
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
}

static VALUE rb_gsl_vector_complex_phasor_singleton(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *vnew;
  double start, step, theta;
  size_t n, i;
  switch (argc) {
  case 1:
    n = (size_t) FIX2INT(argv[0]);
    start = 0;
    step  = 2 * M_PI / n;
    break;
  case 2:
    n = (size_t) FIX2INT(argv[0]);
    start = NUM2DBL(argv[1]);
    step  = 2 * M_PI / n;
    break;
  case 3:
    n = (size_t) FIX2INT(argv[0]);
    start = NUM2DBL(argv[1]);
    step  = NUM2DBL(argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-3)",argc);
  }
  vnew = gsl_vector_complex_alloc(n);
  for (i = 0, theta = start; i < vnew->size; i++, theta += step) {
    gsl_vector_complex_set(vnew, i, gsl_complex_polar(1.0,theta));
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
}

static VALUE rb_gsl_vector_complex_zip(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v0, **vp, *vnew;
  VALUE ary;
  size_t i, j;
  int argc2;
  VALUE *argv2;
  gsl_complex zzero = gsl_complex_rect(0, 0);
  if (VECTOR_COMPLEX_P(obj)) {
    Data_Get_Struct(obj, gsl_vector_complex, v0);
    argc2 = argc;
    argv2 = argv;
  } else {
    if (argc < 1) rb_raise(rb_eArgError, "Too few arguments.");
    Data_Get_Struct(argv[0], gsl_vector_complex, v0);    
    argc2 = argc - 1;
    argv2 = argv + 1;
  }
  for (i = 0; (int) i < argc2; i++) {
    CHECK_VECTOR_COMPLEX(argv2[i]);
  }
  vp = (gsl_vector_complex**) malloc(sizeof(gsl_vector_complex**));
  for (i = 0; (int) i < argc2; i++) {
    Data_Get_Struct(argv2[i], gsl_vector_complex, vp[i]);
  }
  ary = rb_ary_new2(v0->size);
  for (i = 0; i < v0->size; i++) {
    vnew = gsl_vector_complex_alloc(argc2 + 1);
    gsl_vector_complex_set(vnew, 0, gsl_vector_complex_get(v0, i));
    for (j = 0; (int) j < argc2; j++) {
      if (i < vp[j]->size) {
	gsl_vector_complex_set(vnew, j+1, gsl_vector_complex_get(vp[j], i));
      } else {
	gsl_vector_complex_set(vnew, j+1, zzero);
      }
    }
    rb_ary_store(ary, i, Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew));
  }
  
  free((gsl_vector_complex**) vp);
  return ary;
}

// Starting with version 1.15, GSL provides a gsl_vector_complex_equal
// function, but it only determines absolute equality (i.e. is has no epsilon
// argument).
static int gsl_vector_complex_equal_eps(const gsl_vector_complex *v1,
  const gsl_vector_complex *v2, double eps)
{
  gsl_complex z1, z2;
  size_t i;
  if (v1->size != v2->size) return 0;
  for (i = 0; i < v1->size; i++) {
    z1 = gsl_vector_complex_get(v1, i);
    z2 = gsl_vector_complex_get(v2, i);
    if (!rbgsl_complex_equal(&z1, &z2, eps)) return 0;
  }
  return 1;
  
}

static VALUE rb_gsl_vector_complex_equal(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v1, *v2;
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
  Data_Get_Struct(obj, gsl_vector_complex, v1);
  CHECK_VECTOR_COMPLEX(argv[0]);
  Data_Get_Struct(argv[0], gsl_vector_complex, v2);
  ret = gsl_vector_complex_equal_eps(v1, v2, eps);
  if (ret == 1) return Qtrue;
  else return Qfalse;
}

static VALUE rb_gsl_vector_complex_not_equal(int argc, VALUE *argv, VALUE obj)
{
  VALUE ret;
  ret = rb_gsl_vector_complex_equal(argc, argv, obj);
  if (ret == Qtrue) return Qfalse;
  else return Qtrue;
}

void Init_gsl_vector_complex(VALUE module)
{
  rb_define_singleton_method(cgsl_vector_complex, "new", rb_gsl_vector_complex_new, -1);
  rb_define_singleton_method(cgsl_vector_complex, "[]", rb_gsl_vector_complex_new, -1);
  rb_define_singleton_method(cgsl_vector_complex, "alloc", rb_gsl_vector_complex_new, -1);
  rb_define_singleton_method(cgsl_vector_complex, "calloc", rb_gsl_vector_complex_calloc, 1);

  rb_define_singleton_method(cgsl_vector_complex_col, "new", 
			     rb_gsl_vector_complex_row_new, -1);

  rb_define_method(cgsl_vector_complex, "size", rb_gsl_vector_complex_size, 0);
  rb_define_alias(cgsl_vector_complex, "len", "size");
  rb_define_alias(cgsl_vector_complex, "length", "size");
  rb_define_method(cgsl_vector_complex, "stride", rb_gsl_vector_complex_stride, 0);
  rb_define_method(cgsl_vector_complex, "owner", rb_gsl_vector_complex_owner, 0);
  rb_define_method(cgsl_vector_complex, "get", rb_gsl_vector_complex_get, -1);
  rb_define_alias(cgsl_vector_complex, "[]", "get");
  rb_define_method(cgsl_vector_complex, "ptr", rb_gsl_vector_complex_ptr, 1);

  rb_define_method(cgsl_vector_complex, "set", rb_gsl_vector_complex_set, -1);
  rb_define_alias(cgsl_vector_complex, "[]=", "set");
  rb_define_method(cgsl_vector_complex, "set_all", rb_gsl_vector_complex_set_all, -1);

  rb_define_method(cgsl_vector_complex, "each", rb_gsl_vector_complex_each, 0);
  rb_define_method(cgsl_vector_complex, "reverse_each", rb_gsl_vector_complex_reverse_each, 0);
  rb_define_method(cgsl_vector_complex, "each_index", rb_gsl_vector_complex_each_index, 0);
  rb_define_method(cgsl_vector_complex, "reverse_each_index", rb_gsl_vector_complex_reverse_each_index, 0);
  rb_define_method(cgsl_vector_complex, "collect", rb_gsl_vector_complex_collect, 0);
  rb_define_method(cgsl_vector_complex, "collect!", rb_gsl_vector_complex_collect_bang, 0);
  rb_define_alias(cgsl_vector_complex, "map", "collect");
  rb_define_alias(cgsl_vector_complex, "map!", "collect!");

  rb_define_method(cgsl_vector_complex, "set_zero", rb_gsl_vector_complex_set_zero, 0);
  rb_define_method(cgsl_vector_complex, "set_basis", rb_gsl_vector_complex_set_basis, 1);

  rb_define_method(cgsl_vector_complex, "to_s", rb_gsl_vector_complex_to_s, 0);
  rb_define_method(cgsl_vector_complex, "fprintf", rb_gsl_vector_complex_fprintf, -1);
  rb_define_method(cgsl_vector_complex, "printf", rb_gsl_vector_complex_printf, -1);
  rb_define_method(cgsl_vector_complex, "print", rb_gsl_vector_complex_print, 0);
  rb_define_method(cgsl_vector_complex, "inspect", rb_gsl_vector_complex_inspect, 0);
  rb_define_method(cgsl_vector_complex, "fwrite", rb_gsl_vector_complex_fwrite, 1);
  rb_define_method(cgsl_vector_complex, "fread", rb_gsl_vector_complex_fread, 1);
  rb_define_method(cgsl_vector_complex, "fscanf", rb_gsl_vector_complex_fscanf, 1);

  rb_define_method(cgsl_vector_complex, "real", rb_gsl_vector_complex_real, 0);
  rb_define_alias(cgsl_vector_complex, "re", "real");
  rb_define_method(cgsl_vector_complex, "imag", rb_gsl_vector_complex_imag, 0);
  rb_define_alias(cgsl_vector_complex, "im", "imag");

  rb_define_method(cgsl_vector_complex, "set_real", rb_gsl_vector_complex_set_real, 1);
  rb_define_alias(cgsl_vector_complex, "real=", "set_real");
  rb_define_alias(cgsl_vector_complex, "re=", "set_real");
  rb_define_method(cgsl_vector_complex, "set_imag", rb_gsl_vector_complex_set_imag, 1);
  rb_define_alias(cgsl_vector_complex, "imag=", "set_imag");
  rb_define_alias(cgsl_vector_complex, "im=", "set_imag");

  rb_define_method(cgsl_vector_complex, "conj", rb_gsl_vector_complex_conj, 0);
  rb_define_alias(cgsl_vector_complex, "conjugate", "conj");
  rb_define_method(cgsl_vector_complex, "conj!", rb_gsl_vector_complex_conj_bang, 0);
  rb_define_alias(cgsl_vector_complex, "conjugate!", "conj!");

  rb_define_method(cgsl_vector_complex, "to_a", rb_gsl_vector_complex_to_a, 0);
  rb_define_method(cgsl_vector_complex, "to_a2", rb_gsl_vector_complex_to_a2, 0);

  rb_define_method(cgsl_vector_complex, "subvector", rb_gsl_vector_complex_subvector, -1);
  rb_define_alias(cgsl_vector_complex, "view", "subvector");
  rb_define_method(cgsl_vector_complex, "subvector_with_stride", rb_gsl_vector_complex_subvector_with_stride, 3);
  
  rb_define_singleton_method(cgsl_vector_complex, "memcpy", rb_gsl_vector_complex_memcpy, 2);
  rb_define_method(cgsl_vector_complex, "clone", rb_gsl_vector_complex_clone, 0);
  rb_define_alias(cgsl_vector_complex, "duplicate", "clone");
  rb_define_alias(cgsl_vector_complex, "dup", "clone");
  rb_define_method(cgsl_vector_complex, "reverse!", rb_gsl_vector_complex_reverse, 0);
  rb_define_method(cgsl_vector_complex, "reverse", rb_gsl_vector_complex_reverse2, 0);
  rb_define_method(cgsl_vector_complex, "swap_elements", rb_gsl_vector_complex_swap_elements, 2);
  rb_define_method(cgsl_vector_complex, "fftshift!", rb_gsl_vector_complex_fftshift_bang, 0);
  rb_define_method(cgsl_vector_complex, "fftshift", rb_gsl_vector_complex_fftshift, 0);
  rb_define_method(cgsl_vector_complex, "ifftshift!", rb_gsl_vector_complex_ifftshift_bang, 0);
  rb_define_method(cgsl_vector_complex, "ifftshift", rb_gsl_vector_complex_ifftshift, 0);
  rb_define_method(cgsl_vector_complex, "isnull", rb_gsl_vector_complex_isnull, 0);

  rb_define_method(cgsl_vector_complex, "matrix_view", rb_gsl_vector_complex_matrix_view, -1);
  rb_define_method(cgsl_vector_complex, "matrix_view_with_tda", rb_gsl_vector_complex_matrix_view_with_tda, -1);

  rb_define_method(cgsl_vector_complex, "trans", rb_gsl_vector_complex_trans, 0);
  rb_define_alias(cgsl_vector_complex, "transpose", "trans");
  rb_define_method(cgsl_vector_complex, "trans!", rb_gsl_vector_complex_trans2, 0);
  rb_define_alias(cgsl_vector_complex, "transpose!", "trans!");

  /*****/
  rb_define_alias(cgsl_vector_complex, "col", "trans");
  rb_define_alias(cgsl_vector_complex, "col!", "trans!");

  rb_define_alias(cgsl_vector_complex_col, "row", "trans");
  rb_define_alias(cgsl_vector_complex_col, "row!", "trans!");

   /*****/
  rb_define_method(cgsl_vector_complex, "to_real", rb_gsl_vector_complex_to_real, 0);

  rb_define_method(cgsl_vector_complex, "add", rb_gsl_vector_complex_add, 1);  
  rb_define_method(cgsl_vector_complex, "sub", rb_gsl_vector_complex_sub, 1);  
  rb_define_method(cgsl_vector_complex, "mul", rb_gsl_vector_complex_mul, 1);  
  rb_define_method(cgsl_vector_complex, "div", rb_gsl_vector_complex_div, 1);  
  rb_define_method(cgsl_vector_complex, "add!", rb_gsl_vector_complex_add_bang, 1);  
  rb_define_method(cgsl_vector_complex, "sub!", rb_gsl_vector_complex_sub_bang, 1);  
  rb_define_method(cgsl_vector_complex, "mul!", rb_gsl_vector_complex_mul_bang, 1);  
  rb_define_method(cgsl_vector_complex, "div!", rb_gsl_vector_complex_div_bang, 1);  

  rb_define_alias(cgsl_vector_complex, "+", "add");
  rb_define_alias(cgsl_vector_complex, "add_constant", "add");
  rb_define_alias(cgsl_vector_complex, "add_constant!", "add!");
  rb_define_alias(cgsl_vector_complex, "-", "sub");
  rb_define_alias(cgsl_vector_complex, "*", "mul");
  rb_define_alias(cgsl_vector_complex, "scale", "mul");
  rb_define_alias(cgsl_vector_complex, "scale!", "mul!");
  rb_define_alias(cgsl_vector_complex, "/", "div");

  rb_define_method(cgsl_vector_complex, "coerce", rb_gsl_vector_complex_coerce, 1);  

  /* 2.Aug.2004 */
  rb_define_singleton_method(cgsl_vector_complex, "inner_product", rb_gsl_vector_complex_inner_product, -1);
  rb_define_singleton_method(cgsl_vector_complex, "dot", rb_gsl_vector_complex_inner_product, -1);
  rb_define_method(cgsl_vector_complex, "inner_product", rb_gsl_vector_complex_inner_product, -1);
  /*  rb_define_alias(cgsl_vector_complex, "dot", "inner_product");*/

  /*****/
  rb_define_method(cgsl_vector_complex, "-@", rb_gsl_vector_complex_uminus, 0);
  rb_define_method(cgsl_vector_complex, "+@", rb_gsl_vector_complex_uplus, 0);

  rb_define_method(cgsl_vector_complex, "abs2", rb_gsl_vector_complex_abs2, 0);
  rb_define_alias(cgsl_vector_complex, "square", "abs2");
  rb_define_method(cgsl_vector_complex, "abs", rb_gsl_vector_complex_abs, 0);
  rb_define_alias(cgsl_vector_complex, "amp", "abs");
  rb_define_alias(cgsl_vector_complex, "mag", "abs");
  rb_define_method(cgsl_vector_complex, "arg", rb_gsl_vector_complex_arg, 0);
  rb_define_alias(cgsl_vector_complex, "angle", "arg");
  rb_define_alias(cgsl_vector_complex, "phase", "arg");
  rb_define_method(cgsl_vector_complex, "logabs", rb_gsl_vector_complex_logabs, 0);

  rb_define_method(cgsl_vector_complex, "sqrt", rb_gsl_vector_complex_sqrt, 0);
  rb_define_method(cgsl_vector_complex, "sqrt!", rb_gsl_vector_complex_sqrt_bang, 0);
  rb_define_method(cgsl_vector_complex, "exp", rb_gsl_vector_complex_exp, 0);
  rb_define_method(cgsl_vector_complex, "exp!", rb_gsl_vector_complex_exp_bang, 0);
  rb_define_method(cgsl_vector_complex, "pow", rb_gsl_vector_complex_pow, 1);
  rb_define_method(cgsl_vector_complex, "pow!", rb_gsl_vector_complex_pow_bang, 1);
  rb_define_method(cgsl_vector_complex, "log", rb_gsl_vector_complex_log, 0);
  rb_define_method(cgsl_vector_complex, "log!", rb_gsl_vector_complex_log_bang, 0);
  rb_define_method(cgsl_vector_complex, "log10", rb_gsl_vector_complex_log10, 0);
  rb_define_method(cgsl_vector_complex, "log10!", rb_gsl_vector_complex_log10_bang, 0);
  rb_define_method(cgsl_vector_complex, "log_b", rb_gsl_vector_complex_log_b, 1);
  rb_define_method(cgsl_vector_complex, "log_b!", rb_gsl_vector_complex_log_b_bang, 1);

  rb_define_method(cgsl_vector_complex, "sum", rb_gsl_vector_complex_sum, 0);
  rb_define_method(cgsl_vector_complex, "mean", rb_gsl_vector_complex_mean, 0);
  rb_define_method(cgsl_vector_complex, "tss", rb_gsl_vector_complex_tss, 0);
  rb_define_method(cgsl_vector_complex, "tss_m", rb_gsl_vector_complex_tss_m, 1);
  rb_define_method(cgsl_vector_complex, "variance", rb_gsl_vector_complex_variance, 0);
  rb_define_method(cgsl_vector_complex, "variance_m", rb_gsl_vector_complex_variance_m, 1);
  rb_define_method(cgsl_vector_complex, "variance_fm", rb_gsl_vector_complex_variance_fm, 1);
  rb_define_method(cgsl_vector_complex, "sd", rb_gsl_vector_complex_sd, 0);
  rb_define_method(cgsl_vector_complex, "sd_m", rb_gsl_vector_complex_sd_m, 1);
  rb_define_method(cgsl_vector_complex, "sd_fm", rb_gsl_vector_complex_sd_fm, 1);

  rb_define_method(cgsl_vector_complex, "sin", rb_gsl_vector_complex_sin, 0);
  rb_define_method(cgsl_vector_complex, "cos", rb_gsl_vector_complex_cos, 0);
  rb_define_method(cgsl_vector_complex, "tan", rb_gsl_vector_complex_tan, 0);
  rb_define_method(cgsl_vector_complex, "sec", rb_gsl_vector_complex_sec, 0);
  rb_define_method(cgsl_vector_complex, "csc", rb_gsl_vector_complex_csc, 0);
  rb_define_method(cgsl_vector_complex, "cot", rb_gsl_vector_complex_cot, 0);

  rb_define_method(cgsl_vector_complex, "arcsin", rb_gsl_vector_complex_arcsin, 0);
  rb_define_method(cgsl_vector_complex, "arccos", rb_gsl_vector_complex_arccos, 0);
  rb_define_method(cgsl_vector_complex, "arctan", rb_gsl_vector_complex_arctan, 0);
  rb_define_method(cgsl_vector_complex, "arcsec", rb_gsl_vector_complex_arcsec, 0);
  rb_define_method(cgsl_vector_complex, "arccsc", rb_gsl_vector_complex_arccsc, 0);
  rb_define_method(cgsl_vector_complex, "arccot", rb_gsl_vector_complex_arccot, 0);

  rb_define_method(cgsl_vector_complex, "sinh", rb_gsl_vector_complex_sinh, 0);
  rb_define_method(cgsl_vector_complex, "cosh", rb_gsl_vector_complex_cosh, 0);
  rb_define_method(cgsl_vector_complex, "tanh", rb_gsl_vector_complex_tanh, 0);
  rb_define_method(cgsl_vector_complex, "sech", rb_gsl_vector_complex_sech, 0);
  rb_define_method(cgsl_vector_complex, "csch", rb_gsl_vector_complex_csch, 0);
  rb_define_method(cgsl_vector_complex, "coth", rb_gsl_vector_complex_coth, 0);

  rb_define_method(cgsl_vector_complex, "arcsinh", rb_gsl_vector_complex_arcsinh, 0);
  rb_define_method(cgsl_vector_complex, "arccosh", rb_gsl_vector_complex_arccosh, 0);
  rb_define_method(cgsl_vector_complex, "arctanh", rb_gsl_vector_complex_arctanh, 0);
  rb_define_method(cgsl_vector_complex, "arcsech", rb_gsl_vector_complex_arcsech, 0);
  rb_define_method(cgsl_vector_complex, "arccsch", rb_gsl_vector_complex_arccsch, 0);
  rb_define_method(cgsl_vector_complex, "arccoth", rb_gsl_vector_complex_arccoth, 0);

  /*****/
  rb_define_method(cgsl_vector_complex, "concat", rb_gsl_vector_complex_concat, 1);
  rb_define_method(cgsl_vector_complex, "block", rb_gsl_vector_complex_block, 0);

  rb_define_method(cgsl_vector_complex, "indgen", rb_gsl_vector_complex_indgen, -1);
  rb_define_method(cgsl_vector_complex, "indgen!", rb_gsl_vector_complex_indgen_bang, -1);
  rb_define_singleton_method(cgsl_vector_complex, "indgen", rb_gsl_vector_complex_indgen_singleton, -1);
  rb_define_singleton_method(cgsl_vector_complex, "phasor", rb_gsl_vector_complex_phasor_singleton, -1);

  rb_define_method(cgsl_vector_complex, "zip", rb_gsl_vector_complex_zip, -1);
  rb_define_singleton_method(cgsl_vector_complex, "zip", rb_gsl_vector_complex_zip, -1);
  
  rb_define_method(cgsl_vector_complex, "equal?", rb_gsl_vector_complex_equal, -1);
  rb_define_alias(cgsl_vector_complex, "==", "equal?");
  rb_define_method(cgsl_vector_complex, "not_equal?", rb_gsl_vector_complex_not_equal, -1);
  rb_define_alias(cgsl_vector_complex, "!=", "not_equal?");  
}

