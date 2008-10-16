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
void get_range_beg_en_n(VALUE range, int *beg, int *en, size_t *n, int *step);


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
      n = RARRAY(argv[0])->len;
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

static VALUE rb_gsl_vector_complex_get(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL, *vnew;
  gsl_complex *c = NULL, z;
  int i;
  int beg, en, step;
  size_t index, n, j;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  switch (argc) {
  case 1:
    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      CHECK_FIXNUM(argv[0]);
      i = FIX2INT(argv[0]);
      if (i < 0) index = (size_t) (v->size + i);
      else index = (size_t) i;
      c = ALLOC(gsl_complex);
      *c = gsl_vector_complex_get(v, index);
      return Data_Wrap_Struct(cgsl_complex, 0, free, c);
      break;
    case T_ARRAY:
      vnew = gsl_vector_complex_alloc(RARRAY(argv[0])->len);
      for (j = 0; j < vnew->size; j++) {
	i = FIX2INT(rb_ary_entry(argv[0], j));
	if (i < 0) i = v->size + i;
	gsl_vector_complex_set(vnew, j, gsl_vector_complex_get(v, i));
      }
      return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
      break;
    default:
      if (CLASS_OF(argv[0]) == rb_cRange) {
	get_range_beg_en_n(argv[0], &beg, &en, &n, &step);
	vnew = gsl_vector_complex_alloc(n);
	for (j = 0; j < n; j++) {
	  z = gsl_vector_complex_get(v, j+beg);
	  gsl_vector_complex_set(vnew, j, z);
	}
	return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
      }
      break;
    }
    break;
  case 0:
    rb_raise(rb_eArgError, "number of arguments must be > 0");
    break;
  default:
    vnew = gsl_vector_complex_alloc(argc);
    for (j = 0; j < argc; j++) {
      i = FIX2INT(argv[j]);
      if (i >= 0) z = gsl_vector_complex_get(v, i);
      else z = gsl_vector_complex_get(v, v->size+i);
      gsl_vector_complex_set(vnew, j, z);
    }
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
    break;
  }
  return Qnil;
}

static VALUE rb_gsl_vector_complex_set(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_complex tmp, *z = &tmp;
  size_t i;
  if (argc < 2) rb_raise(rb_eArgError, "wrong number of arguments");
  Data_Get_Struct(obj, gsl_vector_complex, v);
  if (TYPE(argv[0]) != T_FIXNUM) {
    for (i = 0; i < argc; i++) {
      if (i >= v->size) break;
      switch (TYPE(argv[i])) {
      case T_ARRAY:
	tmp.dat[0] = NUM2DBL(rb_ary_entry(argv[i], 0));
	tmp.dat[1] = NUM2DBL(rb_ary_entry(argv[i], 1));
	break;
      case T_FLOAT:
      case T_FIXNUM:
      case T_BIGNUM:
	tmp = gsl_complex_rect(NUM2DBL(argv[i]), 0.0);
	break;
      default:
	CHECK_COMPLEX(argv[i]);
	Data_Get_Struct(argv[i], gsl_complex, z);
	tmp = *z;
	break;
      }
      gsl_vector_complex_set(v, i, tmp);
    }
    return obj;
  }
  i = FIX2INT(argv[0]);
  switch (argc) {
  case 2:
    if (rb_obj_is_kind_of(argv[1], cgsl_complex)) {
      Data_Get_Struct(argv[1], gsl_complex, z);
      tmp = *z;
    } else {
      switch(TYPE(argv[1])) {
      case T_ARRAY:
	tmp.dat[0] = NUM2DBL(rb_ary_entry(argv[1], 0));
	tmp.dat[1] = NUM2DBL(rb_ary_entry(argv[1], 1));
	break;
      case T_FLOAT:
      case T_FIXNUM:
      case T_BIGNUM:
	tmp = gsl_complex_rect(NUM2DBL(argv[1]), 0.0);
	break;
      default:
	CHECK_COMPLEX(argv[1]);
	Data_Get_Struct(argv[1], gsl_complex, z);
	tmp = *z;
	break;
      }
    }
    break;
  case 3:
    tmp = gsl_complex_rect(NUM2DBL(argv[1]), NUM2DBL(argv[2]));
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  gsl_vector_complex_set(v, i, tmp);

  return obj;
}

static VALUE rb_gsl_vector_complex_set_all(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_complex *z = NULL, tmp;
  if (argc < 1) rb_raise(rb_eArgError, "wrong number of arguments");
  Data_Get_Struct(obj, gsl_vector_complex, v);
  switch (argc) {
  case 1:
    if (rb_obj_is_kind_of(argv[0], cgsl_complex)) {
      Data_Get_Struct(argv[0], gsl_complex, z);
      tmp = *z;
    } else {
      switch(TYPE(argv[0])) {
      case T_ARRAY:
	tmp.dat[0] = NUM2DBL(rb_ary_entry(argv[0], 0));
	tmp.dat[1] = NUM2DBL(rb_ary_entry(argv[0], 1));
	break;
      case T_FLOAT:
      case T_FIXNUM:
      case T_BIGNUM:
	tmp.dat[0] = NUM2DBL(argv[0]);
	tmp.dat[1] = 0.0;
	break;
      default:
	rb_raise(rb_eTypeError, 
		 "wrong argument type %s", rb_class2name(CLASS_OF(argv[0])));
	break;
      }
    }
    break;
  case 2:
    tmp = gsl_complex_rect(NUM2DBL(argv[0]), NUM2DBL(argv[1]));
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }

  gsl_vector_complex_set_all(v, tmp);

  return obj;
}

static VALUE rb_gsl_vector_complex_each(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_complex *zp = NULL, ztmp, *zp2 = &ztmp;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  for (i = 0; i < v->size; i++) {
    zp = GSL_COMPLEX_AT(v, i);
    ztmp = *zp;  /* create a temporal copy of *zp */
    rb_yield(Data_Wrap_Struct(cgsl_complex, 0, NULL, zp2));
  }
  return obj;
}

static VALUE rb_gsl_vector_complex_reverse_each(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_complex *zp = NULL, ztmp, *zp2 = &ztmp;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  for (i = v->size-1; i >= 0; i--) {
    zp = GSL_COMPLEX_AT(v, i);
    ztmp = *zp;  /* create a temporal copy of *zp */
    rb_yield(Data_Wrap_Struct(cgsl_complex, 0, NULL, zp2));
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
  for (i = v->size-1; i >= 0; i--) {
    rb_yield(INT2FIX(i));
    if (i == 0) break;
  }
  return obj;
}

static VALUE rb_gsl_vector_complex_collect(VALUE obj)
{
  gsl_vector_complex *v = NULL, *vnew;
  size_t size;
  gsl_complex *zp = NULL, ztmp, *zp2 = &ztmp, *zp3;
  size_t i;
  VALUE a;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vnew = gsl_vector_complex_alloc(v->size);
  size = v->size;
  for (i = 0; i < size; i++) {
    zp = GSL_COMPLEX_AT(v, i);
    ztmp = *zp;
    a = rb_yield(Data_Wrap_Struct(cgsl_complex, 0, NULL, zp2));
    CHECK_COMPLEX(a);
    Data_Get_Struct(a, gsl_complex, zp3);
    gsl_vector_complex_set(vnew, i, *zp3);
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
}

static VALUE rb_gsl_vector_complex_collect_bang(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  size_t size;
  gsl_complex *zp = NULL, ztmp, *zp2 = &ztmp, *zp3;
  size_t i;
  VALUE a;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  size = v->size;
  for (i = 0; i < size; i++) {
    zp = GSL_COMPLEX_AT(v, i);
    ztmp = *zp;
    a = rb_yield(Data_Wrap_Struct(cgsl_complex, 0, NULL, zp2));
    CHECK_COMPLEX(a);
    Data_Get_Struct(a, gsl_complex, zp3);
    gsl_vector_complex_set(v, i, *zp3);
  }
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
  size_t offset, n;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  switch (argc) {
  case 0:
    offset = 0;
    n = v->size;
    break;
  case 1:
    CHECK_FIXNUM(argv[0]);
    offset = 0;
    n = FIX2INT(argv[0]);
    break;
  case 2:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    offset = FIX2INT(argv[0]);
    n = FIX2INT(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0, 1 or 2)", argc);
    break;
  }
  vv = gsl_vector_complex_view_alloc();
  *vv = gsl_vector_complex_subvector(v, offset, n);
  if (VECTOR_COMPLEX_ROW_P(obj))
    return Data_Wrap_Struct(cgsl_vector_complex_view, 0, gsl_vector_complex_view_free, vv);
  else
    return Data_Wrap_Struct(cgsl_vector_complex_col_view, 0, gsl_vector_complex_view_free, vv);
}

static VALUE rb_gsl_vector_complex_subvector_with_stride(VALUE obj, VALUE o, VALUE s, VALUE nn)
{
  gsl_vector_complex *v = NULL;
  gsl_vector_complex_view *vv = NULL;
  CHECK_FIXNUM(o); CHECK_FIXNUM(nn); CHECK_FIXNUM(s);
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vv = gsl_vector_complex_view_alloc();
  *vv = gsl_vector_complex_subvector_with_stride(v, FIX2INT(o), FIX2INT(s), FIX2INT(nn));
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
    RBASIC(obj)->klass = cgsl_vector_complex_col;
  else if (CLASS_OF(obj) == cgsl_vector_complex_col) 
    RBASIC(obj)->klass = cgsl_vector_complex;
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
};

static VALUE rb_gsl_vector_complex_arithmetics(int flag, VALUE obj, VALUE bb);

static VALUE rb_gsl_vector_complex_arithmetics(int flag, VALUE obj, VALUE bb) 
{
  gsl_vector *b = NULL;
  gsl_vector_complex *cv = NULL, *cvnew = NULL, *cb = NULL;
  gsl_complex *c = NULL, z;
  Data_Get_Struct(obj, gsl_vector_complex, cv);
  switch (TYPE(bb)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    z = gsl_complex_rect(NUM2DBL(bb), 0.0);
    switch (flag) {
    case GSL_VECTOR_COMPLEX_ADD:
      cvnew = make_vector_complex_clone(cv);
      gsl_vector_complex_add_constant(cvnew, z);
      return Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj), 0, gsl_vector_complex_free, cvnew);
      break;
    case GSL_VECTOR_COMPLEX_SUB:
      cvnew = make_vector_complex_clone(cv);
      gsl_vector_complex_add_constant(cvnew, gsl_complex_negative(z));
      return Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj), 0, gsl_vector_complex_free, cvnew);
      break;
    case GSL_VECTOR_COMPLEX_MUL:
      cvnew = make_vector_complex_clone(cv);
      gsl_vector_complex_scale(cvnew, z);
      return Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj), 0, gsl_vector_complex_free, cvnew);
      break;
    case GSL_VECTOR_COMPLEX_DIV:
      cvnew = make_vector_complex_clone(cv);
      gsl_vector_complex_scale(cvnew, gsl_complex_inverse(z));
      return Data_Wrap_Struct(VECTOR_COMPLEX_ROW_COL(obj), 0, gsl_vector_complex_free, cvnew);
      break;
    default:
      rb_raise(rb_eRuntimeError, "unknown operation");
      break;
    }
    break;
  default:
    if (VECTOR_P(bb)) {
      Data_Get_Struct(bb, gsl_vector, b);
      cb = vector_to_complex(b);
      switch (flag) {
      case GSL_VECTOR_COMPLEX_ADD:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_add(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_SUB:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_sub(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_MUL:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_mul(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_DIV:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_div(cvnew, cb);
	break;
      default:
	rb_raise(rb_eRuntimeError, "unknown operation");
	break;
      }
      gsl_vector_complex_free(cb);
      if (VECTOR_COMPLEX_ROW_P(obj)) 
	return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cvnew);
      else
	return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, cvnew);
    } else if (VECTOR_COMPLEX_P(bb)) {
      Data_Get_Struct(bb, gsl_vector_complex, cb);
      switch (flag) {
      case GSL_VECTOR_COMPLEX_ADD:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_add(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_SUB:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_sub(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_MUL:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_mul(cvnew, cb);
	break;
      case GSL_VECTOR_COMPLEX_DIV:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_div(cvnew, cb);
	break;
      default:
	rb_raise(rb_eRuntimeError, "unknown operation");
	break;
      }
      if (VECTOR_COMPLEX_ROW_P(obj)) 
	return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cvnew);
      else
	return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, cvnew);
    } else if (COMPLEX_P(bb)) {
      Data_Get_Struct(bb, gsl_complex, c);
      switch (flag) {
      case GSL_VECTOR_COMPLEX_ADD:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_add_constant(cvnew, *c);
	break;
      case GSL_VECTOR_COMPLEX_SUB:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_add_constant(cvnew, gsl_complex_negative(*c));
	break;
      case GSL_VECTOR_COMPLEX_MUL:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_scale(cvnew, *c);
	break;
      case GSL_VECTOR_COMPLEX_DIV:
	cvnew = make_vector_complex_clone(cv);
	gsl_vector_complex_scale(cvnew, gsl_complex_inverse(*c));
	break;
      default:
	rb_raise(rb_eRuntimeError, "unknown operation");
	break;
      }
      if (VECTOR_COMPLEX_ROW_P(obj))
	return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cvnew);
      else
	return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, cvnew);
    } else {
      rb_raise(rb_eTypeError, "wrong type argument %s", rb_class2name(CLASS_OF(bb)));
    }
    break;
  }
  return Qnil;
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

static VALUE rb_gsl_vector_complex_scale_bang(VALUE obj, VALUE s)
{
  gsl_vector_complex *m;
  gsl_complex c, *z = &c;
  Data_Get_Struct(obj, gsl_vector_complex, m);
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
  gsl_vector_complex_scale(m, *z);
  return obj;
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

static VALUE rb_gsl_vector_complex_exp(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_exp);
}

static VALUE rb_gsl_vector_complex_pow(VALUE obj, VALUE a)
{
  return rb_gsl_vector_complex_XXXz2(obj, a, gsl_complex_pow);
}

static VALUE rb_gsl_vector_complex_log(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_log);
}

static VALUE rb_gsl_vector_complex_log10(VALUE obj)
{
  return rb_gsl_vector_complex_XXXz(obj, gsl_complex_log10);
}

static VALUE rb_gsl_vector_complex_log_b(VALUE obj, VALUE a)
{
  return rb_gsl_vector_complex_XXXz2(obj, a, gsl_complex_log_b);
}

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

static VALUE rb_gsl_vector_complex_shift(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_complex *z = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  if (v->stride != 1) rb_raise(rb_eRuntimeError, "vector must have stride 1");
  if (v->size == 0) return Qnil;
  z = ALLOC(gsl_complex);
  *z = gsl_vector_complex_get(v, 0);
  v->size -= 1;
  memmove(v->block->data, v->block->data+2, sizeof(double)*v->size*2);
  return Data_Wrap_Struct(cgsl_complex, 0, free, z);
}

static VALUE rb_gsl_vector_complex_pop(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_complex *z = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);
  if (v->stride != 1) rb_raise(rb_eRuntimeError, "vector must have stride 1");
  if (v->size == 0) return Qnil;
  z = ALLOC(gsl_complex);
  *z = gsl_vector_complex_get(v, v->size-1);
  v->size -= 1;
  return Data_Wrap_Struct(cgsl_complex, 0, free, z);
}

static VALUE rb_gsl_vector_complex_unshift(VALUE obj, VALUE x)
{
  gsl_vector_complex *v = NULL;
  gsl_block_complex *b = NULL, *bnew = NULL;
  gsl_complex *z, c;
  switch (TYPE(x)) {
  case T_FIXNUM: case T_FLOAT:
    c.dat[0] = NUM2DBL(x); c.dat[1] = 0.0;
    break;
  case T_ARRAY:
    c.dat[0] = NUM2DBL(rb_ary_entry(x, 0));
    c.dat[1] = NUM2DBL(rb_ary_entry(x, 1));
    break;
  default:
    CHECK_COMPLEX(x);
    Data_Get_Struct(x, gsl_complex, z);
    c = *z;
    break;
  }
  Data_Get_Struct(obj, gsl_vector_complex, v);
  if (v->stride != 1) rb_raise(rb_eRuntimeError, "vector must have stride 1");
  b = v->block;
  if (b->size < (v->size+1)*2) {
    bnew = gsl_block_complex_alloc(v->size + 1);
    memcpy(bnew->data+2, b->data, sizeof(double)*b->size*2);
    gsl_block_complex_free(b);
  } else {
    bnew = b;
    memmove(bnew->data+2, v->data, sizeof(double)*v->size*2);
  }
  v->data = bnew->data;
  v->block = bnew;
  v->size += 1;
  gsl_vector_complex_set(v, 0, c);
  return obj;
}

static VALUE rb_gsl_vector_complex_push(VALUE obj, VALUE x)
{
  gsl_vector_complex *v = NULL;
  gsl_block_complex *b = NULL, *bnew = NULL;
  gsl_complex *z, c;
  switch (TYPE(x)) {
  case T_FIXNUM: case T_FLOAT:
    c.dat[0] = NUM2DBL(x); c.dat[1] = 0.0;
    break;
  case T_ARRAY:
    c.dat[0] = NUM2DBL(rb_ary_entry(x, 0));
    c.dat[1] = NUM2DBL(rb_ary_entry(x, 1));
    break;
  default:
    CHECK_COMPLEX(x);
    Data_Get_Struct(x, gsl_complex, z);
    c = *z;
    break;
  }
  Data_Get_Struct(obj, gsl_vector_complex, v);
  if (v->stride != 1) rb_raise(rb_eRuntimeError, "vector must have stride 1");
  b = v->block;
  if (b->size < (v->size+1)*2) {
    bnew = gsl_block_complex_alloc(v->size + 1);
    memcpy(bnew->data, b->data, sizeof(double)*b->size*2);
    gsl_block_complex_free(b);
  } else {
    bnew = b;
    memmove(bnew->data+2, b->data, sizeof(double)*v->size*2);
  }
  v->block = bnew;
  v->size += 1;
  v->data = bnew->data;
  gsl_vector_complex_set(v, v->size-1, c);
  return obj;
}

static VALUE rb_gsl_vector_complex_concat(VALUE obj, VALUE other)
{
  gsl_vector_complex *v = NULL, *v2 = NULL;
  gsl_block_complex *bnew = NULL;
  CHECK_VECTOR_COMPLEX(other);
  Data_Get_Struct(obj, gsl_vector_complex, v);
  if (v->stride != 1) rb_raise(rb_eRuntimeError, "vector must have stride 1");
  Data_Get_Struct(other, gsl_vector_complex, v2);
  bnew = gsl_block_complex_alloc(v->size + v2->size);
  memcpy(bnew->data, v->block->data, sizeof(double)*v->block->size*2);
  memcpy(bnew->data + v->block->size*2, v2->block->data, 
	 sizeof(double)*v2->block->size*2);
  gsl_block_complex_free(v->block);
  v->size += v2->size;
  v->block = bnew;
  v->data = bnew->data;
  return obj;
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
  size_t start = 0, step = 1;
  size_t i, j;
  switch (argc) {
  case 0:
    break;
  case 1:
    start = FIX2INT(argv[0]);
    break;
  case 2:
    start = FIX2INT(argv[0]);
    step = FIX2INT(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-2)");
  }
  Data_Get_Struct(obj, gsl_vector_complex, v);
  for (i = 0, j = start; i < v->size; i++, j += step) {
    v->data[2*i] = (double) j;
  }
  return obj;
}

static VALUE rb_gsl_vector_complex_indgen(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL, *vnew;
  size_t start = 0, step = 1;
  size_t i, j;
  switch (argc) {
  case 0:
    break;
  case 1:
    start = FIX2INT(argv[0]);
    break;
  case 2:
    start = FIX2INT(argv[0]);
    step = FIX2INT(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-2)");
  }
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vnew = gsl_vector_complex_calloc(v->size);
  for (i = 0, j = start; i < vnew->size; i++, j += step) {
    vnew->data[2*i] = (double) j;
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
}

static VALUE rb_gsl_vector_complex_indgen_singleton(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *vnew;
  size_t n, start = 0, step = 1;
  size_t i, j;
  switch (argc) {
  case 1:
    n = (size_t) FIX2INT(argv[0]);
    break;
  case 2:
    n = (size_t) FIX2INT(argv[0]);
    start = FIX2INT(argv[1]);
    break;
  case 3:
    n = (size_t) FIX2INT(argv[0]);
    start = FIX2INT(argv[1]);
    step = FIX2INT(argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-3)");
  }
  vnew = gsl_vector_complex_calloc(n);
  for (i = 0, j = start; i < vnew->size; i++, j += step) {
    vnew->data[2*i] = (double) j;
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
  for (i = 0; i < argc2; i++) {
    CHECK_VECTOR_COMPLEX(argv2[i]);
  }
  vp = (gsl_vector_complex**) malloc(sizeof(gsl_vector_complex**));
  for (i = 0; i < argc2; i++) {
    Data_Get_Struct(argv2[i], gsl_vector_complex, vp[i]);
  }
  ary = rb_ary_new2(v0->size);
  for (i = 0; i < v0->size; i++) {
    vnew = gsl_vector_complex_alloc(argc2 + 1);
    gsl_vector_complex_set(vnew, 0, gsl_vector_complex_get(v0, i));
    for (j = 0; j < argc2; j++) {
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

static int gsl_vector_complex_equal(const gsl_vector_complex *v1,
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
  ret = gsl_vector_complex_equal(v1, v2, eps);
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

  rb_define_method(cgsl_vector_complex, "set_zero", rb_gsl_vector_complex_set_zero, 0);
  rb_define_method(cgsl_vector_complex, "set_basis", rb_gsl_vector_complex_set_basis, 1);

  rb_define_method(cgsl_vector_complex, "fprintf", rb_gsl_vector_complex_fprintf, -1);
  rb_define_method(cgsl_vector_complex, "printf", rb_gsl_vector_complex_printf, -1);
  rb_define_method(cgsl_vector_complex, "print", rb_gsl_vector_complex_print, 0);
  rb_define_method(cgsl_vector_complex, "inspect", rb_gsl_vector_complex_print, 0);
  rb_define_method(cgsl_vector_complex, "fwrite", rb_gsl_vector_complex_fwrite, 1);
  rb_define_method(cgsl_vector_complex, "fread", rb_gsl_vector_complex_fread, 1);
  rb_define_method(cgsl_vector_complex, "fscanf", rb_gsl_vector_complex_fscanf, 1);

  rb_define_method(cgsl_vector_complex, "real", rb_gsl_vector_complex_real, 0);
  rb_define_alias(cgsl_vector_complex, "re", "real");
  rb_define_method(cgsl_vector_complex, "imag", rb_gsl_vector_complex_imag, 0);
  rb_define_alias(cgsl_vector_complex, "im", "imag");

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

  rb_define_alias(cgsl_vector_complex, "+", "add");
  rb_define_alias(cgsl_vector_complex, "add_constant", "add");
  rb_define_alias(cgsl_vector_complex, "-", "sub");
  rb_define_alias(cgsl_vector_complex, "*", "mul");
  rb_define_alias(cgsl_vector_complex, "scale", "mul");
  rb_define_alias(cgsl_vector_complex, "/", "div");

  rb_define_method(cgsl_vector_complex, "scale!", rb_gsl_vector_complex_scale_bang, 1);
  rb_define_method(cgsl_vector_complex, "coerce", rb_gsl_vector_complex_coerce, 1);  

  /* 2.Aug.2004 */
  rb_define_singleton_method(cgsl_vector, "inner_product", rb_gsl_vector_complex_inner_product, -1);
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
  rb_define_alias(cgsl_vector_complex, "phase", "arg");
  rb_define_method(cgsl_vector_complex, "logabs", rb_gsl_vector_complex_logabs, 0);

  rb_define_method(cgsl_vector_complex, "sqrt", rb_gsl_vector_complex_sqrt, 0);
  rb_define_method(cgsl_vector_complex, "exp", rb_gsl_vector_complex_exp, 0);
  rb_define_method(cgsl_vector_complex, "pow", rb_gsl_vector_complex_pow, 1);
  rb_define_method(cgsl_vector_complex, "log", rb_gsl_vector_complex_log, 0);
  rb_define_method(cgsl_vector_complex, "log10", rb_gsl_vector_complex_log10, 0);
  rb_define_method(cgsl_vector_complex, "log_b", rb_gsl_vector_complex_log_b, 1);

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
  rb_define_method(cgsl_vector_complex, "shift", rb_gsl_vector_complex_shift, 0);
  rb_define_method(cgsl_vector_complex, "pop", rb_gsl_vector_complex_pop, 0);
  rb_define_method(cgsl_vector_complex, "unshift", rb_gsl_vector_complex_unshift, 1);
  rb_define_method(cgsl_vector_complex, "push", rb_gsl_vector_complex_push, 1);
  rb_define_method(cgsl_vector_complex, "concat", rb_gsl_vector_complex_concat, 1);
  rb_define_method(cgsl_vector_complex, "block", rb_gsl_vector_complex_block, 0);

  rb_define_method(cgsl_vector_complex, "indgen", rb_gsl_vector_complex_indgen, -1);
  rb_define_method(cgsl_vector_complex, "indgen!", rb_gsl_vector_complex_indgen_bang, -1);
  rb_define_singleton_method(cgsl_vector_complex, "indgen", rb_gsl_vector_complex_indgen_singleton, -1);

  rb_define_method(cgsl_vector_complex, "zip", rb_gsl_vector_complex_zip, -1);
  rb_define_singleton_method(cgsl_vector_complex, "zip", rb_gsl_vector_complex_zip, -1);
  
  rb_define_method(cgsl_vector_complex, "equal?", rb_gsl_vector_complex_equal, -1);
  rb_define_alias(cgsl_vector_complex, "==", "equal?");
  rb_define_method(cgsl_vector_complex, "not_equal?", rb_gsl_vector_complex_not_equal, -1);
  rb_define_alias(cgsl_vector_complex, "!=", "not_equal?");  
}

