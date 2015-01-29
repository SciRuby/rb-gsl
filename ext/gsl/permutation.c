/*
  permutation.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_array.h"
#include <gsl/gsl_permute.h>
#include <gsl/gsl_permute_vector.h>

VALUE rb_gsl_permutation_alloc(VALUE klass, VALUE nn);

/*
 * Creates a new permutation of size n. The permutation is not initialized
 * and its elements are undefined. Use the method calloc if you want to create
 * a permutation which is initialized to the identity.
 */
VALUE rb_gsl_permutation_alloc(VALUE klass, VALUE nn)
{
  gsl_permutation *p = NULL;
  CHECK_FIXNUM(nn);
  p = gsl_permutation_calloc(FIX2INT(nn));
  return Data_Wrap_Struct(klass, 0, gsl_permutation_free, p);
}

static VALUE rb_gsl_permutation_calloc(VALUE klass, VALUE nn)
{
  gsl_permutation *p = NULL;
  CHECK_FIXNUM(nn);
  p = gsl_permutation_calloc(FIX2INT(nn));
  return Data_Wrap_Struct(klass, 0, gsl_permutation_free, p);
}

static VALUE rb_gsl_permutation_size(VALUE obj)
{
  gsl_permutation *p = NULL;
  Data_Get_Struct(obj, gsl_permutation, p);
  return INT2FIX(p->size);
}

static VALUE rb_gsl_permutation_init(VALUE obj)
{
  gsl_permutation *p = NULL;
  Data_Get_Struct(obj, gsl_permutation, p);
  gsl_permutation_init(p);
  return obj;
}

void get_range_int_beg_en_n(VALUE range, int *beg, int *en, size_t *n, int *step);
static VALUE rb_gsl_permutation_get(int argc, VALUE *argv, VALUE obj)
{
  gsl_permutation *b, *bnew;
  gsl_index *p;
  int beg, en, i, step;
  size_t n, j, k;
  Data_Get_Struct(obj, gsl_permutation, b);
  switch (argc) {
  case 0:
    rb_raise(rb_eArgError, "too few arguments (%d for >= 1)", argc);
    break;
  case 1:
    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      i = FIX2INT(argv[0]);
      if (i < 0) {
        if (i >= -((int) b->size)) j = b->size + i;
        else rb_raise(rb_eRangeError, "offset %d out of range", i);
      } else {
        if (i < ((int) b->size)) j = (size_t) i;
        else rb_raise(rb_eRangeError, "offset %d out of range", i);
      }
      return INT2FIX((int) b->data[j]);
      break;
    case T_ARRAY:
      //      n = RARRAY(argv[0])->len;
      n = RARRAY_LEN(argv[0]);
      bnew = gsl_permutation_alloc(n);
      for (j = 0; j < n; j++) {
        i = FIX2INT(rb_ary_entry(argv[0], j));
        if (i < 0) {
          if (i >= -((int) b->size)) k = b->size + i;
          else rb_raise(rb_eRangeError, "offset %d out of range", i);
        } else {
          if (i < ((int) b->size)) k = i;
          else rb_raise(rb_eRangeError, "offset %d out of range", i);
        }
        bnew->data[j] = b->data[k];
      }
      return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_permutation_free, bnew);
      break;
    default:
      if (PERMUTATION_P(argv[0])) {
  Data_Get_Struct(argv[0], gsl_index, p);
  bnew = gsl_permutation_alloc(p->size);
  for (j = 0; j < p->size; j++) bnew->data[j] = b->data[p->data[j]];
  return Data_Wrap_Struct(CLASS_OF(argv[0]), 0, gsl_permutation_free, bnew);
      } else if (CLASS_OF(argv[0]) == rb_cRange) {
        get_range_int_beg_en_n(argv[0], &beg, &en, &n, &step);
        if (beg < -((int) b->size) || en >= ((int) b->size))
          rb_raise(rb_eRangeError, "range overflow (%d..%d)", beg, en);
        bnew = gsl_permutation_alloc(n);
        for (j = 0; j < n; j++)
          if (beg+(int)j < 0)
            bnew->data[j] = b->data[b->size+beg+j];
          else
            bnew->data[j] = b->data[beg+j];
        return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_permutation_free, bnew);
      } else {
  rb_raise(rb_eArgError, "wrong argument type %s (Fixnum, Array, or Range expected)", rb_class2name(CLASS_OF(argv[0])));
  break;
      }
    }
    break;
  default:
    bnew = gsl_permutation_alloc(argc);
    for (j = 0; j < (size_t) argc; j++) {
      i = FIX2INT(argv[j]);
      if (i < 0) k = b->size + i; else k = i;
      bnew->data[j] = b->data[k];
    }
    return Data_Wrap_Struct(CLASS_OF(argv[0]), 0, gsl_permutation_free, bnew);
    break;
  }
  return Qnil;
}

#ifdef GSL_1_1_LATER
static VALUE rb_gsl_permutation_clone(VALUE obj)
{
  gsl_permutation *p, *p2 = NULL;
  Data_Get_Struct(obj, gsl_permutation, p);
  p2 = gsl_permutation_alloc(p->size);
  gsl_permutation_memcpy(p2, p);
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_permutation_free, p2);
}

/* singleton */
static VALUE rb_gsl_permutation_memcpy(VALUE obj, VALUE pp1, VALUE pp2)
{
  gsl_permutation *p1, *p2 = NULL;
  CHECK_PERMUTATION(pp1);
  CHECK_PERMUTATION(pp2);
  Data_Get_Struct(pp1, gsl_permutation, p1);
  Data_Get_Struct(pp2, gsl_permutation, p2);
  gsl_permutation_memcpy(p1, p2);
  return pp1;
}
#endif

static VALUE rb_gsl_permutation_swap(VALUE obj, VALUE i, VALUE j)
{
  gsl_permutation *p = NULL;
  CHECK_FIXNUM(i); CHECK_FIXNUM(j);
  Data_Get_Struct(obj, gsl_permutation, p);
  gsl_permutation_swap(p, FIX2INT(i), FIX2INT(j));
  return obj;
}

static VALUE rb_gsl_permutation_valid(VALUE obj)
{
  gsl_permutation *p = NULL;
  Data_Get_Struct(obj, gsl_permutation, p);
  return INT2FIX(gsl_permutation_valid(p));
}

static VALUE rb_gsl_permutation_valid2(VALUE obj)
{
  gsl_permutation *p = NULL;
  Data_Get_Struct(obj, gsl_permutation, p);
  if(gsl_permutation_valid(p)) return Qtrue;
  else return Qfalse;
}

/* to array */
static VALUE rb_gsl_permutation_to_a(VALUE obj)
{
  gsl_permutation *p = NULL;
  size_t i;
  VALUE ary;
  Data_Get_Struct(obj, gsl_permutation, p);
  ary = rb_ary_new2(p->size);
  for (i = 0; i < p->size; i++) {
    rb_ary_store(ary, i, INT2FIX(gsl_permutation_get(p, i)));
  }
  return ary;
}

/* to vector */
static VALUE rb_gsl_permutation_to_v(VALUE obj)
{
  gsl_permutation *p = NULL;
  gsl_vector *v;
  size_t size;
  size_t i;
  Data_Get_Struct(obj, gsl_permutation, p);
  size = p->size;
  v = gsl_vector_alloc(size);
  for (i = 0; i < size; i++) {
    gsl_vector_set(v, i, gsl_permutation_get(p, i));
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

static VALUE rb_gsl_permutation_reverse(VALUE obj)
{
  gsl_permutation *p = NULL;
  Data_Get_Struct(obj, gsl_permutation, p);
  gsl_permutation_reverse(p);
  return obj;
}

static VALUE rb_gsl_permutation_inverse(VALUE obj)
{
  gsl_permutation *p, *inv;
  Data_Get_Struct(obj, gsl_permutation, p);
  inv = gsl_permutation_alloc(p->size);
  gsl_permutation_inverse(inv, p);
  return Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, inv);
}

static VALUE rb_gsl_permutation_next(VALUE obj)
{
  gsl_permutation *p = NULL;
  Data_Get_Struct(obj, gsl_permutation, p);
  return INT2FIX(gsl_permutation_next(p));
}

static VALUE rb_gsl_permutation_prev(VALUE obj)
{
  gsl_permutation *p = NULL;
  Data_Get_Struct(obj, gsl_permutation, p);
  return INT2FIX(gsl_permutation_prev(p));
}

static VALUE rb_gsl_permutation_permute_vector(VALUE obj, VALUE vv)
{
  gsl_permutation *p = NULL;
  gsl_vector *v;
  int status;
  CHECK_VECTOR(vv);
  Data_Get_Struct(obj, gsl_permutation, p);
  Data_Get_Struct(vv, gsl_vector, v);
  status = gsl_permute_vector(p, v);
  return INT2FIX(status);
}

static VALUE rb_gsl_permutation_permute_vector_inverse(VALUE obj, VALUE vv)
{
  gsl_permutation *p = NULL;
  gsl_vector *v;
  int status;
  CHECK_VECTOR(vv);
  Data_Get_Struct(obj, gsl_permutation, p);
  Data_Get_Struct(vv, gsl_vector, v);
  status = gsl_permute_vector_inverse(p, v);
  return INT2FIX(status);
}

/* singleton */
static VALUE rb_gsl_permute_vector(VALUE obj, VALUE pp, VALUE vv)
{
  gsl_permutation *p = NULL;
  gsl_vector *v;
  int status;
  CHECK_VECTOR(vv);
  Data_Get_Struct(pp, gsl_permutation, p);
  Data_Get_Struct(vv, gsl_vector, v);
  status = gsl_permute_vector(p, v);
  return INT2FIX(status);
}

/* singleton */
static VALUE rb_gsl_permute_vector_inverse(VALUE obj, VALUE pp, VALUE vv)
{
  gsl_permutation *p = NULL;
  gsl_vector *v;
  int status;
  CHECK_VECTOR(vv);
  Data_Get_Struct(pp, gsl_permutation, p);
  Data_Get_Struct(vv, gsl_vector, v);
  status = gsl_permute_vector_inverse(p, v);
  return INT2FIX(status);
}

/* singleton method */
#ifdef GSL_1_2_LATER
static VALUE rb_gsl_permutation_mul(VALUE obj, VALUE ppa, VALUE ppb)
{
  gsl_permutation *p = NULL;
  gsl_permutation *pa = NULL, *pb = NULL;
  int flag = 0;
  CHECK_PERMUTATION(ppa);
  CHECK_PERMUTATION(ppb);
  Data_Get_Struct(ppa, gsl_permutation, pa);
  Data_Get_Struct(ppb, gsl_permutation, pb);
  if (rb_obj_is_kind_of(obj, cgsl_permutation)) {
    Data_Get_Struct(obj, gsl_permutation, p);
    flag = 1;
  } else {
    p = gsl_permutation_alloc(pa->size);
  }
  gsl_permutation_mul(p, pa, pb);
  if (flag == 1) return obj;
  else return Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
}
#endif

static VALUE rb_gsl_permutation_print(VALUE obj);

static VALUE rb_gsl_permutation_print(VALUE obj)
{
  gsl_permutation *p = NULL;
  size_t size, i;
  Data_Get_Struct(obj, gsl_permutation, p);
  size = p->size;
  for (i = 0; i < size; i++) {
    printf("%3d ", (int) gsl_permutation_get(p, i));
    if ((i+1)%10 == 0) printf("\n");
  }
  printf("\n");
  return obj;
}

static VALUE rb_gsl_permutation_to_s(VALUE obj)
{
  gsl_permutation *v = NULL;
  char buf[16];
  size_t i;
  VALUE str;
  Data_Get_Struct(obj, gsl_permutation, v);
  str = rb_str_new2("[");
  for (i = 0; i < v->size; i++) {
    sprintf(buf,  " %d", (int) gsl_permutation_get(v, i));
    rb_str_cat(str, buf, strlen(buf));
  }
  sprintf(buf, " ]");
  rb_str_cat(str, buf, strlen(buf));
  return str;
}

static VALUE rb_gsl_permutation_inspect(VALUE obj)
{
  VALUE str;
  char buf[64];
  sprintf(buf, "%s\n", rb_class2name(CLASS_OF(obj)));
  str = rb_str_new2(buf);
  return rb_str_concat(str, rb_gsl_permutation_to_s(obj));
}

static VALUE rb_gsl_permutation_fwrite(VALUE obj, VALUE io)
{
  gsl_permutation *h;
  FILE *f = NULL;
  int status, flag = 0;

  Data_Get_Struct(obj, gsl_permutation, h);
  f = rb_gsl_open_writefile(io, &flag);
  status = gsl_permutation_fwrite(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_permutation_fread(VALUE obj, VALUE io)
{
  gsl_permutation *h;
  FILE *f = NULL;
  int status, flag = 0;

  Data_Get_Struct(obj, gsl_permutation, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = gsl_permutation_fread(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_permutation_fprintf(int argc, VALUE *argv, VALUE obj)
{
  gsl_permutation *h;
  FILE *fp = NULL;
  int status, flag = 0;

  if (argc != 1 && argc != 2) rb_raise(rb_eArgError,
               "wrong number of arguments (%d for 1 or 2)", argc);

  Data_Get_Struct(obj, gsl_permutation, h);
  fp = rb_gsl_open_writefile(argv[0], &flag);
  if (argc == 1) {
    status = gsl_permutation_fprintf(fp, h, "%u\n");
  } else {
    Check_Type(argv[1], T_STRING);
    status = gsl_permutation_fprintf(fp, h, STR2CSTR(argv[1]));
  }
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE rb_gsl_permutation_printf(int argc, VALUE *argv, VALUE obj)
{
  gsl_permutation *h;
  int status;

  Data_Get_Struct(obj, gsl_permutation, h);
  if (argc == 0) {
    status = gsl_permutation_fprintf(stdout, h, "%u\n");
  } else {
    Check_Type(argv[0], T_STRING);
    status = gsl_permutation_fprintf(stdout, h, STR2CSTR(argv[0]));
  }
  return INT2FIX(status);
}

static VALUE rb_gsl_permutation_fscanf(VALUE obj, VALUE io)
{
  gsl_permutation *h;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_permutation, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = gsl_permutation_fscanf(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

#ifdef GSL_1_2_LATER
static VALUE rb_gsl_permutation_linear_to_canonical(int argc, VALUE *argv, VALUE obj)
{
  gsl_permutation *p, *q;

  Data_Get_Struct(obj, gsl_permutation, p);
  switch (argc) {
  case 0:
    q = gsl_permutation_alloc(p->size);
    gsl_permutation_linear_to_canonical(q, p);
    return Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, q);
    break;
  case 1:
    CHECK_PERMUTATION(argv[0]);
    Data_Get_Struct(argv[0], gsl_permutation, q);
    gsl_permutation_linear_to_canonical(q, p);
    return obj;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
  }
  return Qtrue;
}

static VALUE rb_gsl_permutation_canonical_to_linear(int argc, VALUE *argv, VALUE obj)
{
  gsl_permutation *p, *q;

  Data_Get_Struct(obj, gsl_permutation, p);
  switch (argc) {
  case 0:
    q = gsl_permutation_alloc(p->size);
    gsl_permutation_canonical_to_linear(q, p);
    return Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, q);
    break;
  case 1:
    CHECK_PERMUTATION(argv[0]);
    Data_Get_Struct(argv[0], gsl_permutation, q);
    gsl_permutation_canonical_to_linear(q, p);
    return obj;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
  }
  return Qtrue;
}

static VALUE rb_gsl_permutation_inversions(VALUE obj)
{
  gsl_permutation *p = NULL;
  Data_Get_Struct(obj, gsl_permutation, p);
  return INT2FIX(gsl_permutation_inversions(p));
}

static VALUE rb_gsl_permutation_linear_cycles(VALUE obj)
{
  gsl_permutation *p = NULL;
  Data_Get_Struct(obj, gsl_permutation, p);
  return INT2FIX(gsl_permutation_linear_cycles(p));
}

static VALUE rb_gsl_permutation_canonical_cycles(VALUE obj)
{
  gsl_permutation *p = NULL;
  Data_Get_Struct(obj, gsl_permutation, p);
  return INT2FIX(gsl_permutation_canonical_cycles(p));
}
#endif

static VALUE rb_gsl_vector_permute(VALUE obj, VALUE pp)
{
  gsl_permutation *p = NULL;
  gsl_vector *v = NULL;
  int status;
  CHECK_PERMUTATION(pp);
  Data_Get_Struct(pp, gsl_permutation, p);
  Data_Get_Struct(obj, gsl_vector, v);
  status = gsl_permute_vector(p, v);
  return INT2FIX(status);
}

static VALUE rb_gsl_vector_permute_inverse(VALUE obj, VALUE pp)
{
  gsl_permutation *p = NULL;
  gsl_vector *v = NULL;
  int status;
  CHECK_PERMUTATION(pp);
  Data_Get_Struct(pp, gsl_permutation, p);
  Data_Get_Struct(obj, gsl_vector, v);
  status = gsl_permute_vector_inverse(p, v);
  return INT2FIX(status);
}

static VALUE rb_gsl_permutation_set(VALUE obj, VALUE ii, VALUE val)
{
  gsl_permutation *p = NULL;
  CHECK_FIXNUM(ii);
  CHECK_FIXNUM(val);
  Data_Get_Struct(obj, gsl_permutation, p);
  p->data[FIX2INT(ii)] = FIX2INT(val);
  return obj;
}

static VALUE rb_gsl_permutation_equal(VALUE obj, VALUE other)
{
  gsl_permutation *p1 = NULL, *p2 = NULL;
  size_t i;
  CHECK_PERMUTATION(other);
  Data_Get_Struct(obj, gsl_permutation, p1);
  Data_Get_Struct(other, gsl_permutation, p2);
  if (p1->size != p2->size) return Qfalse;
  for (i = 0; i < p1->size; i++)
    if (p1->data[i] != p2->data[i]) return Qfalse;
  return Qtrue;
}

void Init_gsl_permutation(VALUE module)
{
  rb_define_singleton_method(cgsl_permutation, "alloc", rb_gsl_permutation_alloc, 1);
  rb_define_singleton_method(cgsl_permutation, "calloc", rb_gsl_permutation_calloc, 1);
  rb_define_method(cgsl_permutation, "size", rb_gsl_permutation_size, 0);
  rb_define_method(cgsl_permutation, "init", rb_gsl_permutation_init, 0);
  rb_define_method(cgsl_permutation, "inspect", rb_gsl_permutation_inspect, 0);
  rb_define_method(cgsl_permutation, "to_s", rb_gsl_permutation_to_s, 0);
  rb_define_method(cgsl_permutation, "get", rb_gsl_permutation_get, -1);
  rb_define_alias(cgsl_permutation, "[]", "get");
  rb_define_method(cgsl_permutation, "set", rb_gsl_permutation_set, 2);
  rb_define_alias(cgsl_permutation, "[]=", "set");
#ifdef GSL_1_1_LATER
  rb_define_singleton_method(cgsl_permutation, "memcpy", rb_gsl_permutation_memcpy, 2);
  rb_define_method(cgsl_permutation, "clone", rb_gsl_permutation_clone, 0);
#endif
  rb_define_method(cgsl_permutation, "swap", rb_gsl_permutation_swap, 2);
  rb_define_method(cgsl_permutation, "valid", rb_gsl_permutation_valid, 0);
  rb_define_method(cgsl_permutation, "valid?", rb_gsl_permutation_valid2, 0);
  rb_define_method(cgsl_permutation, "to_a", rb_gsl_permutation_to_a, 0);
  rb_define_method(cgsl_permutation, "to_v", rb_gsl_permutation_to_v, 0);
  rb_define_method(cgsl_permutation, "reverse", rb_gsl_permutation_reverse, 0);
  rb_define_method(cgsl_permutation, "inverse", rb_gsl_permutation_inverse, 0);
  rb_define_alias(cgsl_permutation, "inv", "inverse");
  rb_define_method(cgsl_permutation, "next", rb_gsl_permutation_next, 0);
  rb_define_method(cgsl_permutation, "prev", rb_gsl_permutation_prev, 0);

  rb_define_method(cgsl_permutation, "permute_vector", rb_gsl_permutation_permute_vector, 1);
  rb_define_alias(cgsl_permutation, "permute", "permute_vector");
  rb_define_method(cgsl_permutation, "permute_vector_inverse", rb_gsl_permutation_permute_vector_inverse, 1);
  rb_define_alias(cgsl_permutation, "permute_inverse", "permute_vector_inverse");

  rb_define_singleton_method(cgsl_permutation, "permute_vector", rb_gsl_permute_vector, 2);
  rb_define_singleton_method(cgsl_permutation, "permute_vector_inverse", rb_gsl_permute_vector_inverse, 2);
  rb_define_module_function(module, "permute_vector", rb_gsl_permute_vector, 2);
  rb_define_module_function(module, "permute_vector_inverse", rb_gsl_permute_vector_inverse, 2);

  rb_define_singleton_method(cgsl_permutation, "permute", rb_gsl_permute_vector, 2);
  rb_define_singleton_method(cgsl_permutation, "permute_inverse", rb_gsl_permute_vector_inverse, 2);
  rb_define_module_function(module, "permute", rb_gsl_permute_vector, 2);
  rb_define_module_function(module, "permute_inverse", rb_gsl_permute_vector_inverse, 2);

  rb_define_method(cgsl_permutation, "fwrite", rb_gsl_permutation_fwrite, 1);
  rb_define_method(cgsl_permutation, "fread", rb_gsl_permutation_fread, 1);
  rb_define_method(cgsl_permutation, "fprintf", rb_gsl_permutation_fprintf, -1);
  rb_define_method(cgsl_permutation, "printf", rb_gsl_permutation_printf, -1);
  rb_define_method(cgsl_permutation, "fscanf", rb_gsl_permutation_fscanf, 1);
  rb_define_method(cgsl_permutation, "print", rb_gsl_permutation_print, 0);

#ifdef GSL_1_2_LATER
  rb_define_singleton_method(cgsl_permutation, "mul", rb_gsl_permutation_mul, 2);
  rb_define_method(cgsl_permutation, "mul", rb_gsl_permutation_mul, 2);
  rb_define_method(cgsl_permutation, "linear_to_canonical", rb_gsl_permutation_linear_to_canonical, -1);
  rb_define_alias(cgsl_permutation, "to_canonical", "linear_to_canonical");
  rb_define_method(cgsl_permutation, "canonical_to_linear", rb_gsl_permutation_canonical_to_linear, -1);
  rb_define_alias(cgsl_permutation, "to_linear", "canonical_to_linear");

  rb_define_method(cgsl_permutation, "inversions", rb_gsl_permutation_inversions, 0);
  rb_define_method(cgsl_permutation, "linear_cycles", rb_gsl_permutation_linear_cycles, 0);
  rb_define_method(cgsl_permutation, "canonical_cycles", rb_gsl_permutation_canonical_cycles, 0);
#endif

  rb_define_method(cgsl_vector, "permute", rb_gsl_vector_permute, 1);
  rb_define_method(cgsl_vector, "permute_inverse", rb_gsl_vector_permute_inverse, 1);

  rb_define_method(cgsl_permutation, "equal?", rb_gsl_permutation_equal, 1);
  rb_define_alias(cgsl_permutation, "==", "equal?");

}
