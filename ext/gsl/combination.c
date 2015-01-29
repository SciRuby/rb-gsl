/*
  combination.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#ifdef GSL_1_1_LATER
#include "include/rb_gsl_common.h"
#include "include/rb_gsl_array.h"

static VALUE cgsl_combination_data;

static VALUE rb_gsl_combination_new(VALUE klass, VALUE n, VALUE k)
{
  gsl_combination *c = NULL;
  CHECK_FIXNUM(n);CHECK_FIXNUM(k);
  c = gsl_combination_alloc(FIX2INT(n), FIX2INT(k));
  return Data_Wrap_Struct(klass, 0, gsl_combination_free, c);
}

static VALUE rb_gsl_combination_calloc(VALUE klass, VALUE n, VALUE k)
{
  gsl_combination *c = NULL;
  CHECK_FIXNUM(n);CHECK_FIXNUM(k);
  c = gsl_combination_calloc(FIX2INT(n), FIX2INT(k));
  return Data_Wrap_Struct(klass, 0, gsl_combination_free, c);
}

static VALUE rb_gsl_combination_init_first(VALUE obj)
{
  gsl_combination *c = NULL;
  Data_Get_Struct(obj, gsl_combination, c);
  gsl_combination_init_first(c);
  return obj;
}

static VALUE rb_gsl_combination_init_last(VALUE obj)
{
  gsl_combination *c = NULL;
  Data_Get_Struct(obj, gsl_combination, c);
  gsl_combination_init_last(c);
  return obj;
}

#ifdef GSL_1_4_LATER
/* singleton */
static VALUE rb_gsl_combination_memcpy(VALUE klass, VALUE dst, VALUE src)
{
  gsl_combination *c, *c2;
  if (!rb_obj_is_kind_of(dst, klass))
    rb_raise(rb_eTypeError, "wrong argument type %s (Combination expected)",
       rb_class2name(CLASS_OF(dst)));
  if (!rb_obj_is_kind_of(src, klass))
    rb_raise(rb_eTypeError, "wrong argument type %s (Combination expected)",
       rb_class2name(CLASS_OF(src)));

  Data_Get_Struct(dst, gsl_combination, c2);
  Data_Get_Struct(src, gsl_combination, c);
  gsl_combination_memcpy(c2, c);
  return dst;
}

static VALUE rb_gsl_combination_clone(VALUE obj)
{
  gsl_combination *c, *c2;
  Data_Get_Struct(obj, gsl_combination, c);
  c2 = gsl_combination_alloc(c->n, c->k);
  gsl_combination_memcpy(c2, c);
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_combination_free, c2);
}
#endif

static VALUE rb_gsl_combination_get(VALUE obj, VALUE ii)
{
  gsl_combination *c = NULL;
  size_t i;
  CHECK_FIXNUM(ii);
  Data_Get_Struct(obj, gsl_combination, c);
  i = FIX2INT(ii);
  if (i > c->n) rb_raise(rb_eIndexError, "index out of range");
  return INT2FIX(gsl_combination_get(c, i));
}

static VALUE rb_gsl_combination_set(VALUE obj, VALUE ii, VALUE val)
{
  gsl_combination *c = NULL;
  size_t i;
  CHECK_FIXNUM(ii);
  CHECK_FIXNUM(val);
  Data_Get_Struct(obj, gsl_combination, c);
  i = FIX2INT(ii);
  c->data[i] = FIX2INT(val);
  return obj;
}

static VALUE rb_gsl_combination_n(VALUE obj)
{
  gsl_combination *c = NULL;
  Data_Get_Struct(obj, gsl_combination, c);
  return INT2FIX(gsl_combination_n(c));

}
static VALUE rb_gsl_combination_k(VALUE obj)
{
  gsl_combination *c = NULL;
  Data_Get_Struct(obj, gsl_combination, c);
  return INT2FIX(gsl_combination_k(c));
}

static VALUE rb_gsl_combination_data(VALUE obj)
{
  gsl_combination *c = NULL;
  gsl_permutation *p = NULL;
  Data_Get_Struct(obj, gsl_combination, c);
  p = ALLOC(gsl_permutation);
  p->size = c->k;
  p->data = c->data;
  return Data_Wrap_Struct(cgsl_combination_data, 0, free, p);
}

static VALUE rb_gsl_combination_valid(VALUE obj)
{
  gsl_combination *c = NULL;
  Data_Get_Struct(obj, gsl_combination, c);
  return INT2FIX(gsl_combination_valid(c));
}

static VALUE rb_gsl_combination_valid2(VALUE obj)
{
  gsl_combination *c = NULL;
  Data_Get_Struct(obj, gsl_combination, c);
  if(gsl_combination_valid(c)) return Qtrue;
  else return Qfalse;
}

static VALUE rb_gsl_combination_next(VALUE obj)
{
  gsl_combination *c = NULL;
  Data_Get_Struct(obj, gsl_combination, c);
  return INT2FIX(gsl_combination_next(c));
}
static VALUE rb_gsl_combination_prev(VALUE obj)
{
  gsl_combination *c = NULL;
  Data_Get_Struct(obj, gsl_combination, c);
  return INT2FIX(gsl_combination_prev(c));
}

static VALUE rb_gsl_combination_fwrite(VALUE obj, VALUE io)
{
  gsl_combination *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_combination, h);
  f = rb_gsl_open_writefile(io, &flag);
  status = gsl_combination_fwrite(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_combination_fread(VALUE obj, VALUE io)
{
  gsl_combination *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;

  Data_Get_Struct(obj, gsl_combination, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = gsl_combination_fread(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_combination_fprintf(int argc, VALUE *argv, VALUE obj)
{
  gsl_combination *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;

  if (argc != 1 && argc != 2) rb_raise(rb_eArgError,
               "wrong number of arguments (%d for 1 or 2)", argc);

  Data_Get_Struct(obj, gsl_combination, h);
  fp = rb_gsl_open_writefile(argv[0], &flag);
  switch (argc) {
  case 1:
    status = gsl_combination_fprintf(fp, h, "%u\n");
    break;
  default:
    Check_Type(argv[1], T_STRING);
    status = gsl_combination_fprintf(fp, h, STR2CSTR(argv[1]));
    break;
  }
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE rb_gsl_combination_printf(int argc, VALUE *argv, VALUE obj)
{
  gsl_combination *h = NULL;
  int status;
  Data_Get_Struct(obj, gsl_combination, h);
  switch (argc) {
  case 0:
    status = gsl_combination_fprintf(stdout, h, "%u\n");
    break;
  default:
    Check_Type(argv[0], T_STRING);
    status = gsl_combination_fprintf(stdout, h, STR2CSTR(argv[0]));
    break;
  }
  return INT2FIX(status);
}


static VALUE rb_gsl_combination_fscanf(VALUE obj, VALUE io)
{
  gsl_combination *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_combination, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = gsl_combination_fscanf(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_combination_equal(VALUE obj, VALUE other)
{
  gsl_combination *p1 = NULL, *p2 = NULL;
  size_t i;
  Data_Get_Struct(obj, gsl_combination, p1);
  Data_Get_Struct(other, gsl_combination, p2);
  if (p1->k != p2->k) return Qfalse;
  for (i = 0; i < p1->k; i++)
    if (p1->data[i] != p2->data[i]) return Qfalse;
  return Qtrue;
}

void Init_gsl_combination(VALUE module)
{
  VALUE cgsl_combination;
  cgsl_combination = rb_define_class_under(module, "Combination", cGSL_Object);
  cgsl_combination_data = rb_define_class_under(cgsl_combination, "Data",
            cgsl_permutation);
  rb_define_singleton_method(cgsl_combination, "new", rb_gsl_combination_new, 2);
  rb_define_singleton_method(cgsl_combination, "alloc", rb_gsl_combination_new, 2);
  rb_define_singleton_method(cgsl_combination, "calloc", rb_gsl_combination_calloc, 2);
  rb_define_method(cgsl_combination, "init_first", rb_gsl_combination_init_first, 0);
  rb_define_method(cgsl_combination, "init_last", rb_gsl_combination_init_last, 0);
#ifdef GSL_1_4_LATER
  rb_define_singleton_method(cgsl_combination, "memcpy", rb_gsl_combination_memcpy, 2);
  rb_define_method(cgsl_combination, "clone", rb_gsl_combination_clone, 0);
#endif
  rb_define_method(cgsl_combination, "get", rb_gsl_combination_get, 1);
  rb_define_alias(cgsl_combination, "[]", "get");
  rb_define_method(cgsl_combination, "set", rb_gsl_combination_set, 2);
  rb_define_alias(cgsl_combination, "[]=", "set");

  rb_define_method(cgsl_combination, "n", rb_gsl_combination_n, 0);
  rb_define_method(cgsl_combination, "k", rb_gsl_combination_k, 0);
  rb_define_method(cgsl_combination, "data", rb_gsl_combination_data, 0);
  rb_define_method(cgsl_combination, "valid", rb_gsl_combination_valid, 0);
  rb_define_method(cgsl_combination, "valid?", rb_gsl_combination_valid2, 0);
  rb_define_method(cgsl_combination, "next", rb_gsl_combination_next, 0);
  rb_define_method(cgsl_combination, "prev", rb_gsl_combination_prev, 0);

  rb_define_method(cgsl_combination, "fwrite", rb_gsl_combination_fwrite, 1);
  rb_define_method(cgsl_combination, "fread", rb_gsl_combination_fread, 1);
  rb_define_method(cgsl_combination, "fprintf", rb_gsl_combination_fprintf, -1);
  rb_define_method(cgsl_combination, "printf", rb_gsl_combination_printf, -1);
  rb_define_method(cgsl_combination, "fscanf", rb_gsl_combination_fscanf, 1);

  rb_define_method(cgsl_combination, "equal?", rb_gsl_combination_equal, 1);
  rb_define_alias(cgsl_combination, "==", "equal?");
}
#endif
