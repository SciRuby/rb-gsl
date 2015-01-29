/*
  common.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_array.h"
#include "include/rb_gsl_histogram.h"
#include <string.h>
#include <ctype.h>

FILE* rb_gsl_open_writefile(VALUE io, int *flag)
{
#ifdef HAVE_RUBY_IO_H
  rb_io_t *fptr = NULL;
#else
  OpenFile *fptr = NULL;
#endif
  FILE *fp = NULL;
  char *name;
  switch (TYPE(io)) {
  case T_STRING:
    name = RSTRING_PTR(io);
    fp = fopen(name, "w");
    *flag = 1;
    break;
  case T_FILE:
    GetOpenFile(io, fptr);
    /*
#ifdef HAVE_RUBY_IO_H
    name = STR2CSTR(fptr->pathv);
#else
    name = fptr->path;
#endif
    */
    rb_io_check_writable(fptr);
#ifdef HAVE_RUBY_IO_H
    fp = rb_io_stdio_file(fptr);
#else
    fp = GetWriteFile(fptr);
#endif
    *flag = 0;
    break;
  default:
    rb_raise(rb_eTypeError, "argv 1 String or File expected");
    break;
  }
  //  if (fp == NULL) rb_raise(rb_eIOError, "Cannot open file %s.", name);
  if (fp == NULL) rb_raise(rb_eIOError, "Cannot open file.");
  return fp;
}

FILE* rb_gsl_open_readfile(VALUE io, int *flag)
{
#ifdef HAVE_RUBY_IO_H
  rb_io_t *fptr = NULL;
#else
  OpenFile *fptr = NULL;
#endif
  FILE *fp = NULL;
  char *name;
  switch (TYPE(io)) {
  case T_STRING:
    name = RSTRING_PTR(io);
    fp = fopen(name, "r");
    *flag = 1;
    break;
  case T_FILE:
    GetOpenFile(io, fptr);
    /*
#ifdef HAVE_RUBY_IO_H
    name = STR2CSTR(fptr->pathv);
#else
    name = fptr->path;
#endif
    */
    rb_io_check_readable(fptr);
#ifdef HAVE_RUBY_IO_H
    fp = rb_io_stdio_file(fptr);
#else
    fp = fptr->f;
#endif
    *flag = 0;
    break;
  default:
    rb_raise(rb_eTypeError, "argv 1 String or File expected");
    break;
  }
  //  if (fp == NULL) rb_raise(rb_eIOError, "Cannot open file %s.", name);
  if (fp == NULL) rb_raise(rb_eIOError, "Cannot open file");
  return fp;
}

VALUE rb_gsl_obj_read_only(int argc, VALUE *argv, VALUE obj)
{
  rb_raise(rb_eRuntimeError, "Read only object.");
}

int str_tail_grep(const char *s0, const char *s1)
{
  int len0, len1;
  char *p = NULL;
  len0 = strlen(s0);
  len1 = strlen(s1);
  p = (char *) s0 + len0 - len1;
  return strcmp(p, s1);
}

int str_head_grep(const char *s0, const char *s1)
{
  int len0, len1;
  size_t i, len;
  char *p0, *p1;
  len0 = strlen(s0);
  len1 = strlen(s1);
  len = (size_t) GSL_MIN_INT(len0, len1);
  p0 = (char *) s0;
  p1 = (char *) s1;
  for (i = 0; i < len; i++) if (*p0++ != *p1++) return 1;
  return 0;
}

size_t count_columns(const char *str)
{
  size_t n = 0;
  int flag = 1;
  char *p;
  p = (char *) str;
  do {
    if (isspace(*p)) {
      flag = 1;
    } else {
      if (flag == 1) {
        flag = 0;
        n++;
      }
    }
    p++;
  } while (*p != '\0' && *p != '\n');
  return n;
}

char* str_scan_double(const char *str, double *val)
{
  char buf[256];
  char *p, *q;
  double x;
  int flag = 0;
  p = (char *) str;
  q = buf;
  do {
    if (isspace(*p)) {
      if (flag == 0) {
        /* do nothing */
      } else {
        break;
      }
    } else {
      *q++ = *p;
      flag = 1;
    }
    p++;
  } while (*p != '\0' && *p != '\n');
  if (flag == 0) {
    *val = 0;
    return NULL;
  }
  *q = '\0';
  flag = sscanf(buf, "%lf", &x);
  if (flag == 1) {
    *val = x;
    return p;
  } else {
    *val = 0;
    return NULL;
  }
}

char* str_scan_int(const char *str, int *val)
{
  char buf[256];
  char *p, *q;
  int x;
  int flag = 0;
  p = (char *) str;
  q = buf;
  do {
    if (isspace(*p)) {
      if (flag == 0) {
        /* do nothing */
      } else {
        break;
      }
    } else {
      *q++ = *p;
      flag = 1;
    }
    p++;
  } while (*p != '\0' && *p != '\n');
  if (flag == 0) {
    *val = 0;
    return NULL;
  }
  *q = '\0';
  flag = sscanf(buf, "%d", &x);
  if (flag == 1) {
    *val = x;
    return p;
  } else {
    *val = 0;
    return NULL;
  }
}

double* get_ptr_double3(VALUE obj, size_t *size, size_t *stride, int *flag)
{
  gsl_vector *v;
#ifdef HAVE_NARRAY_H
  double *ptr;
  struct NARRAY *na;
  if (NA_IsNArray(obj)) {
    obj = na_change_type(obj, NA_DFLOAT);
    GetNArray(obj, na);
    ptr = (double *) na->ptr;
    *size = na->total;
    *stride = 1;
    *flag = 1;
    return ptr;
  }
#endif
  CHECK_VECTOR(obj);
  Data_Get_Struct(obj, gsl_vector, v);
  *size = v->size;
  *stride = v->stride;
  *flag = 0;
  return v->data;
}

gsl_complex ary2complex(VALUE obj)
{
  gsl_complex *z, c;
  switch (TYPE(obj)) {
  case T_ARRAY:
    GSL_SET_REAL(&c, NUM2DBL(rb_ary_entry(obj, 0)));
    GSL_SET_IMAG(&c, NUM2DBL(rb_ary_entry(obj, 1)));
    break;
  default:
    if (COMPLEX_P(obj)) {
      Data_Get_Struct(obj, gsl_complex, z);
      c = *z;
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (Array or Complex expected)",
         rb_class2name(CLASS_OF(obj)));
    }
    break;
  }
  return c;
}

VALUE vector_eval_create(VALUE obj, double (*func)(double))
{
  gsl_vector *vnew;
  size_t i, size, stride;
  double *ptr;
  ptr = get_vector_ptr(obj, &stride, &size);
  vnew = gsl_vector_alloc(size);
  for (i = 0; i < size; i++) {
    gsl_vector_set(vnew, i, (*func)(ptr[i*stride]));
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
}

VALUE matrix_eval_create(VALUE obj, double (*func)(double))
{
  gsl_matrix *m, *mnew;
  size_t i, j;
  Data_Get_Struct(obj, gsl_matrix, m);
  mnew = gsl_matrix_alloc(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      gsl_matrix_set(mnew, i, j, (*func)(gsl_matrix_get(m, i, j)));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
}

VALUE rb_gsl_ary_eval1(VALUE ary, double (*f)(double))
{
  VALUE ary2;
  size_t i, n;
  double val;
  //  n = RARRAY(ary)->len;
  n = RARRAY_LEN(ary);
  ary2 = rb_ary_new2(n);
  for (i = 0; i < n; i++) {
    val = (*f)(NUM2DBL(rb_ary_entry(ary, i)));
    rb_ary_store(ary2, i, rb_float_new(val));
  }
  return ary2;
}

#ifdef HAVE_NARRAY_H
VALUE rb_gsl_nary_eval1(VALUE ary, double (*f)(double))
{
  VALUE ary2;
  struct NARRAY *na;
  double *ptr1, *ptr2;
  size_t i, n;
  ary = na_change_type(ary, NA_DFLOAT);
  GetNArray(ary, na);
  ptr1 = (double *) na->ptr;
  n = na->total;
  ary2 = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(ary));
  ptr2 = NA_PTR_TYPE(ary2, double*);
  for (i = 0; i < n; i++) ptr2[i] = (*f)(ptr1[i]);
  return ary2;
}
#endif
