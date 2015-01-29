#include "include/rb_gsl.h"

#ifdef GSL_1_14_LATER

VALUE cMultiset;

VALUE rb_gsl_multiset_alloc(VALUE klass, VALUE nn, VALUE kk)
{

  gsl_multiset *m;
  m = gsl_multiset_alloc(FIX2INT(nn), FIX2INT(kk));
  return Data_Wrap_Struct(klass, 0, gsl_multiset_free, m);
}

VALUE rb_gsl_multiset_calloc(VALUE klass, VALUE nn, VALUE kk)
{
  gsl_multiset *m;
  m = gsl_multiset_alloc(FIX2INT(nn), FIX2INT(kk));
  return Data_Wrap_Struct(klass, 0, gsl_multiset_free, m);
}

VALUE rb_gsl_multiset_init_first(VALUE mm)
{
  gsl_multiset *m;
  Data_Get_Struct(mm, gsl_multiset, m);
  gsl_multiset_init_first(m);
  return Qnil;
}

VALUE rb_gsl_multiset_init_last(VALUE mm)
{
  gsl_multiset *m;
  Data_Get_Struct(mm, gsl_multiset, m);
  gsl_multiset_init_last(m);
  return Qnil;
}

VALUE rb_gsl_multiset_memcpy(VALUE klass, VALUE m1, VALUE m2)
{
  gsl_multiset *mm1, *mm2;
  if (!rb_obj_is_kind_of(m1, klass)) {
    rb_raise(rb_eTypeError, "Wrong type %s (GSL::Multiset expected)", rb_class2name(CLASS_OF(m1)));
  }
  if (!rb_obj_is_kind_of(m2, klass)) {
    rb_raise(rb_eTypeError, "Wrong type %s (GSL::Multiset expected)", rb_class2name(CLASS_OF(m2)));
  }
  Data_Get_Struct(m1, gsl_multiset, mm1);
  Data_Get_Struct(m2, gsl_multiset, mm2);
  return FIX2INT(gsl_multiset_memcpy(mm1, mm2));
}

VALUE rb_gsl_multiset_get(VALUE mm, VALUE i)
{
  gsl_multiset *m;
  Data_Get_Struct(mm, gsl_multiset, m);
  return INT2FIX(gsl_multiset_get(m, FIX2INT(i)));
}

VALUE rb_gsl_multiset_n(VALUE mm)
{
  gsl_multiset *m;
  Data_Get_Struct(mm, gsl_multiset, m);
  return INT2FIX(gsl_multiset_n(m));
}

VALUE rb_gsl_multiset_k(VALUE mm)
{
  gsl_multiset *m;
  Data_Get_Struct(mm, gsl_multiset, m);
  return INT2FIX(gsl_multiset_k(m));
}

VALUE rb_gsl_multiset_valid(VALUE mm)
{
  gsl_multiset *m;
  Data_Get_Struct(mm, gsl_multiset, m);
  return INT2FIX(gsl_multiset_valid(m));
}

VALUE rb_gsl_multiset_valid2(VALUE mm)
{
  gsl_multiset *m;
  Data_Get_Struct(mm, gsl_multiset, m);
  if (gsl_multiset_valid(m)) return Qtrue;
  else return Qfalse;
}

VALUE rb_gsl_multiset_next(VALUE mm)
{
  gsl_multiset *m;
  Data_Get_Struct(mm, gsl_multiset, m);
  return INT2FIX(gsl_multiset_next(m));
}

VALUE rb_gsl_multiset_prev(VALUE mm)
{
  gsl_multiset *m;
  Data_Get_Struct(mm, gsl_multiset, m);
  return INT2FIX(gsl_multiset_prev(m));
}

VALUE rb_gsl_multiset_fwrite(VALUE mm, VALUE name)
{
  gsl_multiset *m;
  FILE *fp = NULL;
  int ret;
  Data_Get_Struct(mm, gsl_multiset, m);
  fp = fopen(STR2CSTR(name), "wb");
  if (fp == NULL) {
    rb_raise(rb_eIOError, "Cannot open file %s", STR2CSTR(name));
  }
  ret = gsl_multiset_fwrite(fp, m);
  fclose(fp);
  return INT2FIX(ret);
}

VALUE rb_gsl_multiset_fread(VALUE mm, VALUE name)
{
  gsl_multiset *m;
  FILE *fp = NULL;
  int ret;
  Data_Get_Struct(mm, gsl_multiset, m);
  fp = fopen(STR2CSTR(name), "wb");
  if (fp == NULL) {
    rb_raise(rb_eIOError, "Cannot open file %s", STR2CSTR(name));
  }
  ret = gsl_multiset_fread(fp, m);
  fclose(fp);
  return INT2FIX(ret);
}

VALUE rb_gsl_multiset_fprintf(VALUE mm, VALUE name, VALUE format)
{
  gsl_multiset *m;
  FILE *fp = NULL;
  int ret;
  Data_Get_Struct(mm, gsl_multiset, m);
  fp = fopen(STR2CSTR(name), "w");
  if (fp == NULL) {
    rb_raise(rb_eIOError, "Cannot open file %s", STR2CSTR(name));
  }
  ret = gsl_multiset_fprintf(fp, m, STR2CSTR(format));
  fclose(fp);
  return INT2FIX(ret);
}

VALUE rb_gsl_multiset_fscanf(VALUE mm, VALUE name)
{
  gsl_multiset *m;
  FILE *fp = NULL;
  int ret;
  Data_Get_Struct(mm, gsl_multiset, m);
  fp = fopen(STR2CSTR(name), "r");
  if (fp == NULL) {
    rb_raise(rb_eIOError, "Cannot open file %s", STR2CSTR(name));
  }
  ret = gsl_multiset_fscanf(fp, m);
  fclose(fp);
  return INT2FIX(ret);
}

VALUE rb_gsl_multiset_data(VALUE mm)
{
  gsl_multiset *m;
  size_t *p;
  gsl_vector_int *v;
  size_t i;
  Data_Get_Struct(mm, gsl_multiset, m);
  p = gsl_multiset_data(m);
  v = gsl_vector_int_alloc(m->k);
  for (i = 0; i < v->size; i++) gsl_vector_int_set(v, i, p[i]);
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, v);
}

VALUE rb_gsl_multiset_data2(VALUE mm, VALUE i)
{
  gsl_multiset *m;
  size_t *p;
  Data_Get_Struct(mm, gsl_multiset, m);
  p = gsl_multiset_data(m);
  return INT2FIX(p[i]);
}

void Init_multiset(VALUE module)
{
  cMultiset = rb_define_class_under(module, "Multiset", cGSL_Object);
  rb_define_singleton_method(cMultiset, "alloc", rb_gsl_multiset_alloc, 2);
  rb_define_singleton_method(cMultiset, "calloc", rb_gsl_multiset_calloc, 2);
  rb_define_singleton_method(cMultiset, "memcpy", rb_gsl_multiset_memcpy, 2);

  rb_define_method(cMultiset, "init_first", rb_gsl_multiset_init_first, 0);
  rb_define_method(cMultiset, "init_last", rb_gsl_multiset_init_last, 0);

  rb_define_method(cMultiset, "get", rb_gsl_multiset_get, 1);
  rb_define_alias(cMultiset, "[]", "get");

  rb_define_method(cMultiset, "n", rb_gsl_multiset_n, 0);
  rb_define_method(cMultiset, "k", rb_gsl_multiset_k, 0);
  rb_define_method(cMultiset, "data", rb_gsl_multiset_data, 0);
  rb_define_method(cMultiset, "data[]", rb_gsl_multiset_data2, 1);
  rb_define_method(cMultiset, "valid", rb_gsl_multiset_valid, 0);
  rb_define_method(cMultiset, "valid?", rb_gsl_multiset_valid2, 0);

  rb_define_method(cMultiset, "next", rb_gsl_multiset_next, 0);
  rb_define_method(cMultiset, "prev", rb_gsl_multiset_prev, 0);

  rb_define_method(cMultiset, "fwrite", rb_gsl_multiset_fwrite, 1);
  rb_define_method(cMultiset, "fread", rb_gsl_multiset_fread, 1);
  rb_define_method(cMultiset, "fprintf", rb_gsl_multiset_fprintf, 2);
  rb_define_method(cMultiset, "fscanf", rb_gsl_multiset_fscanf, 1);
}

#endif

