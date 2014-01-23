/*
  block_source.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2005 by Yoshiki Tsunesada
                               Cameron McBride

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#ifdef BASE_DOUBLE
#define CHECK_BL CHECK_BLOCK
#define C_TO_VALUE rb_float_new
#define BL_P BLOCK_P
#elif defined(BASE_INT)
#define C_TO_VALUE INT2FIX
#define CHECK_BL CHECK_BLOCK_INT
#define BL_P BLOCK_INT_P
#else
#define C_TO_VALUE INT2FIX
#define CHECK_BL CHECK_BLOCK_UCHAR
#define BL_P BLOCK_UCHAR_P
#endif

static VALUE FUNCTION(rb_gsl_block,new)(VALUE klass, VALUE nn)
{
  GSL_TYPE(gsl_block) *block = NULL;
  CHECK_FIXNUM(nn);
  block = FUNCTION(gsl_block,alloc)(FIX2INT(nn));
  return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_block,free), block);
}

static VALUE FUNCTION(rb_gsl_block,calloc)(VALUE klass, VALUE nn)
{
  GSL_TYPE(gsl_block) *block = NULL;
  CHECK_FIXNUM(nn);
  block = FUNCTION(gsl_block,calloc)(FIX2INT(nn));
  return Data_Wrap_Struct(klass, 0, FUNCTION(gsl_block,free), block);
}

static VALUE FUNCTION(rb_gsl_block,size)(VALUE obj)
{
  GSL_TYPE(gsl_block) *block = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), block);
  return INT2FIX(block->size);
}

static VALUE FUNCTION(rb_gsl_block,fwrite)(VALUE obj, VALUE io)
{
  GSL_TYPE(gsl_block) *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), h);
  f = rb_gsl_open_writefile(io, &flag);
  status = FUNCTION(gsl_block,fwrite)(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_gsl_block,fread)(VALUE obj, VALUE io)
{
  GSL_TYPE(gsl_block) *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), h);
  f = rb_gsl_open_readfile(io, &flag);
  status = FUNCTION(gsl_block,fread)(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

#ifdef BASE_DOUBLE
#define FORMAT_DEFAULT "%g"
#else
#define FORMAT_DEFAULT "%d"
#endif

static VALUE FUNCTION(rb_gsl_block,fprintf)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_block) *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  if (argc != 1 && argc != 2) 
    rb_raise(rb_eArgError, 
	     "wrong number of arguments (%d for 1 or 2)", argc);
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), h);
  fp = rb_gsl_open_writefile(argv[0], &flag);
  if (argc == 2) {
    Check_Type(argv[1], T_STRING);
    status = FUNCTION(gsl_block,fprintf)(fp, h, STR2CSTR(argv[1]));
  } else {
    status = FUNCTION(gsl_block,fprintf)(fp, h, FORMAT_DEFAULT);
  }
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_gsl_block,printf)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_block) *h = NULL;
  int status;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), h);
  if (argc == 1) {
    Check_Type(argv[0], T_STRING);
    status = FUNCTION(gsl_block,fprintf)(stdout, h, STR2CSTR(argv[0]));
  } else {
    status = FUNCTION(gsl_block,fprintf)(stdout, h, FORMAT_DEFAULT);
  }
  return INT2FIX(status);
}

#undef FORMAT_DEFAULT

static VALUE FUNCTION(rb_gsl_block,fscanf)(VALUE obj, VALUE io)
{
  GSL_TYPE(gsl_block) *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), h);
  fp = rb_gsl_open_readfile(io, &flag);
  status = FUNCTION(gsl_block,fscanf)(fp, h);
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

#ifdef BASE_DOUBLE
#define SHOW_ELM 6
#define PRINTF_FORMAT "%4.3e "
#define TYPE2 double
#else
#define SHOW_ELM 15
#define PRINTF_FORMAT "%d "
#define TYPE2 int
#endif

static VALUE FUNCTION(rb_gsl_block,to_s)(VALUE obj)
{
  GSL_TYPE(gsl_block) *v = NULL;
  char buf[32];
  size_t i, n;
  VALUE str;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), v);
  str = rb_str_new2("[ ");
  n = v->size;
  if (rb_obj_is_kind_of(obj, cgsl_block_complex)) n *= 2;
  for (i = 0; i < n; i++) {
    sprintf(buf,  PRINTF_FORMAT, (TYPE2) v->data[i]);
    rb_str_cat(str, buf, strlen(buf));
    if (i == SHOW_ELM && i != v->size-1) {
      strcpy(buf, "... ");
      rb_str_cat(str, buf, strlen(buf));
      break;
    }
  }
  sprintf(buf, "]");
  rb_str_cat(str, buf, strlen(buf));
  return str;
}
#undef SHOW_ELM
#undef PRINTF_FORMAT

static VALUE FUNCTION(rb_gsl_block,inspect)(VALUE obj)
{
  VALUE str;
  char buf[64];
  sprintf(buf, "%s\n", rb_class2name(CLASS_OF(obj)));
  str = rb_str_new2(buf);
  return rb_str_concat(str, FUNCTION(rb_gsl_block,to_s)(obj));
}

#ifdef BASE_DOUBLE
#define C_TO_VALUE rb_float_new
#define NUMCONV NUM2DBL
#else
#define C_TO_VALUE INT2FIX
#define NUMCONV FIX2INT
#endif

void get_range_int_beg_en_n(VALUE range, int *beg, int *en, size_t *n, int *step);
static VALUE FUNCTION(rb_gsl_block,get)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_block) *b, *bnew;
  gsl_index *p;
  int beg, en, i, step;
  size_t n, j, k;

  Data_Get_Struct(obj, GSL_TYPE(gsl_block), b);
  switch (argc) {
  case 0:
    rb_raise(rb_eArgError, "too few arguments (%d for >= 1)", argc);
    break;
  case 1:
    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      i = FIX2INT(argv[0]);
      if (i < 0) j = b->size + i; else j = (size_t) i;
      return C_TO_VALUE(b->data[j]);
      break;
    case T_ARRAY:
      //      n = RARRAY(argv[0])->len;
      n = RARRAY_LEN(argv[0]);
      bnew = FUNCTION(gsl_block,alloc)(n);
      for (j = 0; j < n; j++) {
	i = FIX2INT(rb_ary_entry(argv[0], j));
	if (i < 0) k = b->size + i; else k = i;
	bnew->data[j] = b->data[k];
      }
      return Data_Wrap_Struct(GSL_TYPE(cgsl_block), 0, FUNCTION(gsl_block,free), bnew);
      break;
    default:
      if (PERMUTATION_P(argv[0])) {
	Data_Get_Struct(argv[0], gsl_index, p);
	bnew = FUNCTION(gsl_block,alloc)(p->size);
	for (j = 0; j < p->size; j++) bnew->data[j] = b->data[p->data[j]];
	return Data_Wrap_Struct(GSL_TYPE(cgsl_block), 0, FUNCTION(gsl_block,free), bnew);
      } else if (CLASS_OF(argv[0]) == rb_cRange) {
	get_range_int_beg_en_n(argv[0], &beg, &en, &n, &step);
	bnew = FUNCTION(gsl_block,alloc)(n);
	for (j = 0; j < n; j++) 
	  bnew->data[j] = b->data[beg+j];
	return Data_Wrap_Struct(GSL_TYPE(cgsl_block), 0, FUNCTION(gsl_block,free),
				bnew);
      } else {
	rb_raise(rb_eArgError, "wrong argument type %s (Fixnum, Array, or Range expected)", rb_class2name(CLASS_OF(argv[0])));
	break;
      }
    }
    break;
  default:
    bnew = FUNCTION(gsl_block,alloc)(argc);
    for (j = 0; (int) j < argc; j++) {
      i = FIX2INT(argv[j]);
      if (i < 0) k = b->size + i; else k = i;
      bnew->data[j] = b->data[k];
    }
    return Data_Wrap_Struct(GSL_TYPE(cgsl_block), 0, FUNCTION(gsl_block,free), bnew);


    break;
  }
  return Qnil;
}

static VALUE FUNCTION(rb_gsl_block,set)(VALUE obj, VALUE ii, VALUE xx)
{
  GSL_TYPE(gsl_block) *b;
  BASE x;
  size_t i;
  CHECK_FIXNUM(ii);
  i = FIX2INT(ii);
  x = (BASE) NUMCONV(xx);
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), b);
  b->data[i] = x;
  return obj;
}


static int FUNCTION(gsl_block,eq)(const GSL_TYPE(gsl_block) *a, 
				   const GSL_TYPE(gsl_block) *b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b->data[i];
    c->data[i] = (x > y || x < y) ? 0 : 1;
  }
  return 0;
}

static int FUNCTION(gsl_block,ne)(const GSL_TYPE(gsl_block) *a, 
				   const GSL_TYPE(gsl_block) *b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b->data[i];
    c->data[i] = (x > y || x < y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,gt)(const GSL_TYPE(gsl_block) *a, 
				   const GSL_TYPE(gsl_block) *b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b->data[i];
    c->data[i] =  (x > y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,ge)(const GSL_TYPE(gsl_block) *a, 
				   const GSL_TYPE(gsl_block) *b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b->data[i];
    c->data[i] =  (x >= y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,lt)(const GSL_TYPE(gsl_block) *a, 
				   const GSL_TYPE(gsl_block) *b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b->data[i];
    c->data[i] =  (x < y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,le)(const GSL_TYPE(gsl_block) *a, 
				   const GSL_TYPE(gsl_block) *b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b->data[i];
    c->data[i] =  (x <= y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,and)(const GSL_TYPE(gsl_block) *a, 
				   const GSL_TYPE(gsl_block) *b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b->data[i];
    c->data[i] = (x != 0 && y != 0) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,or)(const GSL_TYPE(gsl_block) *a, 
				   const GSL_TYPE(gsl_block) *b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b->data[i];
    c->data[i] = (x != 0 || y != 0) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,xor)(const GSL_TYPE(gsl_block) *a, 
				    const GSL_TYPE(gsl_block) *b,
				    gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != b->size) return -1;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b->data[i];
    c->data[i] = ((x != 0) == (y != 0)) ? 0 : 1;
  }
  return 0;
}

static int FUNCTION(gsl_block,eq2)(const GSL_TYPE(gsl_block) *a, 
				    BASE b, gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b;
    c->data[i] = (x > y || x < y) ? 0 : 1;
  }
  return 0;
}

static int FUNCTION(gsl_block,ne2)(const GSL_TYPE(gsl_block) *a, 
				   BASE b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b;
    c->data[i] = (x > y || x < y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,gt2)(const GSL_TYPE(gsl_block) *a, 
				   BASE b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b;
    c->data[i] =  (x > y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,ge2)(const GSL_TYPE(gsl_block) *a, 
				   BASE b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b;
    c->data[i] =  (x >= y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,lt2)(const GSL_TYPE(gsl_block) *a, 
				   BASE b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b;
    c->data[i] =  (x < y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,le2)(const GSL_TYPE(gsl_block) *a, 
				    BASE b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b;
    c->data[i] =  (x <= y) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,and2)(const GSL_TYPE(gsl_block) *a, 
				    BASE b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b;
    c->data[i] = (x != 0 && y != 0) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,or2)(const GSL_TYPE(gsl_block) *a, 
				    BASE b,
				   gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b;
    c->data[i] = (x != 0 || y != 0) ? 1 : 0;
  }
  return 0;
}

static int FUNCTION(gsl_block,xor2)(const GSL_TYPE(gsl_block) *a, 
				     BASE b,
				    gsl_block_uchar *c)
{
  size_t i;
  BASE x, y;
  if (a->size != c->size) return -2;
  for (i = 0; i < a->size; i++) {
    x = a->data[i];
    y = b;
    c->data[i] = ((x != 0) == (y != 0)) ? 0 : 1;
  }
  return 0;
}


static VALUE FUNCTION(rb_gsl_block,compare)(VALUE aa, VALUE bb,
					     int (*cmp)(const GSL_TYPE(gsl_block)*,
							const GSL_TYPE(gsl_block)*,
							gsl_block_uchar*),
					     int (*cmp2)(const GSL_TYPE(gsl_block)*,
							 BASE,
							 gsl_block_uchar*))
{
  GSL_TYPE(gsl_block) *a, *b;
  /*  gsl_block_int *c;*/
  gsl_block_uchar *c;
  // local variable "status" declared and set, but never used
  //int status;
  Data_Get_Struct(aa, GSL_TYPE(gsl_block), a);
  c = gsl_block_uchar_alloc(a->size);
  if (BL_P(bb)) {
    Data_Get_Struct(bb, GSL_TYPE(gsl_block), b);
    if (a->size != b->size) 
      rb_raise(rb_eRuntimeError, "Block size mismatch, %d and %d", (int) a->size, 
	       (int) b->size);
    /*status =*/ (*cmp)(a, b, c);
  } else {
    /*status =*/ (*cmp2)(a, NUMCONV(bb), c);
  }
  return Data_Wrap_Struct(cgsl_block_uchar, 0, gsl_block_uchar_free, c);
}

static VALUE FUNCTION(rb_gsl_block,eq)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_block,compare)(aa, bb, FUNCTION(gsl_block,eq),
					 FUNCTION(gsl_block,eq2));
}

static VALUE FUNCTION(rb_gsl_block,ne)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_block,compare)(aa, bb, FUNCTION(gsl_block,ne),
					 FUNCTION(gsl_block,ne2));
}

static VALUE FUNCTION(rb_gsl_block,gt)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_block,compare)(aa, bb, FUNCTION(gsl_block,gt),
					 FUNCTION(gsl_block,gt2));
}

static VALUE FUNCTION(rb_gsl_block,ge)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_block,compare)(aa, bb, FUNCTION(gsl_block,ge),
					 FUNCTION(gsl_block,ge2));
}

static VALUE FUNCTION(rb_gsl_block,lt)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_block,compare)(aa, bb, FUNCTION(gsl_block,lt),
					 FUNCTION(gsl_block,lt2));
}

static VALUE FUNCTION(rb_gsl_block,le)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_block,compare)(aa, bb, FUNCTION(gsl_block,le),
					 FUNCTION(gsl_block,le2));
}

static VALUE FUNCTION(rb_gsl_block,and)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_block,compare)(aa, bb, FUNCTION(gsl_block,and),
					 FUNCTION(gsl_block,and2));
}

static VALUE FUNCTION(rb_gsl_block,or)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_block,compare)(aa, bb, FUNCTION(gsl_block,or),
					 FUNCTION(gsl_block,or2));
}

static VALUE FUNCTION(rb_gsl_block,xor)(VALUE aa, VALUE bb)
{
  return FUNCTION(rb_gsl_block,compare)(aa, bb, FUNCTION(gsl_block,xor),
					 FUNCTION(gsl_block,xor2));
}

static VALUE FUNCTION(rb_gsl_block,not)(VALUE obj)
{
  GSL_TYPE(gsl_block) *v;
  gsl_block_uchar *vv;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), v);
  vv = gsl_block_uchar_alloc(v->size);
  for (i = 0; i < v->size; i++) vv->data[i] = (v->data[i] != 0) ? 0 : 1;
  return Data_Wrap_Struct(cgsl_block_uchar, 0, gsl_block_uchar_free, vv);
}

static VALUE FUNCTION(rb_gsl_block,any)(VALUE obj)
{
  GSL_TYPE(gsl_block) *v = NULL;
  size_t i;

  Data_Get_Struct(obj, GSL_TYPE(gsl_block), v);
  if (rb_block_given_p()) {
    for (i = 0; i < v->size; i++) {
      if (rb_yield(C_TO_VALUE(v->data[i]))) return INT2FIX(1);
    }
    return INT2FIX(0);
  } else {
    for (i = 0; i < v->size; i++) if (v->data[i]) return INT2FIX(1);
    return INT2FIX(0);
  }
}

static VALUE FUNCTION(rb_gsl_block,any2)(VALUE obj)
{
  GSL_TYPE(gsl_block) *v = NULL;
  size_t i;

  Data_Get_Struct(obj, GSL_TYPE(gsl_block), v);
  if (rb_block_given_p()) {
    for (i = 0; i < v->size; i++) {
      if (rb_yield(C_TO_VALUE(v->data[i]))) return Qtrue;
    }
    return Qfalse;
  } else {
    for (i = 0; i < v->size; i++) if (v->data[i]) return Qtrue;
    return Qfalse;
  }
}

static VALUE FUNCTION(rb_gsl_block,all)(VALUE obj)
{
  GSL_TYPE(gsl_block) *v;
  size_t i;

  Data_Get_Struct(obj, GSL_TYPE(gsl_block), v);
  if (rb_block_given_p()) {
    for (i = 0; i < v->size; i++) 
      if (!rb_yield(C_TO_VALUE(v->data[i]))) return Qfalse;
    return Qtrue;
  } else {
    for (i = 0; i < v->size; i++) 
      if (!v->data[i]) return Qfalse;
    return Qtrue;
  }
}

static VALUE FUNCTION(rb_gsl_block,none)(VALUE obj)
{
  GSL_TYPE(gsl_block) *v;
  size_t i;

  Data_Get_Struct(obj, GSL_TYPE(gsl_block), v);
  if (rb_block_given_p()) {
    for (i = 0; i < v->size; i++) 
      if (rb_yield(C_TO_VALUE(v->data[i]))) return Qfalse;
    return Qtrue;
  } else {
    for (i = 0; i < v->size; i++) 
      if (v->data[i]) return Qfalse;
    return Qtrue;
  }
}

static VALUE FUNCTION(rb_gsl_block,where)(VALUE obj)
{
  GSL_TYPE(gsl_block) *v;
  gsl_index *vv;
  gsl_block_uchar *btmp = NULL;
  size_t i, j, n = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), v);
  /* count true elements */
  if (rb_block_given_p()) {
    btmp = gsl_block_uchar_alloc(v->size);
    for (i = 0; i < v->size; i++) {
      if (rb_yield(C_TO_VALUE(v->data[i]))) {
	btmp->data[i] = 1;
	n++;
      } else {
	btmp->data[i] = 0;
      }
    }  /* for */
  } else { /* block is not given */
    for (i = 0; i < v->size; i++) { if (v->data[i]) n++; }
  }
  if (n == 0) {
    if (btmp) gsl_block_uchar_free(btmp);
    return Qnil;
  }
  vv = gsl_permutation_alloc(n);
  for (i = 0, j = 0; i < v->size; i++) {
    if ((!btmp && v->data[i]) || (btmp && btmp->data[i])) {
      vv->data[j++] = i;
    }
  }
  if (btmp) gsl_block_uchar_free(btmp);
  return Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, vv);
}

static VALUE FUNCTION(rb_gsl_block,where2)(VALUE obj)
{
  GSL_TYPE(gsl_block) *v;
  gsl_index *v1, *v2;
  gsl_block_uchar *btmp = NULL;
  VALUE vv1, vv2;
  size_t i, j, k, n = 0;

  Data_Get_Struct(obj, GSL_TYPE(gsl_block), v);
  if (rb_block_given_p()) {
    btmp = gsl_block_uchar_alloc(v->size);
    for (i = 0; i < v->size; i++) {
      if (rb_yield(C_TO_VALUE(v->data[i]))) {
	btmp->data[i] = 1;
	n++;
      } else {
	btmp->data[i] = 0;
      }
    } /* for */
  } else {  /* block is not given */
    for (i = 0; i < v->size; i++) { if (v->data[i]) n++; }
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
      if ((!btmp && v->data[i]) || (btmp && btmp->data[i])) v1->data[j++] = i;
      else v2->data[k++] = i;
    }
    vv1 = Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, v1);
    vv2 = Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, v2);
  }
  if (btmp) gsl_block_uchar_free(btmp);
  return rb_ary_new3(2, vv1, vv2);
}

static VALUE FUNCTION(rb_gsl_block,each)(VALUE obj)
{
  GSL_TYPE(gsl_block) *b = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), b);
  for (i = 0; i < b->size; i++) {
    rb_yield(C_TO_VALUE(b->data[i]));
  }
  return obj;
}

static VALUE FUNCTION(rb_gsl_block,each_index)(VALUE obj)
{
  GSL_TYPE(gsl_block) *b = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), b);
  for (i = 0; i < b->size; i++) {
    rb_yield(INT2FIX(i));
  }
  return obj;
}

static VALUE FUNCTION(rb_gsl_block,collect)(VALUE obj)
{
  GSL_TYPE(gsl_block) *b = NULL, *bnew;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), b);
  bnew = FUNCTION(gsl_block,alloc)(b->size);
  for (i = 0; i < b->size; i++) {
    bnew->data[i] = NUMCONV(rb_yield(C_TO_VALUE(b->data[i])));
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_block), 0, FUNCTION(gsl_block,free), bnew);
}

static VALUE FUNCTION(rb_gsl_block,collect_bang)(VALUE obj)
{
  GSL_TYPE(gsl_block) *b = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_block), b);
  for (i = 0; i < b->size; i++) {
    b->data[i] = NUMCONV(rb_yield(C_TO_VALUE(b->data[i])));
  }
  return obj;
}

#undef C_TO_VALUE
#undef NUMCONV
#undef TYPE2
#undef CHECK_BL
#undef BL_P

void FUNCTION(Init_gsl_block,init)(VALUE module)
{
  rb_define_singleton_method(GSL_TYPE(cgsl_block), "new", 
			     FUNCTION(rb_gsl_block,new), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_block), "alloc", 
			     FUNCTION(rb_gsl_block,new), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_block), "calloc", 
			     FUNCTION(rb_gsl_block,calloc), 1);

  rb_define_method(GSL_TYPE(cgsl_block), "size", FUNCTION(rb_gsl_block,size), 0);
  rb_define_alias(GSL_TYPE(cgsl_block), "length", "size");
  rb_define_method(GSL_TYPE(cgsl_block), "fwrite", FUNCTION(rb_gsl_block,fwrite), 1);
  rb_define_method(GSL_TYPE(cgsl_block), "fread", FUNCTION(rb_gsl_block,fread), 1);
  rb_define_method(GSL_TYPE(cgsl_block), "fprintf", FUNCTION(rb_gsl_block,fprintf), -1);
  rb_define_method(GSL_TYPE(cgsl_block), "printf", FUNCTION(rb_gsl_block,printf), -1);
  rb_define_method(GSL_TYPE(cgsl_block), "fscanf", FUNCTION(rb_gsl_block,fscanf), 1);

  rb_define_method(GSL_TYPE(cgsl_block), "inspect", FUNCTION(rb_gsl_block,inspect), 0);
  rb_define_method(GSL_TYPE(cgsl_block), "to_s", FUNCTION(rb_gsl_block,to_s), 0);
  rb_define_method(GSL_TYPE(cgsl_block), "get", FUNCTION(rb_gsl_block,get), -1);
  rb_define_alias(GSL_TYPE(cgsl_block), "[]", "get");
  rb_define_method(GSL_TYPE(cgsl_block), "set", FUNCTION(rb_gsl_block,set), 2);
  rb_define_alias(GSL_TYPE(cgsl_block), "[]=", "set");
  /*****/
  rb_define_method(GSL_TYPE(cgsl_block), "eq", FUNCTION(rb_gsl_block,eq), 1);
  rb_define_method(GSL_TYPE(cgsl_block), "ne", FUNCTION(rb_gsl_block,ne), 1);
  rb_define_method(GSL_TYPE(cgsl_block), "gt", FUNCTION(rb_gsl_block,gt), 1);
  rb_define_alias(GSL_TYPE(cgsl_block), ">", "gt");
  rb_define_method(GSL_TYPE(cgsl_block), "ge", FUNCTION(rb_gsl_block,ge), 1);
  rb_define_alias(GSL_TYPE(cgsl_block), ">=", "ge");
  rb_define_method(GSL_TYPE(cgsl_block), "lt", FUNCTION(rb_gsl_block,lt), 1);
  rb_define_alias(GSL_TYPE(cgsl_block), "<", "lt");
  rb_define_method(GSL_TYPE(cgsl_block), "le", FUNCTION(rb_gsl_block,le), 1);
  rb_define_alias(GSL_TYPE(cgsl_block), "<=", "le");

  rb_define_method(GSL_TYPE(cgsl_block), "and", FUNCTION(rb_gsl_block,and), 1);
  rb_define_method(GSL_TYPE(cgsl_block), "or", FUNCTION(rb_gsl_block,or), 1);
  rb_define_method(GSL_TYPE(cgsl_block), "xor", FUNCTION(rb_gsl_block,xor), 1);
  rb_define_method(GSL_TYPE(cgsl_block), "not", FUNCTION(rb_gsl_block,not), 0);

  rb_define_method(GSL_TYPE(cgsl_block), "all?", FUNCTION(rb_gsl_block,all), 0);
  rb_define_method(GSL_TYPE(cgsl_block), "none?", FUNCTION(rb_gsl_block,none), 0);
  rb_define_method(GSL_TYPE(cgsl_block), "any",
		   FUNCTION(rb_gsl_block,any), 0);
  rb_define_method(GSL_TYPE(cgsl_block), "any?",
		   FUNCTION(rb_gsl_block,any2), 0);

  rb_define_method(GSL_TYPE(cgsl_block), "where", FUNCTION(rb_gsl_block,where), 0);
  rb_define_method(GSL_TYPE(cgsl_block), "where2", FUNCTION(rb_gsl_block,where2), 0);

  rb_define_method(GSL_TYPE(cgsl_block), "each", FUNCTION(rb_gsl_block,each), 0);
  rb_define_method(GSL_TYPE(cgsl_block), "each_index", FUNCTION(rb_gsl_block,each_index), 0);
  rb_define_method(GSL_TYPE(cgsl_block), "collect", FUNCTION(rb_gsl_block,collect), 0);
  rb_define_method(GSL_TYPE(cgsl_block), "collect!", FUNCTION(rb_gsl_block,collect_bang), 0);
  rb_define_alias(GSL_TYPE(cgsl_block), "map", "collect");
  rb_define_alias(GSL_TYPE(cgsl_block), "map!", "collect!");
}
