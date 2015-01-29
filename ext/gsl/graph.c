/*
  graph.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_graph.h"

void make_graphcommand(char *command, VALUE hash)
{
  VALUE val;
  if (TYPE(hash) == T_STRING) {
    sprintf(command, "graph -T X -g 3 %s", STR2CSTR(hash));
    return;
  }

  strcpy(command, "graph");
  if (TYPE(hash) != T_HASH) rb_raise(rb_eTypeError,
             "wrong argument type %s (Hash expected)",
             rb_class2name(CLASS_OF(hash)));
  if ((val = rb_hash_aref(hash, rb_str_new2("T"))) != Qnil)
    sprintf(command, "%s -T %s", command, STR2CSTR(val));
  else
    sprintf(command, "%s -T X", command);

  val = rb_hash_aref(hash, rb_str_new2("C"));
  if (val == Qtrue)
    sprintf(command, "%s -C", command);

  if ((val = rb_hash_aref(hash, rb_str_new2("g"))) != Qnil)
    sprintf(command, "%s -g %d", command, (int) FIX2INT(val));
  else
    sprintf(command, "%s -g 3", command);

  if ((val = rb_hash_aref(hash, rb_str_new2("B"))) == Qtrue)
    sprintf(command, "%s -B", command);
  if ((val = rb_hash_aref(hash, rb_str_new2("E"))) != Qnil)
    sprintf(command, "%s -E %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("f"))) != Qnil)
    sprintf(command, "%s -f %f", command, NUM2DBL(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("F"))) != Qnil)
    sprintf(command, "%s -F %s", command, STR2CSTR(val));

  if ((val = rb_hash_aref(hash, rb_str_new2("h"))) != Qnil)
    sprintf(command, "%s -h %f", command, NUM2DBL(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("k"))) != Qnil)
    sprintf(command, "%s -k %f", command, NUM2DBL(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("K"))) != Qnil)
    sprintf(command, "%s -K %d", command, (int) FIX2INT(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("l"))) != Qnil) {
    if (str_tail_grep(STR2CSTR(val), "xy") || str_tail_grep(STR2CSTR(val), "x/y"))
      sprintf(command, "%s -l x -l y", command);
    else
      sprintf(command, "%s -l %s", command, STR2CSTR(val));
  }

  if ((val = rb_hash_aref(hash, rb_str_new2("L"))) != Qnil)
    sprintf(command, "%s -L \"%s\"", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("N"))) != Qnil)
    sprintf(command, "%s -N %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("r"))) != Qnil)
    sprintf(command, "%s -r %f", command, NUM2DBL(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("R"))) != Qnil)
    sprintf(command, "%s -R %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("s"))) == Qtrue)
    sprintf(command, "%s -s", command);
  if ((val = rb_hash_aref(hash, rb_str_new2("t"))) == Qtrue)
    sprintf(command, "%s -t", command);
  if ((val = rb_hash_aref(hash, rb_str_new2("u"))) != Qnil)
    sprintf(command, "%s -u %f", command, NUM2DBL(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("w"))) != Qnil)
    sprintf(command, "%s -w %f", command, NUM2DBL(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("x"))) != Qnil)
    sprintf(command, "%s -x %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("X"))) != Qnil)
    sprintf(command, "%s -X \"%s\"", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("y"))) != Qnil)
    sprintf(command, "%s -y %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("Y"))) != Qnil)
    sprintf(command, "%s -Y \"%s\"", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("bg-color"))) != Qnil)
    sprintf(command, "%s --bg-color %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("bitmap-size"))) != Qnil)
    sprintf(command, "%s --bitmap-size %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("frame-color"))) != Qnil)
    sprintf(command, "%s --frame-color %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("frame-line-width"))) != Qnil)
    sprintf(command, "%s --frame-line-width %f", command, NUM2DBL(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("max-line-length"))) != Qnil)
    sprintf(command, "%s --max-line-length %f", command, NUM2DBL(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("page-size"))) != Qnil)
    sprintf(command, "%s --page-size %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("pen-colors"))) != Qnil)
    sprintf(command, "%s --pen-colors %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("rotation"))) != Qnil)
    sprintf(command, "%s --rotation %f", command, NUM2DBL(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("title-font-name"))) != Qnil)
    sprintf(command, "%s --title-font-name %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("title-font-size"))) != Qnil)
    sprintf(command, "%s --title-font-size %f", command, NUM2DBL(val));

  if ((val = rb_hash_aref(hash, rb_str_new2("toggle-rotate-y-label"))) == Qtrue)
    sprintf(command, "%s --toggle-rotate-y-label", command);

  if ((val = rb_hash_aref(hash, rb_str_new2("m"))) != Qnil)
    sprintf(command, "%s -m %d", command, (int) FIX2INT(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("S"))) != Qnil)
    sprintf(command, "%s -S %d", command, (int) FIX2INT(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("W"))) != Qnil)
    sprintf(command, "%s -W %f", command, NUM2DBL(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("q"))) != Qnil)
    sprintf(command, "%s -q %f", command, NUM2DBL(val));

  if ((val = rb_hash_aref(hash, rb_str_new2("symbol-font-name"))) != Qnil)
    sprintf(command, "%s --symbol-font-name %s", command, STR2CSTR(val));

  if ((val = rb_hash_aref(hash, rb_str_new2("reposition"))) != Qnil)
    sprintf(command, "%s --reposition %s", command, STR2CSTR(val));
  if ((val = rb_hash_aref(hash, rb_str_new2("blankout"))) != Qnil)
    sprintf(command, "%s --blankout %s", command, STR2CSTR(val));

  if ((val = rb_hash_aref(hash, rb_str_new2("O"))) == Qtrue)
    sprintf(command, "%s -O", command);
}

static void gsl_graph_init(gsl_graph *g);
gsl_graph* gsl_graph_new()
{
  gsl_graph *g = NULL;
  g = ALLOC(gsl_graph);
  gsl_graph_init(g);
  return g;
}

static void gsl_graph_init(gsl_graph *g)
{
  g->xdata = Qnil;
  g->ydata = Qnil;

  g->T = Qnil;
  g->E = Qnil;
  g->f = Qnil;
  g->F = Qnil;
  g->g = Qnil;
  g->h = Qnil;
  g->k = Qnil;
  g->K = Qnil;
  g->l = Qnil;
  g->L = Qnil;
  g->N = Qnil;
  g->r = Qnil;
  g->R = Qnil;
  g->u = Qnil;
  g->w = Qnil;
  g->x = Qnil;
  g->y = Qnil;
  g->X = Qnil;
  g->Y = Qnil;

  g->bg = Qnil;
  g->bitmap_size = Qnil;
  g->frame = Qnil;
  g->frame_line_width = Qnil;
  g->max_line_length = Qnil;
  g->page_size = Qnil;
  g->pen_colors = Qnil;
  g->rotation = Qnil;
  g->title_font_size = Qnil;
  g->title_font_name = Qnil;
  g->rotate_y_label = Qfalse;

  g->I = Qnil;
  g->m = Qnil;
  g->S = Qnil;
  g->W = Qnil;
  g->q = Qnil;
  g->symbol_font_name = Qnil;
  g->reposition = Qnil;
  g->blankout = Qnil;

  g->s = Qfalse;
  g->t = Qfalse;

  g->B = Qfalse;
  g->C = Qfalse;
  g->O = Qfalse;
}

static void gsl_graph_mark(gsl_graph *g)
{
  rb_gc_mark(g->xdata);
  rb_gc_mark(g->ydata);
  rb_gc_mark(g->T);
  rb_gc_mark(g->E);
  rb_gc_mark(g->f);
  rb_gc_mark(g->F);
  rb_gc_mark(g->g);
  rb_gc_mark(g->h);
  rb_gc_mark(g->k);
  rb_gc_mark(g->K);
  rb_gc_mark(g->l);
  rb_gc_mark(g->L);
  rb_gc_mark(g->N);
  rb_gc_mark(g->r);
  rb_gc_mark(g->R);
  rb_gc_mark(g->s);
  rb_gc_mark(g->t);
  rb_gc_mark(g->u);
  rb_gc_mark(g->w);
  rb_gc_mark(g->x);
  rb_gc_mark(g->y);
  rb_gc_mark(g->X);
  rb_gc_mark(g->Y);
  rb_gc_mark(g->bg);
  rb_gc_mark(g->bitmap_size);
  rb_gc_mark(g->frame);
  rb_gc_mark(g->frame_line_width);
  rb_gc_mark(g->max_line_length);
  rb_gc_mark(g->page_size);
  rb_gc_mark(g->pen_colors);
  rb_gc_mark(g->rotation);
  rb_gc_mark(g->title_font_name);
  rb_gc_mark(g->title_font_size);
  rb_gc_mark(g->rotate_y_label);
  rb_gc_mark(g->I);
  rb_gc_mark(g->B);
  rb_gc_mark(g->m);
  rb_gc_mark(g->S);
  rb_gc_mark(g->W);
  rb_gc_mark(g->q);
  rb_gc_mark(g->C);
  rb_gc_mark(g->symbol_font_name);
  rb_gc_mark(g->reposition);
  rb_gc_mark(g->blankout);
  rb_gc_mark(g->O);
}

void gsl_graph_free(gsl_graph *g)
{
  free((gsl_graph *) g);
}

static VALUE rb_gsl_graph_set_xdata(VALUE obj, VALUE xx);
static VALUE rb_gsl_graph_set_ydata(VALUE obj, VALUE yy);
static VALUE rb_gsl_graph_set_xydata(VALUE obj, VALUE xx, VALUE yy);
static VALUE rb_gsl_graph_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_graph *g = NULL;
  VALUE obj;
  g = gsl_graph_new();
  obj = Data_Wrap_Struct(klass, gsl_graph_mark, gsl_graph_free, g);
  switch (argc) {
  case 1:
    rb_gsl_graph_set_xdata(obj, argv[0]);
    break;
  case 2:
    rb_gsl_graph_set_xydata(obj, argv[0], argv[1]);
    break;
  }
  return obj;
}

static VALUE rb_gsl_graph_init(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  gsl_graph_init(g);
  return obj;
}

static VALUE rb_gsl_graph_xdata(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->xdata;
}

static VALUE rb_gsl_graph_ydata(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->ydata;
}

static VALUE rb_gsl_graph_xydata(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return rb_ary_new3(2, g->xdata, g->ydata);
}

static VALUE rb_gsl_graph_set_xdata(VALUE obj, VALUE xx)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  CHECK_VECTOR(xx);
  g->xdata = xx;
  return obj;
}

static VALUE rb_gsl_graph_set_ydata(VALUE obj, VALUE yy)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  CHECK_VECTOR(yy);
  g->ydata = yy;
  return obj;
}

static VALUE rb_gsl_graph_set_xydata(VALUE obj, VALUE xx, VALUE yy)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  CHECK_VECTOR(xx); CHECK_VECTOR(yy);
  g->xdata = xx;
  g->ydata = yy;
  return obj;
}

static VALUE rb_gsl_graph_set_T(VALUE obj, VALUE T)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  Check_Type(T, T_STRING);
  g->T = T;
  return T;
}

static VALUE rb_gsl_graph_T(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->T;
}

static VALUE rb_gsl_graph_set_E(VALUE obj, VALUE E)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  Check_Type(E, T_STRING);
  g->E = E;
  return E;
}

static VALUE rb_gsl_graph_E(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->E;
}

static VALUE rb_gsl_graph_set_f(VALUE obj, VALUE f)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->f = f;
  return f;
}

static VALUE rb_gsl_graph_f(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->f;
}

static VALUE rb_gsl_graph_set_F(VALUE obj, VALUE F)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->F = F;
  return F;
}

static VALUE rb_gsl_graph_F(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->F;
}

static VALUE rb_gsl_graph_set_g(VALUE obj, VALUE gg)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  Check_Type(gg, T_FIXNUM);
  g->g = gg;
  return gg;
}

static VALUE rb_gsl_graph_g(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->g;
}

static VALUE rb_gsl_graph_set_h(VALUE obj, VALUE h)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->h = h;
  return h;
}

static VALUE rb_gsl_graph_h(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->h;
}

static VALUE rb_gsl_graph_set_k(VALUE obj, VALUE k)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->k = k;
  return k;
}

static VALUE rb_gsl_graph_k(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->k;
}

static VALUE rb_gsl_graph_set_K(VALUE obj, VALUE K)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->K = K;
  return K;
}

static VALUE rb_gsl_graph_K(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->K;
}

static VALUE rb_gsl_graph_set_l(VALUE obj, VALUE l)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->l = l;
  return l;
}

static VALUE rb_gsl_graph_l(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->l;
}

static VALUE rb_gsl_graph_set_L(VALUE obj, VALUE L)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->L = L;
  return L;
}

static VALUE rb_gsl_graph_L(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->L;
}

static VALUE rb_gsl_graph_set_N(VALUE obj, VALUE N)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->N = N;
  return N;
}

static VALUE rb_gsl_graph_N(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->N;
}

static VALUE rb_gsl_graph_set_r(VALUE obj, VALUE r)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->r = r;
  return r;
}

static VALUE rb_gsl_graph_r(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->r;
}

static VALUE rb_gsl_graph_set_R(VALUE obj, VALUE R)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->R = R;
  return R;
}

static VALUE rb_gsl_graph_R(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->R;
}

static VALUE rb_gsl_graph_set_s(VALUE obj, VALUE s)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->s = s;
  return s;
}

static VALUE rb_gsl_graph_s(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->s;
}

static VALUE rb_gsl_graph_set_t(VALUE obj, VALUE t)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->t = t;
  return t;
}

static VALUE rb_gsl_graph_t(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->t;
}

static VALUE rb_gsl_graph_set_u(VALUE obj, VALUE u)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->u = u;
  return u;
}

static VALUE rb_gsl_graph_u(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->u;
}

static VALUE rb_gsl_graph_set_w(VALUE obj, VALUE w)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->w = w;
  return w;
}

static VALUE rb_gsl_graph_w(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->w;
}

static VALUE rb_gsl_graph_set_x(VALUE obj, VALUE x)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->x = x;
  return x;
}

static VALUE rb_gsl_graph_x(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->x;
}

static VALUE rb_gsl_graph_set_y(VALUE obj, VALUE y)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->y = y;
  return y;
}

static VALUE rb_gsl_graph_y(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->y;
}
static VALUE rb_gsl_graph_set_X(VALUE obj, VALUE X)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  Check_Type(X, T_STRING);
  g->X = X;
  return X;
}

static VALUE rb_gsl_graph_set_Y(VALUE obj, VALUE Y)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  Check_Type(Y, T_STRING);
  g->Y = Y;
  return Y;
}

static VALUE rb_gsl_graph_X(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->X;
}

static VALUE rb_gsl_graph_Y(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->Y;
}

static VALUE rb_gsl_graph_set_bg(VALUE obj, VALUE bg)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->bg = bg;
  return bg;
}

static VALUE rb_gsl_graph_bg(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->bg;
}

static VALUE rb_gsl_graph_set_bitmap_size(VALUE obj, VALUE bitmap_size)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->bitmap_size = bitmap_size;
  return bitmap_size;
}

static VALUE rb_gsl_graph_bitmap_size(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->bitmap_size;
}

static VALUE rb_gsl_graph_set_frame(VALUE obj, VALUE frame)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->frame = frame;
  return frame;
}

static VALUE rb_gsl_graph_frame(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->frame;
}

static VALUE rb_gsl_graph_set_frame_line_width(VALUE obj, VALUE frame_line_width)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->frame_line_width = frame_line_width;
  return frame_line_width;
}

static VALUE rb_gsl_graph_frame_line_width(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->frame_line_width;
}

static VALUE rb_gsl_graph_set_max_line_length(VALUE obj, VALUE max_line_length)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->max_line_length = max_line_length;
  return max_line_length;
}

static VALUE rb_gsl_graph_max_line_length(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->max_line_length;
}

static VALUE rb_gsl_graph_set_page_size(VALUE obj, VALUE page_size)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->page_size = page_size;
  return page_size;
}

static VALUE rb_gsl_graph_page_size(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->page_size;
}

static VALUE rb_gsl_graph_set_pen_colors(VALUE obj, VALUE pen_colors)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->pen_colors = pen_colors;
  return pen_colors;
}

static VALUE rb_gsl_graph_pen_colors(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->pen_colors;
}

static VALUE rb_gsl_graph_set_rotation(VALUE obj, VALUE rotation)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->rotation = rotation;
  return rotation;
}

static VALUE rb_gsl_graph_rotation(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->rotation;
}

static VALUE rb_gsl_graph_set_title_font_name(VALUE obj, VALUE title_font_name)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->title_font_name = title_font_name;
  return title_font_name;
}

static VALUE rb_gsl_graph_title_font_name(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->title_font_name;
}

static VALUE rb_gsl_graph_set_title_font_size(VALUE obj, VALUE title_font_size)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->title_font_size = title_font_size;
  return title_font_size;
}

static VALUE rb_gsl_graph_title_font_size(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->title_font_size;
}

static VALUE rb_gsl_graph_set_rotate_y_label(VALUE obj, VALUE rotate_y_label)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->rotate_y_label = rotate_y_label;
  return rotate_y_label;
}

static VALUE rb_gsl_graph_rotate_y_label(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->rotate_y_label;
}

static VALUE rb_gsl_graph_set_I(VALUE obj, VALUE I)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->I = I;
  return I;
}

static VALUE rb_gsl_graph_I(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->I;
}

static VALUE rb_gsl_graph_set_B(VALUE obj, VALUE B)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->B = B;
  return B;
}

static VALUE rb_gsl_graph_B(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->B;
}

static VALUE rb_gsl_graph_set_m(VALUE obj, VALUE m)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->m = m;
  return m;
}

static VALUE rb_gsl_graph_m(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->m;
}

static VALUE rb_gsl_graph_set_S(VALUE obj, VALUE S)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->S = S;
  return S;
}

static VALUE rb_gsl_graph_S(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->S;
}

static VALUE rb_gsl_graph_set_W(VALUE obj, VALUE W)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->W = W;
  return W;
}

static VALUE rb_gsl_graph_W(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->W;
}

static VALUE rb_gsl_graph_set_q(VALUE obj, VALUE q)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->q = q;
  return q;
}

static VALUE rb_gsl_graph_q(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->q;
}

static VALUE rb_gsl_graph_set_C(VALUE obj, VALUE C)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->C = C;
  return C;
}

static VALUE rb_gsl_graph_C(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->C;
}

static VALUE rb_gsl_graph_set_symbol_font_name(VALUE obj, VALUE symbol_font_name)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->symbol_font_name = symbol_font_name;
  return symbol_font_name;
}

static VALUE rb_gsl_graph_symbol_font_name(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->symbol_font_name;
}

static VALUE rb_gsl_graph_set_reposition(VALUE obj, VALUE r)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->reposition = r;
  return r;
}

static VALUE rb_gsl_graph_reposition(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->reposition;
}

static VALUE rb_gsl_graph_set_blankout(VALUE obj, VALUE r)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->blankout = r;
  return r;
}

static VALUE rb_gsl_graph_blankout(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->blankout;
}

static VALUE rb_gsl_graph_set_O(VALUE obj, VALUE O)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  g->O = O;
  return O;
}

static VALUE rb_gsl_graph_O(VALUE obj)
{
  gsl_graph *g = NULL;
  Data_Get_Struct(obj, gsl_graph, g);
  return g->O;
}

#ifdef HAVE_GNU_GRAPH
static void gsl_graph_set_command(gsl_graph *g, char *command)
{
  char str[256];
  size_t i, len;
  VALUE val;
  strcpy(command, "graph");
  if (g->T == Qnil)
    sprintf(command, "%s -T X", command);
  else
    sprintf(command, "%s -T %s", command, STR2CSTR(g->T));

  if (TYPE(g->E) == T_STRING) {
    strcpy(str, STR2CSTR(g->E));
    if (strcmp(str, "x") == 0)
      sprintf(command, "%s -E x", command);
    else if (strcmp(str, "y") == 0)
      sprintf(command, "%s -E y", command);
    else if (strcmp(str, "xy")*strcmp(str, "x/y") == 0)
      sprintf(command, "%s -E x -E y", command);
    else
      rb_raise(rb_eRuntimeError, "unrecognized -E option %s", str);
  }

  if (g->f != Qnil)
    sprintf(command, "%s -f %f", command, NUM2DBL(g->f));

  if (TYPE(g->F) == T_STRING)
    sprintf(command, "%s -F %s", command, STR2CSTR(g->F));

  if (TYPE(g->g) == T_FIXNUM)
    sprintf(command, "%s -g %d", command, (int) FIX2INT(g->g));

  if (g->h != Qnil)
    sprintf(command, "%s -h %f", command, NUM2DBL(g->h));

  if (g->k != Qnil)
    sprintf(command, "%s -k %f", command, NUM2DBL(g->k));

  if (TYPE(g->K) == T_FIXNUM)
    sprintf(command, "%s -K %d", command, (int) FIX2INT(g->K));

  if (TYPE(g->l) == T_STRING) {
    strcpy(str, STR2CSTR(g->l));
    if (strcmp(str, "x") == 0)
      sprintf(command, "%s -l x", command);
    else if (strcmp(str, "y") == 0)
      sprintf(command, "%s -l y", command);
    else if (strcmp(str, "xy")*strcmp(str, "x/y") == 0)
      sprintf(command, "%s -l x -l y", command);
    else
      rb_raise(rb_eRuntimeError, "unrecognized -l option %s", str);
  }

  if (TYPE(g->L) == T_STRING)
    sprintf(command, "%s -L \"%s\"", command, STR2CSTR(g->L));

  if (TYPE(g->N) == T_STRING) {
    strcpy(str, STR2CSTR(g->N));
    if (strcmp(str, "x") == 0)
      sprintf(command, "%s -N x", command);
    else if (strcmp(str, "y") == 0)
      sprintf(command, "%s -N y", command);
    else if (strcmp(str, "xy")*strcmp(str, "x/y") == 0)
      sprintf(command, "%s -N x -N y", command);
    else
      rb_raise(rb_eRuntimeError, "unrecognized -N option %s", str);
  }

  if (g->r != Qnil)
    sprintf(command, "%s -r %f", command, NUM2DBL(g->r));

  if (TYPE(g->R) == T_STRING) {
    strcpy(str, STR2CSTR(g->R));
    if (strcmp(str, "x") == 0)
      sprintf(command, "%s -R x", command);
    else if (strcmp(str, "y") == 0)
      sprintf(command, "%s -R y", command);
    else if (strcmp(str, "xy")*strcmp(str, "x/y") == 0)
      sprintf(command, "%s -R x -R y", command);
    else
      rb_raise(rb_eRuntimeError, "unrecognized -R option %s", str);
  }

  if (g->u != Qnil)
    sprintf(command, "%s -u %f", command, NUM2DBL(g->u));

  if (g->w != Qnil)
    sprintf(command, "%s -w %f", command, NUM2DBL(g->w));

  switch (TYPE(g->x)) {
  case T_STRING:
    sprintf(command, "%s -x %s", command, STR2CSTR(g->x));
    break;
  case T_ARRAY:
    sprintf(command, "%s -x", command);
    //    len = RARRAY(g->x)->len;
    len = RARRAY_LEN(g->x);
    for (i = 0; i < len; i++) {
      val = rb_ary_entry(g->x, i);
      Need_Float(val);
      sprintf(command, "%s %f", command, NUM2DBL(val));
    }
    break;
  default:
    /* do nothing */
    break;
  }

  switch (TYPE(g->y)) {
  case T_STRING:
    sprintf(command, "%s -y %s", command, STR2CSTR(g->y));
    break;
  case T_ARRAY:
    sprintf(command, "%s -y", command);
    //    len = RARRAY(g->y)->len;
    len = RARRAY_LEN(g->y);
    for (i = 0; i < len; i++) {
      val = rb_ary_entry(g->y, i);
      Need_Float(val);
      sprintf(command, "%s %f", command, NUM2DBL(val));
    }
    break;
  default:
    /* do nothing */
    break;
  }

  if (g->X != Qnil)
    sprintf(command, "%s -X \"%s\"", command, STR2CSTR(g->X));
  if (g->Y != Qnil)
    sprintf(command, "%s -Y \"%s\"", command, STR2CSTR(g->Y));

  if (TYPE(g->bg) == T_STRING)
    sprintf(command, "%s --bg-color %s", command, STR2CSTR(g->bg));

  if (TYPE(g->bitmap_size) == T_STRING)
    sprintf(command, "%s --bitmap-size %s", command, STR2CSTR(g->bitmap_size));

  if (TYPE(g->frame) == T_STRING)
    sprintf(command, "%s --frame-color %s", command, STR2CSTR(g->frame));

  if (g->frame_line_width != Qnil)
    sprintf(command, "%s --frame-line-width %f", command, NUM2DBL(g->frame_line_width));

 if (g->max_line_length != Qnil)
   sprintf(command, "%s --max_line_length %d", command,
     (int) FIX2INT(g->max_line_length));

  if (g->page_size != Qnil)
    sprintf(command, "%s --page-size %s", command, STR2CSTR(g->page_size));

  if (g->pen_colors != Qnil)
    sprintf(command, "%s --pen-colors %s", command, STR2CSTR(g->pen_colors));

  if (g->rotation != Qnil)
    sprintf(command, "%s --rotation %d", command, (int) FIX2INT(g->rotation));

  if (g->title_font_name != Qnil)
    sprintf(command, "%s --title-font-name %s", command, STR2CSTR(g->title_font_name));  if (g->title_font_size != Qnil)
    sprintf(command, "%s --title-font-size %f", command, NUM2DBL(g->title_font_size));

  if (g->rotate_y_label == Qtrue)
    sprintf(command, "%s --toggle-rotate-y-label", command);

  if (g->I != Qnil)
    sprintf(command, "%s -I %s", command, STR2CSTR(g->I));

  if (g->s == Qtrue)
    sprintf(command, "%s -s", command);
  if (g->t == Qtrue)
    sprintf(command, "%s -t", command);
  if (g->B == Qtrue)
    sprintf(command, "%s -B", command);

  if (g->m != Qnil)
    sprintf(command, "%s -m %d", command, (int) FIX2INT(g->m));

  switch (TYPE(g->S)) {
  case T_STRING:
    sprintf(command, "%s -S %s", command, STR2CSTR(g->S));
    break;
  case T_ARRAY:
    //    if (RARRAY(g->S)->len == 2)
    if (RARRAY_LEN(g->S) == 2)
      sprintf(command, "%s -S %d %f", command, (int) FIX2INT(rb_ary_entry(g->S, 0)),
        NUM2DBL(rb_ary_entry(g->S, 1)));
    break;
  default:
    /* do nothing */
    break;
  }

  if (g->W != Qnil)
    sprintf(command, "%s -W %f", command, NUM2DBL(g->W));

  if (g->q != Qnil)
    sprintf(command, "%s -q %f", command, NUM2DBL(g->q));

  if (g->C == Qtrue)
    sprintf(command, "%s -C", command);

  if (g->symbol_font_name != Qnil)
    sprintf(command, "%s --symbol_font_name %s", command, STR2CSTR(g->symbol_font_name));

  switch (TYPE(g->reposition)) {
  case T_STRING:
    sprintf(command, "%s --reposition %s", command, STR2CSTR(g->reposition));
    break;
  case T_ARRAY:
    sprintf(command, "%s --reposition", command);
    //    len =  RARRAY(g->reposition)->len;
    len =  RARRAY_LEN(g->reposition);
    for (i = 0; i <len; i++) {
      val = rb_ary_entry(g->reposition, i);
      Need_Float(val);
      sprintf(command, "%s %f", command, NUM2DBL(val));
    }
    break;
  default:
    /* do nothing */
    break;
  }

  if (g->blankout != Qnil)
    sprintf(command, "%s --blankout %f", command, NUM2DBL(g->blankout));

  if (g->O == Qtrue)
    sprintf(command, "%s -O", command);
}
#endif

static VALUE rb_gsl_graph_graph(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_GNU_GRAPH
  gsl_graph *g = NULL;
  gsl_histogram *h = NULL;
  gsl_vector *x = NULL, *y = NULL;
  size_t i, size;
  FILE *fp;
  char command[1024];
  Data_Get_Struct(obj, gsl_graph, g);

  gsl_graph_set_command(g, command);
  switch (argc) {
  case 3:
    Check_Type(argv[2], T_STRING);
    sprintf(command, "%s %s", command, STR2CSTR(argv[2]));
    /* no break */
  case 2:
    if (TYPE(argv[1]) == T_STRING) {
      sprintf(command, "%s %s", command, STR2CSTR(argv[1]));
    } else if (VECTOR_P(argv[1])) {
      g->ydata = argv[1];
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (Vector or String expected)",
         rb_class2name(CLASS_OF(argv[1])));
    }
    /* no break */
  case 1:
    if (TYPE(argv[0]) == T_STRING) {
      sprintf(command, "%s %s", command, STR2CSTR(argv[0]));
    } else if (VECTOR_P(argv[0])) {
      g->xdata = argv[0];
    } else if (HISTOGRAM_P(argv[0])) {
      Data_Get_Struct(argv[0], gsl_histogram, h);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (Vector or String expected)",
         rb_class2name(CLASS_OF(argv[0])));
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of argumeuts (%d for 1-3)", argc);
    break;
  }
  if (VECTOR_P(g->xdata)) Data_Get_Struct(g->xdata, gsl_vector, x);
  if (VECTOR_P(g->ydata)) Data_Get_Struct(g->ydata, gsl_vector, y);

  if (x == NULL && h == NULL)
    rb_raise(rb_eRuntimeError, "data is not given");

  if (h) size = h->n;
  else size = x->size;

  fp = popen(command, "w");
  if (fp == NULL) rb_raise(rb_eIOError, "GNU graph not found.");
  for (i = 0; i < size; i++) {
    if (h)
      fprintf(fp, "%g %g\n%g %g\n", h->range[i], h->bin[i], h->range[i+1], h->bin[i]);
    else if (y == NULL)
      fprintf(fp, "%d %g\n", (int) i, gsl_vector_get(x, i));
    else
      fprintf(fp, "%g %g\n", gsl_vector_get(x, i), gsl_vector_get(y, i));
  }
  fflush(fp);
  pclose(fp);
  fp = NULL;
  return Qtrue;
#else
  rb_raise(rb_eNoMethodError, "GNU plotutils required");
  return Qfalse;
#endif
}

static VALUE rb_gsl_graph_step(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_GNU_GRAPH
  gsl_graph *g = NULL;
  gsl_vector *x = NULL, *y = NULL;
  size_t i, size;
  FILE *fp;
  char command[1024];
  Data_Get_Struct(obj, gsl_graph, g);

  gsl_graph_set_command(g, command);
  switch (argc) {
  case 3:
    Check_Type(argv[2], T_STRING);
    sprintf(command, "%s %s", command, STR2CSTR(argv[2]));
    /* no break */
  case 2:
    if (TYPE(argv[1]) == T_STRING) {
      sprintf(command, "%s %s", command, STR2CSTR(argv[1]));
    } else if (VECTOR_P(argv[1])) {
      g->ydata = argv[1];
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (Vector or String expected)",
         rb_class2name(CLASS_OF(argv[1])));
    }
    /* no break */
  case 1:
    if (TYPE(argv[0]) == T_STRING) {
      sprintf(command, "%s %s", command, STR2CSTR(argv[0]));
    } else if (VECTOR_P(argv[0])) {
      g->xdata = argv[0];
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (Vector or String expected)",
         rb_class2name(CLASS_OF(argv[0])));
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of argumeuts (%d for 1-3)", argc);
    break;
  }
  if (VECTOR_P(g->xdata)) Data_Get_Struct(g->xdata, gsl_vector, x);
  if (VECTOR_P(g->ydata)) Data_Get_Struct(g->ydata, gsl_vector, y);

  if (x == NULL)
    rb_raise(rb_eRuntimeError, "data is not given");

  size = x->size;

  fp = popen(command, "w");
  if (fp == NULL) rb_raise(rb_eIOError, "GNU graph not found.");
  for (i = 0; i < size; i++) {
    if (y == NULL) {
      fprintf(fp, "%d %g\n%d %g\n", (int) i, gsl_vector_get(x, i),
        (int) (i+1), gsl_vector_get(x, i));
    } else {
      if (i != size-1)
        fprintf(fp, "%g %g\n%g %g\n", gsl_vector_get(x, i), gsl_vector_get(y, i),
          gsl_vector_get(x, i+1), gsl_vector_get(y, i));
      else
        fprintf(fp, "%g %g\n%g %g", gsl_vector_get(x, i), gsl_vector_get(y, i),
          2.0*gsl_vector_get(x, i)-gsl_vector_get(x, i-1), gsl_vector_get(y, i));
    }
  }
  fflush(fp);
  pclose(fp);
  fp = NULL;
  return Qtrue;
#else
  rb_raise(rb_eNoMethodError, "GNU plotutils required");
  return Qfalse;
#endif
}

void Init_gsl_graph(VALUE module)
{
  VALUE cgsl_graph;

  cgsl_graph = rb_define_class_under(module, "Graph", cGSL_Object);

  rb_define_singleton_method(cgsl_graph, "new", rb_gsl_graph_new, -1);
  rb_define_singleton_method(cgsl_graph, "alloc", rb_gsl_graph_new, -1);

  /*****/

  rb_define_method(cgsl_graph, "init", rb_gsl_graph_init, 0);

  rb_define_method(cgsl_graph, "set_xdata", rb_gsl_graph_set_xdata, 1);
  rb_define_method(cgsl_graph, "set_ydata", rb_gsl_graph_set_ydata, 1);
  rb_define_method(cgsl_graph, "set_xydata", rb_gsl_graph_set_xydata, 2);

  rb_define_method(cgsl_graph, "xdata", rb_gsl_graph_xdata, 0);
  rb_define_method(cgsl_graph, "ydata", rb_gsl_graph_ydata, 0);
  rb_define_method(cgsl_graph, "xydata", rb_gsl_graph_xydata, 0);

  rb_define_method(cgsl_graph, "graph", rb_gsl_graph_graph, -1);
  rb_define_alias(cgsl_graph, "draw", "graph");
  rb_define_alias(cgsl_graph, "plot", "graph");

  rb_define_method(cgsl_graph, "graph_step", rb_gsl_graph_step, -1);
  rb_define_alias(cgsl_graph, "step", "graph_step");

  /*****/
  rb_define_method(cgsl_graph, "set_T", rb_gsl_graph_set_T, 1);
  rb_define_alias(cgsl_graph, "T=", "set_T");
  rb_define_alias(cgsl_graph, "display_type=", "set_T");
  rb_define_method(cgsl_graph, "T", rb_gsl_graph_T, 0);
  rb_define_alias(cgsl_graph, "display_type", "T");

  rb_define_method(cgsl_graph, "set_E", rb_gsl_graph_set_E, 1);
  rb_define_alias(cgsl_graph, "E=", "set_E");
  rb_define_alias(cgsl_graph, "axis_end=", "set_E");
  rb_define_alias(cgsl_graph, "toggle_axis_end=", "set_E");
  rb_define_method(cgsl_graph, "E", rb_gsl_graph_E, 0);
  rb_define_alias(cgsl_graph, "axis_end", "E");
  rb_define_alias(cgsl_graph, "toggle_axis_end", "E");

  rb_define_method(cgsl_graph, "set_f", rb_gsl_graph_set_f, 1);
  rb_define_alias(cgsl_graph, "f=", "set_f");
  rb_define_alias(cgsl_graph, "font_size=", "set_f");
  rb_define_method(cgsl_graph, "f", rb_gsl_graph_f, 0);
  rb_define_alias(cgsl_graph, "font_size", "f");

  rb_define_method(cgsl_graph, "set_F", rb_gsl_graph_set_F, 1);
  rb_define_alias(cgsl_graph, "F=", "set_F");
  rb_define_alias(cgsl_graph, "font_name=", "set_F");
  rb_define_method(cgsl_graph, "F", rb_gsl_graph_F, 0);
  rb_define_alias(cgsl_graph, "font_name", "F");

  rb_define_method(cgsl_graph, "set_g", rb_gsl_graph_set_g, 1);
  rb_define_alias(cgsl_graph, "g=", "set_g");
  rb_define_alias(cgsl_graph, "grid_style=", "set_g");
  rb_define_method(cgsl_graph, "g", rb_gsl_graph_g, 0);
  rb_define_alias(cgsl_graph, "grid_style", "g");

  rb_define_method(cgsl_graph, "set_h", rb_gsl_graph_set_h, 1);
  rb_define_alias(cgsl_graph, "h=", "set_h");
  rb_define_alias(cgsl_graph, "height=", "set_h");
  rb_define_alias(cgsl_graph, "height_of_plot=", "set_h");
  rb_define_method(cgsl_graph, "h", rb_gsl_graph_h, 0);
  rb_define_alias(cgsl_graph, "height", "h");
  rb_define_alias(cgsl_graph, "height_of_plot", "h");

  rb_define_method(cgsl_graph, "set_k", rb_gsl_graph_set_k, 1);
  rb_define_alias(cgsl_graph, "k=", "set_k");
  rb_define_alias(cgsl_graph, "tick_size=", "set_k");
  rb_define_method(cgsl_graph, "k", rb_gsl_graph_k, 0);
  rb_define_alias(cgsl_graph, "tick_size", "k");

  rb_define_method(cgsl_graph, "set_K", rb_gsl_graph_set_K, 1);
  rb_define_alias(cgsl_graph, "K=", "set_K");
  rb_define_alias(cgsl_graph, "clip_mode=", "set_K");
  rb_define_method(cgsl_graph, "K", rb_gsl_graph_K, 0);
  rb_define_alias(cgsl_graph, "clip_mode", "K");

  rb_define_method(cgsl_graph, "set_l", rb_gsl_graph_set_l, 1);
  rb_define_alias(cgsl_graph, "l=", "set_l");
  rb_define_alias(cgsl_graph, "log_axis=", "set_l");
  rb_define_alias(cgsl_graph, "toggle_log_axis=", "set_l");
  rb_define_method(cgsl_graph, "l", rb_gsl_graph_l, 0);
  rb_define_alias(cgsl_graph, "log_axis", "l");
  rb_define_alias(cgsl_graph, "toggle_log_axis", "l");

  rb_define_method(cgsl_graph, "set_L", rb_gsl_graph_set_L, 1);
  rb_define_alias(cgsl_graph, "L=", "set_L");
  rb_define_alias(cgsl_graph, "top_label=", "set_L");
  rb_define_method(cgsl_graph, "L", rb_gsl_graph_L, 0);
  rb_define_alias(cgsl_graph, "top_label", "L");

  rb_define_method(cgsl_graph, "set_N", rb_gsl_graph_set_N, 1);
  rb_define_alias(cgsl_graph, "N=", "set_N");
  rb_define_alias(cgsl_graph, "no_tics=", "set_N");
  rb_define_alias(cgsl_graph, "toggle_no_tics=", "set_N");
  rb_define_method(cgsl_graph, "N", rb_gsl_graph_N, 0);
  rb_define_alias(cgsl_graph, "no_tics", "N");
  rb_define_alias(cgsl_graph, "toggle_no_tics", "N");

  rb_define_method(cgsl_graph, "set_r", rb_gsl_graph_set_r, 1);
  rb_define_alias(cgsl_graph, "r=", "set_r");
  rb_define_alias(cgsl_graph, "right_shift=", "set_r");
  rb_define_method(cgsl_graph, "r", rb_gsl_graph_r, 0);
  rb_define_alias(cgsl_graph, "right_shift", "r");

  rb_define_method(cgsl_graph, "set_R", rb_gsl_graph_set_R, 1);
  rb_define_alias(cgsl_graph, "R=", "set_R");
  rb_define_alias(cgsl_graph, "round-to-next-tick=", "set_R");
  rb_define_alias(cgsl_graph, "toggle_round-to-next-tick=", "set_R");
  rb_define_method(cgsl_graph, "R", rb_gsl_graph_R, 0);
  rb_define_alias(cgsl_graph, "round-to-next-tick", "R");
  rb_define_alias(cgsl_graph, "toggle_round-to-next-tick", "R");

  rb_define_method(cgsl_graph, "set_s", rb_gsl_graph_set_s, 1);
  rb_define_alias(cgsl_graph, "s=", "set_s");
  rb_define_alias(cgsl_graph, "save_screen=", "set_s");
  rb_define_method(cgsl_graph, "s", rb_gsl_graph_s, 0);
  rb_define_alias(cgsl_graph, "save_screen", "s");

  rb_define_method(cgsl_graph, "set_t", rb_gsl_graph_set_t, 1);
  rb_define_alias(cgsl_graph, "t=", "set_t");
  rb_define_alias(cgsl_graph, "transpose_axes=", "set_t");
  rb_define_alias(cgsl_graph, "toggle_transpose_axes=", "set_t");
  rb_define_method(cgsl_graph, "t", rb_gsl_graph_t, 0);
  rb_define_alias(cgsl_graph, "transpose_axes", "t");
  rb_define_alias(cgsl_graph, "toggle_transpose_axes", "t");

  rb_define_method(cgsl_graph, "set_u", rb_gsl_graph_set_u, 1);
  rb_define_alias(cgsl_graph, "u=", "set_u");
  rb_define_alias(cgsl_graph, "upward_shift=", "set_u");
  rb_define_method(cgsl_graph, "u", rb_gsl_graph_u, 0);
  rb_define_alias(cgsl_graph, "upward_shift", "u");

  rb_define_method(cgsl_graph, "set_w", rb_gsl_graph_set_w, 1);
  rb_define_alias(cgsl_graph, "w=", "set_w");
  rb_define_alias(cgsl_graph, "width=", "set_w");
  rb_define_alias(cgsl_graph, "width_of_plot=", "set_w");
  rb_define_method(cgsl_graph, "w", rb_gsl_graph_w, 0);
  rb_define_alias(cgsl_graph, "width", "w");
  rb_define_alias(cgsl_graph, "width_of_plot", "w");

  rb_define_method(cgsl_graph, "set_x", rb_gsl_graph_set_x, 1);
  rb_define_alias(cgsl_graph, "x=", "set_x");
  rb_define_alias(cgsl_graph, "x_limits=", "set_x");
  rb_define_method(cgsl_graph, "x", rb_gsl_graph_x, 0);
  rb_define_alias(cgsl_graph, "x_limits", "x");

  rb_define_method(cgsl_graph, "set_y", rb_gsl_graph_set_y, 1);
  rb_define_alias(cgsl_graph, "y=", "set_y");
  rb_define_alias(cgsl_graph, "y_limits=", "set_y");
  rb_define_method(cgsl_graph, "y", rb_gsl_graph_y, 0);
  rb_define_alias(cgsl_graph, "y_limits", "y");

  rb_define_method(cgsl_graph, "set_X", rb_gsl_graph_set_X, 1);
  rb_define_alias(cgsl_graph, "X=", "set_X");
  rb_define_alias(cgsl_graph, "x_label=", "set_X");
  rb_define_method(cgsl_graph, "X", rb_gsl_graph_X, 0);
  rb_define_alias(cgsl_graph, "x_label", "X");
  rb_define_method(cgsl_graph, "set_Y", rb_gsl_graph_set_Y, 1);
  rb_define_alias(cgsl_graph, "Y=", "set_Y");
  rb_define_alias(cgsl_graph, "y_label=", "set_Y");
  rb_define_method(cgsl_graph, "Y", rb_gsl_graph_Y, 0);
  rb_define_alias(cgsl_graph, "y_label", "Y");

  rb_define_method(cgsl_graph, "set_bg", rb_gsl_graph_set_bg, 1);
  rb_define_alias(cgsl_graph, "bg=", "set_bg");
  rb_define_alias(cgsl_graph, "bg_color=", "set_bg");
  rb_define_alias(cgsl_graph, "set_bg_color", "set_bg");
  rb_define_method(cgsl_graph, "bg", rb_gsl_graph_bg, 0);
  rb_define_alias(cgsl_graph, "bg_color", "bg");

  rb_define_method(cgsl_graph, "set_bitmap_size", rb_gsl_graph_set_bitmap_size, 1);
  rb_define_alias(cgsl_graph, "bitmap_size=", "set_bitmap_size");
  rb_define_alias(cgsl_graph, "set_bitmap_size", "set_bitmap_size");
  rb_define_method(cgsl_graph, "bitmap_size", rb_gsl_graph_bitmap_size, 0);

  rb_define_method(cgsl_graph, "set_frame", rb_gsl_graph_set_frame, 1);
  rb_define_alias(cgsl_graph, "frame=", "set_frame");
  rb_define_alias(cgsl_graph, "frame_color=", "set_frame");
  rb_define_alias(cgsl_graph, "set_frame_color", "set_frame");
  rb_define_method(cgsl_graph, "frame", rb_gsl_graph_frame, 0);
  rb_define_alias(cgsl_graph, "frame_color", "frame");

  rb_define_method(cgsl_graph, "set_frame_line_width", rb_gsl_graph_set_frame_line_width, 1);
  rb_define_alias(cgsl_graph, "frame_line_width=", "set_frame_line_width");
  rb_define_method(cgsl_graph, "frame_line_width", rb_gsl_graph_frame_line_width, 0);

  rb_define_method(cgsl_graph, "set_max_line_length", rb_gsl_graph_set_max_line_length, 1);
  rb_define_alias(cgsl_graph, "max_line_length=", "set_max_line_length");
  rb_define_method(cgsl_graph, "max_line_length", rb_gsl_graph_max_line_length, 0);

  rb_define_method(cgsl_graph, "set_page_size", rb_gsl_graph_set_page_size, 1);
  rb_define_alias(cgsl_graph, "page_size=", "set_page_size");
  rb_define_method(cgsl_graph, "page_size", rb_gsl_graph_page_size, 0);

  rb_define_method(cgsl_graph, "set_pen_colors", rb_gsl_graph_set_pen_colors, 1);
  rb_define_alias(cgsl_graph, "pen_colors=", "set_pen_colors");
  rb_define_alias(cgsl_graph, "pen=", "set_pen_colors");
  rb_define_method(cgsl_graph, "pen_colors", rb_gsl_graph_pen_colors, 0);
  rb_define_alias(cgsl_graph, "pen", "pen_colors");

  rb_define_method(cgsl_graph, "set_rotation", rb_gsl_graph_set_rotation, 1);
  rb_define_alias(cgsl_graph, "rotation=", "set_rotation");
  rb_define_method(cgsl_graph, "rotation", rb_gsl_graph_rotation, 0);

  rb_define_method(cgsl_graph, "set_title_font_name", rb_gsl_graph_set_title_font_name, 1);
  rb_define_alias(cgsl_graph, "title_font_name=", "set_title_font_name");
  rb_define_method(cgsl_graph, "title_font_name", rb_gsl_graph_title_font_name, 0);

  rb_define_method(cgsl_graph, "set_title_font_size", rb_gsl_graph_set_title_font_size, 1);
  rb_define_alias(cgsl_graph, "title_font_size=", "set_title_font_size");
  rb_define_method(cgsl_graph, "title_font_size", rb_gsl_graph_title_font_size, 0);

  rb_define_method(cgsl_graph, "set_rotate_y_label", rb_gsl_graph_set_rotate_y_label, 1);
  rb_define_alias(cgsl_graph, "rotate_y_label=", "set_rotate_y_label");
  rb_define_alias(cgsl_graph, "toggle_rotate_y_label=", "set_rotate_y_label");
  rb_define_method(cgsl_graph, "rotate_y_label", rb_gsl_graph_rotate_y_label, 0);
  rb_define_alias(cgsl_graph, "toggle_rotate_y_label", "rotate_y_label");

  rb_define_method(cgsl_graph, "set_B", rb_gsl_graph_set_B, 1);
  rb_define_alias(cgsl_graph, "B=", "set_B");
  rb_define_alias(cgsl_graph, "auto_dump=", "set_B");
  rb_define_alias(cgsl_graph, "toggle_auto_dump=", "set_B");
  rb_define_method(cgsl_graph, "B", rb_gsl_graph_B, 0);
  rb_define_alias(cgsl_graph, "auto_dump", "B");
  rb_define_alias(cgsl_graph, "toggle_auto_dump", "B");

  rb_define_method(cgsl_graph, "set_I", rb_gsl_graph_set_I, 1);
  rb_define_alias(cgsl_graph, "I=", "set_I");
  rb_define_alias(cgsl_graph, "input_format=", "set_I");
  rb_define_method(cgsl_graph, "I", rb_gsl_graph_I, 0);
  rb_define_alias(cgsl_graph, "input_format", "I");

  rb_define_method(cgsl_graph, "set_m", rb_gsl_graph_set_m, 1);
  rb_define_alias(cgsl_graph, "m=", "set_m");
  rb_define_alias(cgsl_graph, "line_mode=", "set_m");
  rb_define_method(cgsl_graph, "m", rb_gsl_graph_m, 0);
  rb_define_alias(cgsl_graph, "line_mode", "m");

  rb_define_method(cgsl_graph, "set_S", rb_gsl_graph_set_S, 1);
  rb_define_alias(cgsl_graph, "S=", "set_S");
  rb_define_alias(cgsl_graph, "symbol=", "set_S");
  rb_define_method(cgsl_graph, "S", rb_gsl_graph_S, 0);
  rb_define_alias(cgsl_graph, "symbol", "S");

  rb_define_method(cgsl_graph, "set_W", rb_gsl_graph_set_W, 1);
  rb_define_alias(cgsl_graph, "W=", "set_W");
  rb_define_alias(cgsl_graph, "line_width=", "set_W");
  rb_define_method(cgsl_graph, "W", rb_gsl_graph_W, 0);
  rb_define_alias(cgsl_graph, "line_width", "W");

  rb_define_method(cgsl_graph, "set_q", rb_gsl_graph_set_q, 1);
  rb_define_alias(cgsl_graph, "q=", "set_q");
  rb_define_alias(cgsl_graph, "fill_fraction=", "set_q");
  rb_define_method(cgsl_graph, "q", rb_gsl_graph_q, 0);
  rb_define_alias(cgsl_graph, "fill_fraction", "q");

  rb_define_method(cgsl_graph, "set_C", rb_gsl_graph_set_C, 1);
  rb_define_alias(cgsl_graph, "C=", "set_C");
  rb_define_alias(cgsl_graph, "use_color=", "set_C");
  rb_define_alias(cgsl_graph, "toggle_use_color=", "set_C");
  rb_define_method(cgsl_graph, "C", rb_gsl_graph_C, 0);
  rb_define_alias(cgsl_graph, "use_color", "C");
  rb_define_alias(cgsl_graph, "toggle_use_color", "C");

  rb_define_method(cgsl_graph, "set_symbol_font_name",
       rb_gsl_graph_set_symbol_font_name, 1);
  rb_define_alias(cgsl_graph, "symbol_font_name=", "set_symbol_font_name");
  rb_define_method(cgsl_graph, "symbol_font_name", rb_gsl_graph_symbol_font_name, 0);

  rb_define_method(cgsl_graph, "set_reposition",
       rb_gsl_graph_set_reposition, 1);
  rb_define_alias(cgsl_graph, "reposition=", "set_reposition");
  rb_define_method(cgsl_graph, "reposition", rb_gsl_graph_reposition, 0);

  rb_define_method(cgsl_graph, "set_blankout",
       rb_gsl_graph_set_blankout, 1);
  rb_define_alias(cgsl_graph, "blankout=", "set_blankout");
  rb_define_method(cgsl_graph, "blankout", rb_gsl_graph_blankout, 0);

  rb_define_method(cgsl_graph, "set_O", rb_gsl_graph_set_O, 1);
  rb_define_alias(cgsl_graph, "O=", "set_O");
  rb_define_alias(cgsl_graph, "portable_output=", "set_O");
  rb_define_method(cgsl_graph, "O", rb_gsl_graph_O, 0);
  rb_define_alias(cgsl_graph, "portable_output", "O");
}
