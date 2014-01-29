/*
  rb_gsl_graph.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or CHEBNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef ___RB_GSL_GRAPH_H___
#define ___RB_GSL_GRAPH_H___

#include "rb_gsl.h"

typedef struct __rb_gsl_plot {
  VALUE xdata, ydata;
  VALUE T;
  VALUE E;
  VALUE f;
  VALUE F;
  VALUE g;
  VALUE h;
  VALUE k;
  VALUE K;
  VALUE l;
  VALUE L;
  VALUE N;
  VALUE r;
  VALUE R;
  VALUE s;
  VALUE t;
  VALUE u;
  VALUE w;
  VALUE x;
  VALUE y;
  VALUE bg;
  VALUE bitmap_size;
  VALUE frame;
  VALUE frame_line_width;
  VALUE max_line_length;
  VALUE page_size;
  VALUE pen_colors;
  VALUE rotation;
  VALUE title_font_name;
  VALUE title_font_size;
  VALUE rotate_y_label;
  VALUE I;
  VALUE B;
  VALUE m;
  VALUE S;
  VALUE W;
  VALUE q;
  VALUE C;
  VALUE symbol_font_name;
  VALUE reposition;
  VALUE blankout;
  VALUE O;
  VALUE X, Y; 
  
} gsl_graph;

gsl_graph* gsl_graph_new();
void gsl_graph_free(gsl_graph *g);

#endif
