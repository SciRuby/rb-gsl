#include "include/rb_gsl.h"

#ifdef HAVE_TAMU_ANOVA_TAMU_ANOVA_H
VALUE rb_tamu_anova_alloc(int argc, VALUE *argv, VALUE klass)
{
  gsl_vector *data;
  gsl_vector_long *factor;
  long I, J;
  struct tamu_anova_table *table;
  switch (argc) {
  case 3:
  case 4:
    Data_Get_Struct(argv[0], gsl_vector, data);
    Data_Get_Struct(argv[1], gsl_vector_long, factor);
    if (argc == 3) {
      I = data->size;
      J = NUM2INT(argv[2]);
    } else {
      I = NUM2INT(argv[2]);
      J = NUM2INT(argv[3]);
    }
    table = (struct tamu_anova_table *) malloc(sizeof(struct tamu_anova_table));
    *table = tamu_anova(data->data, factor->data, I, J);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 3 or 4)", argc);
    break;
  }
  return Data_Wrap_Struct(klass, 0, free, table);
}

VALUE rb_tamu_anova_printtable(VALUE *vTable)
{
  struct tamu_anova_table *table;
  Data_Get_Struct(vTable, struct tamu_anova_table, table);
  tamu_anova_printtable(*table);
  return Qtrue;
}

#endif

void Init_tamu_anova(VALUE module)
{
#ifdef HAVE_TAMU_ANOVA_TAMU_ANOVA_H
  VALUE mTAMU_ANOVA;
  VALUE cTable;

  mTAMU_ANOVA = rb_define_module_under(module, "TAMU_ANOVA");
  cTable = rb_define_class_under(mTAMU_ANOVA, "Table", cGSL_Object);

  rb_define_singleton_method(cTable, "alloc", rb_tamu_anova_alloc, -1);
  rb_define_singleton_method(cTable, "oneway", rb_tamu_anova_alloc, -1);

  rb_define_method(cTable, "print", rb_tamu_anova_printtable, 0);
#endif
}
