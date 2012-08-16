const int8_t Upcast[15][15] = {
  { NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE },
  { NM_NONE, NM_BYTE, NM_INT8, NM_INT16, NM_INT32, NM_INT64, NM_FLOAT32, NM_FLOAT64, NM_COMPLEX64, NM_COMPLEX128, NM_RATIONAL32, NM_RATIONAL64, NM_RATIONAL128, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_INT8, NM_INT8, NM_INT16, NM_INT32, NM_INT64, NM_FLOAT32, NM_FLOAT64, NM_COMPLEX64, NM_COMPLEX128, NM_RATIONAL32, NM_RATIONAL64, NM_RATIONAL128, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_INT16, NM_INT16, NM_INT16, NM_INT32, NM_INT64, NM_FLOAT32, NM_FLOAT64, NM_COMPLEX64, NM_COMPLEX128, NM_RATIONAL32, NM_RATIONAL64, NM_RATIONAL128, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_INT32, NM_INT32, NM_INT32, NM_INT32, NM_INT64, NM_FLOAT32, NM_FLOAT64, NM_COMPLEX64, NM_COMPLEX128, NM_RATIONAL32, NM_RATIONAL64, NM_RATIONAL128, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_INT64, NM_INT64, NM_INT64, NM_INT64, NM_INT64, NM_FLOAT32, NM_FLOAT64, NM_COMPLEX64, NM_COMPLEX128, NM_RATIONAL32, NM_RATIONAL64, NM_RATIONAL128, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_FLOAT32, NM_FLOAT32, NM_FLOAT32, NM_FLOAT32, NM_FLOAT32, NM_FLOAT32, NM_FLOAT64, NM_COMPLEX64, NM_COMPLEX128, NM_FLOAT64, NM_FLOAT64, NM_FLOAT64, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_FLOAT64, NM_FLOAT64, NM_FLOAT64, NM_FLOAT64, NM_FLOAT64, NM_FLOAT64, NM_FLOAT64, NM_COMPLEX128, NM_COMPLEX128, NM_FLOAT64, NM_FLOAT64, NM_FLOAT64, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_COMPLEX64, NM_COMPLEX64, NM_COMPLEX64, NM_COMPLEX64, NM_COMPLEX64, NM_COMPLEX64, NM_COMPLEX128, NM_COMPLEX64, NM_COMPLEX128, NM_COMPLEX64, NM_COMPLEX64, NM_COMPLEX64, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_COMPLEX128, NM_COMPLEX128, NM_COMPLEX128, NM_COMPLEX128, NM_COMPLEX128, NM_COMPLEX128, NM_COMPLEX128, NM_COMPLEX128, NM_COMPLEX128, NM_COMPLEX128, NM_COMPLEX128, NM_COMPLEX128, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_RATIONAL32, NM_RATIONAL32, NM_RATIONAL32, NM_RATIONAL32, NM_RATIONAL32, NM_FLOAT64, NM_FLOAT64, NM_COMPLEX64, NM_COMPLEX128, NM_RATIONAL32, NM_RATIONAL64, NM_RATIONAL128, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_RATIONAL64, NM_RATIONAL64, NM_RATIONAL64, NM_RATIONAL64, NM_RATIONAL64, NM_FLOAT64, NM_FLOAT64, NM_COMPLEX64, NM_COMPLEX128, NM_RATIONAL64, NM_RATIONAL64, NM_RATIONAL128, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_RATIONAL128, NM_RATIONAL128, NM_RATIONAL128, NM_RATIONAL128, NM_RATIONAL128, NM_FLOAT64, NM_FLOAT64, NM_COMPLEX64, NM_COMPLEX128, NM_RATIONAL128, NM_RATIONAL128, NM_RATIONAL128, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_ROBJ, NM_ROBJ, NM_ROBJ, NM_ROBJ, NM_ROBJ, NM_ROBJ, NM_ROBJ, NM_ROBJ, NM_ROBJ, NM_ROBJ, NM_ROBJ, NM_ROBJ, NM_ROBJ, NM_NONE },
  { NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE, NM_NONE },
};


nm_delete_t DeleteFuncs = {
  delete_dense_storage,
  delete_list_storage,
  delete_yale_storage
};

nm_delete_t DeleteFuncsRef = {
  delete_dense_storage_ref,
  delete_list_storage,
  delete_yale_storage
};

nm_matrix_multiply_op_t CastedMultiplyFuncs = {
  multiply_matrix_dense_casted,
  multiply_matrix_list_casted,
  multiply_matrix_yale_casted
};


//static NMATRIX* elementwise_dense_casted(STORAGE_PAIR casted_storage, char op) {
static NMATRIX* ew_hom_dense_casted(STORAGE_PAIR casted_storage, char op) {
  DENSE_STORAGE *left  = (DENSE_STORAGE*)(casted_storage.left),
                *right = (DENSE_STORAGE*)(casted_storage.right),
                *result;

  // We can safely get dtype from the casted matrices; post-condition of binary_storage_cast_alloc is that dtype is the
  // same for left and right.
  size_t i;
  int8_t dtype = left->dtype;

  // Setup matrix shape for result
  size_t* shape = ALLOC_N(size_t, left->rank);
  for (i = 0; i < left->rank; ++i) shape[i] = left->shape[i];

  // Create result storage.
  result = create_dense_storage(dtype, shape, left->rank, NULL, 0);

  // Do the operation
  EwDenseHom[dtype](left->elements, right->elements, result->elements, count_dense_storage_elements(result), op);

  return nm_create(S_DENSE, result);
}


static NMATRIX* ew_hom_list_casted(STORAGE_PAIR casted_storage, char op) {
  rb_raise(rb_eNotImpError, "elementwise operations not implemented for list-of-list matrices");
  return NULL;
}


//static NMATRIX* elementwise_yale_casted(STORAGE_PAIR casted_storage, char op) {
static NMATRIX* ew_hom_yale_casted(STORAGE_PAIR casted_storage, enum MathHomOps op) {
  YALE_STORAGE *left  = (YALE_STORAGE*)(casted_storage.left),
               *right = (YALE_STORAGE*)(casted_storage.right);
  YALE_STORAGE *result = yale_storage_create_merged(left, right);

  fprintf(stderr, "result: %d, %d\n", result->dtype, result->itype);

  EwYaleHom[result->dtype][result->itype](left->a, right->a, result->a, left->ija, right->ija, result->ija, result->shape[0], result->shape[1], op);

  return nm_create(S_YALE, result);
}


nm_elementwise_binary_op_casted_t CastedElementwiseFuncs = {
  ew_hom_dense_casted,
  ew_hom_list_casted,
  ew_hom_yale_casted
};


nm_compare_t EqEqFuncs = {
  dense_storage_eqeq,
  list_storage_eqeq,
  yale_storage_eqeq
};


inline bool eqeq_generic(const void* x, const void* y, const int len, const int dtype_size) {
  return (!memcmp(x, y, len * dtype_size));
}

/*
 * TODO: Make this auto-generated in some future version of CSquare.
 *
 * element eqeq -- like memcmp but handles 0.0 == -0.0 for complex and floating
 * points. Second dimension is for hermitians -- complex conjugate. Use 0 for
 * regular equality and 1 for conjugate equality.
 *
 * FIXME: Each comparison here requires a function call, even if it is a
 * simple integer comparison.  They also can't be inlined because the
 * dtype isn't know until runtime.  This might be something that can be
 * solved with templating.
 */
bool (*ElemEqEq[NM_TYPES][2])(const void*, const void*, const int, const int) = {
  {NULL, NULL},
  {eqeq_generic, eqeq_generic}, // byte
  {eqeq_generic, eqeq_generic}, // int8
  {eqeq_generic, eqeq_generic}, // int16
  {eqeq_generic, eqeq_generic}, // int32
  {eqeq_generic, eqeq_generic}, // int64
  {eqeq_f32, eqeq_f32}, // float32
  {eqeq_f64, eqeq_f64}, // float64
  {eqeq_c64, conjeq_c64}, // complex64
  {eqeq_c128, conjeq_c128}, // complex128
  {eqeq_generic, eqeq_generic}, // rational32
  {eqeq_generic, eqeq_generic}, // rational64
  {eqeq_generic, eqeq_generic}, // rational128
  {eqeq_generic, eqeq_generic}  // Ruby object
};


static STORAGE* nm_dense_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, size_t init_val_len, VALUE self) {
  return (STORAGE*)(create_dense_storage(dtype, shape, rank, init_val, init_val_len));
}


static STORAGE* nm_list_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, size_t init_val_len, VALUE self) {
  if (init_val_len > 1) {
    rb_raise(rb_eArgError, "list storage needs initial size, not initial value\n");
    return NULL;
  }
  return (STORAGE*)(create_list_storage(dtype, shape, rank, init_val));
}


static STORAGE* nm_yale_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, size_t init_val_len, VALUE self) {
  YALE_STORAGE* s;

  if (init_val_len > 1) {
    rb_raise(rb_eArgError, "list storage needs initial size, not initial value\n");
    return NULL;
  }

  s = create_yale_storage(dtype, shape, rank, *(size_t*)init_val);
  init_yale_storage(s);
  free(init_val);

  if (!s) rb_raise(rb_eNoMemError, "Yale allocation failed");

  return (STORAGE*)(s);
  //return Data_Wrap_Struct(self, NULL, nm_delete, matrix);
}


nm_create_storage_t CreateFuncs = {
  nm_dense_new,
  nm_list_new,
  nm_yale_new
};


nm_cast_copy_storage_t CastCopyFuncs = {
  cast_copy_dense_storage,
  cast_copy_list_storage,
  cast_copy_yale_storage
};



nm_scast_copy_storage_t ScastCopyFuncs = {
  {cast_copy_dense_storage, scast_copy_dense_list, scast_copy_dense_yale},
  {scast_copy_list_dense, cast_copy_list_storage, scast_copy_list_yale},
  {scast_copy_yale_dense, scast_copy_yale_list, cast_copy_yale_storage}
};


nm_stype_slice_t GetFuncs = {
  dense_storage_get,
  list_storage_get,
  yale_storage_get
};


nm_stype_slice_t RefFuncs = {
  dense_storage_ref,
  list_storage_ref,
  yale_storage_ref
};


VALUE nm_dense_set(STORAGE* s, SLICE* slice, VALUE val) {
  void* v = ALLOCA_N(char, nm_sizeof[s->dtype]);
  SetFuncs[s->dtype][NM_ROBJ](1, v, 0, &val, 0);
  dense_storage_set( (DENSE_STORAGE*)s, slice, v );
  return val;
}


// Should work exactly the same as nm_dense_set.
VALUE nm_yale_set(STORAGE* s, SLICE* slice, VALUE val) {
  void* v = ALLOCA_N(char, nm_sizeof[s->dtype]);
  SetFuncs[s->dtype][NM_ROBJ](1, v, 0, &val, 0);
  yale_storage_set( (YALE_STORAGE*)s, slice, v );
  return val;
}


// TODO: Why can't you be more like your brothers, nm_dense_set and nm_yale_set?
VALUE nm_list_set(STORAGE* s, SLICE* slice, VALUE val) {
  void *v = ALLOC_N(char, nm_sizeof[s->dtype]), *rm;
  LIST_STORAGE* ls = (LIST_STORAGE*)s;

  //fprintf(stderr, "    create_val: %p\n", v);

  SetFuncs[s->dtype][NM_ROBJ](1, v, 0, &val, 0);

  if (!memcmp(ls->default_val, v, nm_sizeof[s->dtype])) {
    // User asked to insert default_value, which is actually node *removal*.
    // So let's do that instead.

    rm = list_storage_remove( ls, slice);

    //if (rm) fprintf(stderr, "    remove_val: %p\n", rm);

    if (rm) free(rm);
    return val;

  } else if (list_storage_insert( ls, slice, v ))    return val;
  return Qnil;
  // No need to free; the list keeps v.
}



nm_stype_ins_t InsFuncs = {
  nm_dense_set,
  nm_list_set,
  nm_yale_set,
};




/*
 * Equality operator. Returns a single true or false value indicating whether the matrices are equivalent.
 *
 * For elementwise, use == instead.
 *
 * This method will raise an exception if dimensions do not match.
 */
static VALUE nm_eqeq(VALUE left, VALUE right) {
  bool result;
  NMATRIX *l, *r;
  STORAGE_PAIR casted;

  CheckNMatrixType(left);
  CheckNMatrixType(right);

  UnwrapNMatrix(left, l);
  UnwrapNMatrix(right, r);

  if (l->stype != r->stype) //rb_raise(nm_eStorageTypeError, "wrong storage type");
    rb_raise(rb_eNotImpError, "comparison between different matrix stypes not yet implemented");

  casted = binary_storage_cast_alloc(l, r);

  result = EqEqFuncs[l->stype](casted.left, casted.right);

  // Free any casted-storage we created for the comparison.
  // TODO: Can we make the Ruby GC take care of this stuff now that we're using it?
  // If we did that, we night not have to re-create these every time, right? Or wrong? Need to do
  // more research.
  if (l->storage != casted.left)   DeleteFuncs[l->stype](casted.left);
  if (r->storage != casted.right)  DeleteFuncs[l->stype](casted.right);

  return result ? Qtrue : Qfalse;
}


static VALUE nm_elementwise(VALUE leftv, VALUE rightv, char op) {
  ///TODO: multiplication for non-dense and/or non-decimal matrices
  NMATRIX *result, *left, *right;
  STORAGE_PAIR casted;

  CheckNMatrixType(leftv);
  CheckNMatrixType(rightv);

  UnwrapNMatrix(rightv, right);
  UnwrapNMatrix(leftv, left);

  // Make sure both of our matrices are of the correct type.
  casted = binary_storage_cast_alloc(left, right);

  result = CastedElementwiseFuncs[left->stype](casted, op);

  // Free up temporary casted matrices
  if (left->storage != casted.left)   DeleteFuncs[left->stype](casted.left);
  if (right->storage != casted.right) DeleteFuncs[left->stype](casted.right);

  if (result) return Data_Wrap_Struct(cNMatrix, MarkFuncs[result->stype], nm_delete, result);
  return Qnil; // Only if we try to multiply list matrices should we return Qnil.
}


/*
 * Matrix element-wise addition.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 */
static VALUE nm_ew_add(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_MATHOP_ADD);
}

/*
 * Matrix element-wise subtraction.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 */
static VALUE nm_ew_subtract(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_MATHOP_SUB);
}

/*
 * Matrix element-wise multiplication.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * For dot product, use +dot+ instead.
 */
static VALUE nm_ew_multiply(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_MATHOP_MUL);
}

/*
 * Matrix element-wise division.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 */
static VALUE nm_ew_divide(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_MATHOP_DIV);
}


/*
 * Matrix element-wise modulus/norm.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 */
static VALUE nm_ew_mod(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_MATHOP_MOD);
}


/*
 * Matrix element-wise comparison (equality) operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_eqeq(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_MATHOP_EQEQ);
}

/*
 * Matrix element-wise less-than-or-equals operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_leq(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_MATHOP_LTE);
}


/*
 * Matrix element-wise greater-than-or-equals operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_geq(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_MATHOP_GTE);
}


/*
 * Matrix element-wise strictly-less-than operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_lt(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_MATHOP_LT);
}


/*
 * Matrix element-wise strictly-greater-than operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_gt(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_MATHOP_GT);
}


/*
 * Matrix element-wise inequality operator.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 *
 * Not available for list matrices. You should cast to a yale or dense matrix first.
 *
 * Note that the matrix returned will be of the same dtype as the upcast of the input matrices. If that's not what you
 * want, use +cast+.
 */
static VALUE nm_ew_neq(VALUE left, VALUE right) {
  return nm_elementwise(left, right, NM_MATHOP_NEQ);
}




