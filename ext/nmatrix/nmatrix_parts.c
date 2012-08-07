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
  YALE_STORAGE *result = create_merged_yale_storage(left, right);

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




// You must create the storage manually. NMatrix will clean it up when the matrix itself is destroyed.
NMATRIX* nm_create(int8_t stype, void* storage) {
  NMATRIX* mat = ALLOC(NMATRIX);

  mat->stype   = stype;
  mat->storage = storage;

  return mat;
}



/* Interprets cblas argument which could be any of false/:no_transpose, :transpose, or :complex_conjugate,
 * into an enum recognized by cblas.
 *
 * Called by nm_cblas_gemm -- basically inline.
 *
 */
static char gemm_op_sym(VALUE op) {
  if (op == false || rb_to_id(op) == nm_id_no_transpose) return CblasNoTrans;
  else if (rb_to_id(op) == nm_id_transpose) return CblasTrans;
  else if (rb_to_id(op) == nm_id_complex_conjugate) return CblasConjTrans;
  else rb_raise(rb_eArgError, "Expected false, :transpose, or :complex_conjugate");
  return CblasNoTrans;
}


/* Call any of the cblas_xgemm functions as directly as possible.
 *
 * The cblas_xgemm functions (dgemm, sgemm, cgemm, and zgemm) define the following operation:
 *
 *    C = alpha*op(A)*op(B) + beta*C
 *
 * where op(X) is one of <tt>op(X) = X</tt>, <tt>op(X) = X**T</tt>, or the complex conjugate of X.
 *
 * Note that this will only work for dense matrices that are of types :float32, :float64, :complex64, and :complex128.
 * Other types are not implemented in BLAS, and while they exist in NMatrix, this method is intended only to
 * expose the ultra-optimized ATLAS versions.
 *
 * == Arguments
 * See: http://www.netlib.org/blas/dgemm.f
 *
 * You probably don't want to call this function. Instead, why don't you try cblas_gemm, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_cblas_gemm(VALUE self,
                           VALUE trans_a, VALUE trans_b,
                           VALUE m, VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{
  char *pAlpha = ALLOCA_N(char, nm_sizeof[NM_DTYPE(c)]),
       *pBeta  = ALLOCA_N(char, nm_sizeof[NM_DTYPE(c)]);

  SetFuncs[NM_DTYPE(c)][NM_ROBJ](1, pAlpha, 0, &alpha, 0);
  SetFuncs[NM_DTYPE(c)][NM_ROBJ](1, pBeta, 0, &beta, 0);

  Gemm[NM_DTYPE(c)](gemm_op_sym(trans_a), gemm_op_sym(trans_b), FIX2INT(n), FIX2INT(m), FIX2INT(k), pAlpha,
                    ((DENSE_STORAGE*)(NM_STORAGE(b)))->elements, FIX2INT(ldb),
                    ((DENSE_STORAGE*)(NM_STORAGE(a)))->elements, FIX2INT(lda), pBeta,
                    ((DENSE_STORAGE*)(NM_STORAGE(c)))->elements, FIX2INT(ldc));

  return Qtrue;
}


/* Call any of the cblas_xgemv functions as directly as possible.
 *
 * The cblas_xgemv functions (dgemv, sgemv, cgemv, and zgemv) define the following operation:
 *
 *    y = alpha*op(A)*x + beta*y
 *
 * where op(A) is one of <tt>op(A) = A</tt>, <tt>op(A) = A**T</tt>, or the complex conjugate of A.
 *
 * Note that this will only work for dense matrices that are of types :float32, :float64, :complex64, and :complex128.
 * Other types are not implemented in BLAS, and while they exist in NMatrix, this method is intended only to
 * expose the ultra-optimized ATLAS versions.
 *
 * == Arguments
 * See: http://www.netlib.org/blas/dgemm.f
 *
 * You probably don't want to call this function. Instead, why don't you try cblas_gemv, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_cblas_gemv(VALUE self,
                           VALUE trans_a,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE x, VALUE incx,
                           VALUE beta,
                           VALUE y, VALUE incy)
{

  char *pAlpha = ALLOCA_N(char, nm_sizeof[NM_DTYPE(y)]),
       *pBeta  = ALLOCA_N(char, nm_sizeof[NM_DTYPE(y)]);

  SetFuncs[NM_DTYPE(y)][NM_ROBJ](1, pAlpha, 0, &alpha, 0);
  SetFuncs[NM_DTYPE(y)][NM_ROBJ](1, pBeta, 0, &beta, 0);

  Gemv[NM_ROBJ](gemm_op_sym(trans_a),
                FIX2INT(m), FIX2INT(n), pAlpha,
                ((DENSE_STORAGE*)(NM_STORAGE(a)))->elements, FIX2INT(lda),
                ((DENSE_STORAGE*)(NM_STORAGE(x)))->elements, FIX2INT(incx), pBeta,
                ((DENSE_STORAGE*)(NM_STORAGE(y)))->elements, FIX2INT(incy));

  return Qtrue;
}



/*
 * Get the size of a Yale matrix (the number of elements actually stored).
 *
 * For capacity (the maximum number of elements that can be stored without a resize), use capacity instead.
 */
static VALUE nm_yale_size(VALUE self) {
  VALUE sz;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  SetFuncs[NM_ROBJ][s->itype](1, &sz, 0, (YALE_SIZE_PTR((s), nm_sizeof[s->itype])), 0);
  return sz;
}


/*
 * Get the A array of a Yale matrix (which stores the diagonal and the LU portions of the matrix).
 */
static VALUE nm_yale_a(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*sz);

  SetFuncs[NM_ROBJ][s->dtype](sz, vals, nm_sizeof[NM_ROBJ], s->a, nm_sizeof[s->dtype]);
  ary = rb_ary_new4(sz, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


/*
 * Get the diagonal ("D") portion of the A array of a Yale matrix.
 */
static VALUE nm_yale_d(VALUE self) {
  y_size_t sz;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*s->shape[0]);

  SetFuncs[NM_ROBJ][s->dtype](s->shape[0], vals, nm_sizeof[NM_ROBJ], s->a, nm_sizeof[s->dtype]);
  ary = rb_ary_new4(s->shape[0], vals);

  return ary;
}


/*
 * Get the non-diagonal ("LU") portion of the A array of a Yale matrix.
 */
static VALUE nm_yale_lu(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*(s->capacity - s->shape[0]));

  SetFuncs[NM_ROBJ][s->dtype](sz - s->shape[0] - 1, vals, nm_sizeof[NM_ROBJ], (char*)(s->a) + (s->shape[0] + 1)*nm_sizeof[s->dtype], nm_sizeof[s->dtype]);
  ary = rb_ary_new4(sz - s->shape[0] - 1, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


/*
 * Get the IA portion of the IJA array of a Yale matrix. This gives the start and end positions of rows in the
 * JA and LU portions of the IJA and A arrays, respectively.
 */
static VALUE nm_yale_ia(VALUE self) {
  y_size_t sz;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*(s->shape[0]+1));

  SetFuncs[NM_ROBJ][s->itype](s->shape[0]+1, vals, nm_sizeof[NM_ROBJ], s->ija, nm_sizeof[s->itype]);
  ary = rb_ary_new4(s->shape[0]+1, vals);

  return ary;
}


/*
 * Get the JA portion of the IJA array of a Yale matrix. This gives the column indices for entries in corresponding
 * positions in the LU portion of the A array.
 */
static VALUE nm_yale_ja(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*(s->capacity - s->shape[0]));

  SetFuncs[NM_ROBJ][s->itype](sz - s->shape[0] - 1, vals, nm_sizeof[NM_ROBJ], (char*)(s->ija) + (s->shape[0] + 1)*nm_sizeof[s->itype], nm_sizeof[s->itype]);
  ary = rb_ary_new4(sz - s->shape[0] - 1, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


/*
 * Get the IJA array of a Yale matrix.
 */
static VALUE nm_yale_ija(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*s->capacity);

  SetFuncs[NM_ROBJ][s->itype](sz, vals, nm_sizeof[NM_ROBJ], s->ija, nm_sizeof[s->itype]);
  ary = rb_ary_new4(sz, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}

/*
 * Calculate the exact determinant of a dense matrix.
 *
 * Returns nil for dense matrices which are not square or rank other than 2.
 *
 * Note: Currently only implemented for 2x2 and 3x3 matrices.
 */
static VALUE nm_det_exact(VALUE self) {
  NMATRIX* m;
  void* result = ALLOCA_N(char, nm_sizeof[NM_DTYPE(self)]);
  VALUE ret;

  if (NM_STYPE(self) != S_DENSE) rb_raise(nm_eStorageTypeError, "can only calculate exact determinant for dense matrices");

  UnwrapNMatrix(self, m);
  if (m->storage->rank != 2 || m->storage->shape[0] != m->storage->shape[1]) return Qnil;

  // Calculate the determinant and then assign it to the return value
  DetExact[NM_DTYPE(self)](m->storage->shape[0], ((DENSE_STORAGE*)(m->storage))->elements, m->storage->shape[0], result);
  SetFuncs[NM_ROBJ][NM_DTYPE(self)](1, &ret, 0, result, 0);

  return ret;
}


// This is probably faster and smaller than writing an array of transpose functions. But if you want to see what it would look like,
// see transp.template.c (not the yale one).
//
// Note that this is a copy-transpose. In-place transpose is a whole different operation and bag of worms.
static void dense_transpose_generic(const unsigned int M, const unsigned int N, const char* A, const int lda, char* B, const int ldb, size_t dtype_size) {
  unsigned int i, j;

  for (i = 0; i < N; ++i) {
    for (j = 0; j < M; ++j) {
      memcpy(B + (i*ldb+j)*dtype_size, A + (j*lda+i)*dtype_size, dtype_size);
    }
  }
}


static NMATRIX* transpose_new_dense(NMATRIX* self_m, size_t* shape) {
  NMATRIX* result = nm_create(S_DENSE, create_dense_storage(self_m->storage->dtype, shape, 2, NULL, 0));

  dense_transpose_generic(self_m->storage->shape[0],
                          self_m->storage->shape[1],
                          ((DENSE_STORAGE*)(self_m->storage))->elements,
                          self_m->storage->shape[1],
                          ((DENSE_STORAGE*)(result->storage))->elements,
                          result->storage->shape[1],
                          nm_sizeof[self_m->storage->dtype]);
  return result;
}


static NMATRIX* transpose_new_yale(NMATRIX* self_m, size_t* shape) {
  YALE_PARAM A, B;
  NMATRIX* result;
  size_t sz;

  YaleGetSize(sz, (YALE_STORAGE*)(self_m->storage)); // size of new matrix is going to be size of old matrix
  result = nm_create(S_YALE, create_yale_storage(self_m->storage->dtype, shape, 2, sz));

  // TODO: Do we really need to initialize the whole thing? Or just the A portion?
  init_yale_storage((YALE_STORAGE*)(result->storage));

  A.ia = A.ja = ((YALE_STORAGE*)(self_m->storage))->ija;
  B.ia = B.ja = ((YALE_STORAGE*)(result->storage))->ija;
  A.a  = ((YALE_STORAGE*)(self_m->storage))->a;
  B.a  = ((YALE_STORAGE*)(result->storage))->a;
  A.diag = true;

  // call the appropriate function pointer
  Transp[ self_m->storage->dtype ][ ((YALE_STORAGE*)(self_m->storage))->itype ](shape[0], shape[1], A.ia, A.ja, A.a, A.diag, B.ia, B.ja, B.a, true);

  return result;
}


static NMATRIX* transpose_new_err(NMATRIX* self_m, size_t* shape) {
  free(shape);
  rb_raise(rb_eNotImpError, "no transpose written for this type");
  return self_m;
}

nm_transpose_t TransposeFuncs = {
  transpose_new_dense,
  transpose_new_err,
  transpose_new_yale
};

/*
 * Transform the matrix (in-place) to its complex conjugate. Only works on complex matrices.
 */
static VALUE nm_complex_conjugate_bang(VALUE self) {
  NMATRIX* m;
  void* elem;
  size_t sz, p;

  UnwrapNMatrix(self, m);

  if (m->stype == S_DENSE) sz = count_storage_max_elements(m->storage);
  else if (m->stype == S_YALE) YaleGetSize(sz, m->storage);
  else rb_raise(rb_eNotImpError, "please cast to yale or dense (complex) first");

  elem = m->storage->elements; // this gets A array or elements array from dense and yale

  // Walk through and negate the imaginary component
  if (NM_DTYPE(self) == NM_COMPLEX64) {
    for (p = 0; p < sz; ++p)      ((complex64*)elem)[p].i = -((complex64*)elem)[p].i;
  } else if (NM_DTYPE(self) == NM_COMPLEX128) {
    for (p = 0; p < sz; ++p)      ((complex128*)elem)[p].i = -((complex128*)elem)[p].i;
  } else {
    rb_raise(rb_eNotImpError, "can only calculate in-place complex conjugate on matrices of type :complex64 or :complex128");
  }

  return self;
}


/*
 * Create a transposed copy of this matrix.
 */
static VALUE nm_transpose_new(VALUE self) {
  NMATRIX *self_m, *result;
  size_t* shape   = ALLOC_N(size_t, 2);

  UnwrapNMatrix( self, self_m );

  // switch the dimensions
  shape[1] = self_m->storage->shape[0];
  shape[0] = self_m->storage->shape[1];

  result = TransposeFuncs[self_m->stype](self_m, shape);

  return Data_Wrap_Struct(cNMatrix, MarkFuncs[result->stype], nm_delete, result);
}


/*
 * Given a binary operation between types t1 and t2, what type will be returned?
 */
static VALUE nm_upcast(VALUE self, VALUE t1, VALUE t2) {
  int8_t dtype = Upcast[nm_dtypesymbol_to_dtype(t1)][nm_dtypesymbol_to_dtype(t2)];

  // The actual Upcast table returns NM_NONE if the types are unrecognized. If NM_NONE is the result, nil will be
  // returned instead.
  if (dtype == NM_NONE) return Qnil;

  return ID2SYM(rb_intern(nm_dtypestring[dtype]));
}



