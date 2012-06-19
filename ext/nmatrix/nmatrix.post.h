// These have to come after enumerators
typedef void     (*nm_setfunc_t[NM_TYPES][NM_TYPES])(); // copy functions
typedef void     (*nm_incfunc_t[NM_TYPES])();           // increment functions
typedef void*    (*nm_stype_ref_t[S_TYPES])(STORAGE*, SLICE*);        // get/ref
typedef VALUE    (*nm_stype_ins_t[S_TYPES])(STORAGE*, SLICE*, VALUE); // insert
typedef STORAGE* (*nm_create_storage_t[S_TYPES])();
typedef STORAGE* (*nm_cast_copy_storage_t[S_TYPES])();
typedef STORAGE* (*nm_scast_copy_storage_t[S_TYPES][S_TYPES])();
typedef NMATRIX* (*nm_matrix_multiply_op_t[S_TYPES])();
typedef NMATRIX* (*nm_elementwise_binary_op_casted_t[S_TYPES])();
typedef int      (*nm_d_elementwise_binary_op_t[NM_TYPES])();
typedef int      (*nm_y_elementwise_binary_op_t[NM_TYPES][NM_INDEX_TYPES])();
typedef bool     (*nm_compare_t[S_TYPES])();
typedef void     (*nm_delete_t[S_TYPES])();
typedef void     (*nm_mark_t[S_TYPES])(void*);
typedef void     (*nm_gemm_t[NM_TYPES])();           // general matrix/matrix multiply
typedef void     (*nm_det_t[NM_TYPES])(const int, const void*, const int, void*);            // determinant
typedef NMATRIX* (*nm_transpose_t[S_TYPES])();
typedef void     (*nm_dense_transpose_t[NM_TYPES])(); // dense transpose
typedef void     (*nm_gemv_t[NM_TYPES])();           // general matrix/vector multiply
typedef void     (*nm_smmp_t[NM_TYPES][NM_INDEX_TYPES])(); // sparse (yale) multiply
typedef void     (*nm_smmp_transpose_t[NM_TYPES][NM_INDEX_TYPES])(y_size_t, y_size_t, YALE_PARAM, YALE_PARAM, bool); // sparse (yale) transpose

extern nm_setfunc_t SetFuncs;
extern nm_incfunc_t Increment;
extern ID nm_id_real, nm_id_imag;
extern ID nm_id_denom, nm_id_numer;
extern ID nm_id_mult, nm_id_multeq, nm_id_add;

//TODO: Auto-generate this
extern int (*EwDenseHom[15])(const void *, const void *, void *, const int, enum MathHomOps);
extern int (*EwDenseBool[15])(const void *, const void *, void *, const int, const enum MathBoolOps);
extern int (*EwDenseBit[15])(const void *, const void *, void *, const int, const enum MathBitOps);
extern int (*EwYaleHom[15][7])(const void *, const void *, void *, const void *, const void *, void *, const unsigned int, const unsigned int, const enum MathHomOps);
extern int (*EwYaleBool[15][7])(const void *, const void *, void *, const void *, const void *, void *, const unsigned int, const unsigned int, const enum MathBoolOps);
extern int (*EwYaleBit[15][7])(const void *, const void *, void *, const void *, const void *, void *, const unsigned int, const unsigned int, const enum MathBitOps);

const char *nm_dtypestring[] = {
  "none",
  "byte",
  "int8",
  "int16",
  "int32",
  "int64",
  "float32",
  "float64",
  "complex64",
  "complex128",
  "rational32",
  "rational64",
  "rational128",
  "object",
  "dtypes",
};

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

/* cblas.c */
DENSE_PARAM init_cblas_params_for_nm_multiply_matrix(int8_t dtype);
void cblas_bgemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_bgemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i8gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i8gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i16gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i16gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i32gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i32gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i64gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i64gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_sgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_sgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_dgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_dgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_cgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_cgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_zgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_zgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_r32gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_r32gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_r64gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_r64gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_r128gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_r128gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_vgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);


/* dense.c */
DENSE_STORAGE*  create_dense_storage(int8_t dtype, size_t* shape, size_t rank, void* elements, size_t elements_length);
void            delete_dense_storage_ref(DENSE_STORAGE* s);
void            delete_dense_storage(DENSE_STORAGE* s);
void            mark_dense_storage(void* s);
DENSE_STORAGE*  cast_copy_dense_storage(DENSE_STORAGE* rhs, int8_t new_dtype);

size_t          count_dense_storage_elements(const DENSE_STORAGE* s);
bool            dense_storage_eqeq(const DENSE_STORAGE*, const DENSE_STORAGE*);
bool            dense_is_symmetric(const DENSE_STORAGE*, int, bool);

size_t          dense_storage_pos(DENSE_STORAGE* s, SLICE* slice);
void*           dense_storage_get(DENSE_STORAGE* s, SLICE* slice);
void            dense_storage_set(DENSE_STORAGE* s, SLICE* slice, void* val);

/* list.c */
LIST_STORAGE*   create_list_storage(int8_t dtype, size_t* shape, size_t rank, void* init_val);
void            delete_list_storage(LIST_STORAGE* s);
void            mark_list_storage(void* s);
LIST_STORAGE*   cast_copy_list_storage(LIST_STORAGE* rhs, int8_t new_dtype);
size_t          count_storage_max_elements(const STORAGE*);

void*           list_storage_get(LIST_STORAGE* s, SLICE* slice);
void*           list_storage_insert(LIST_STORAGE* s, SLICE* slice, void* val);
void*           list_storage_remove(LIST_STORAGE* s, SLICE* slice);
bool            list_storage_eqeq(const LIST_STORAGE*, const LIST_STORAGE*);

/* yale.c */
void print_vectors(YALE_STORAGE* s);
YALE_STORAGE*   create_yale_storage(int8_t dtype, size_t* shape, size_t rank, size_t init_capacity);
YALE_STORAGE*   create_yale_storage_from_old_yale(int8_t dtype, size_t* shape, char* ia, char* ja, char* a, int8_t from_dtype, int8_t from_index_dtype);
void            init_yale_storage(YALE_STORAGE* s);
void            delete_yale_storage(YALE_STORAGE* s);
void            mark_yale_storage(void* s);
YALE_STORAGE*   cast_copy_yale_storage(YALE_STORAGE* rhs, int8_t new_dtype);
bool            yale_storage_eqeq(const YALE_STORAGE*, const YALE_STORAGE*);

void*           yale_storage_ref(YALE_STORAGE* s, SLICE* slice);
char            yale_storage_set(YALE_STORAGE* s, SLICE* slice, void* v);

YALE_STORAGE*   create_merged_yale_storage(const YALE_STORAGE*, const YALE_STORAGE*);

size_t          count_list_storage_nd_elements(const LIST_STORAGE*);
size_t          count_list_storage_elements(const LIST_STORAGE*);


/* stype casts */
DENSE_STORAGE* scast_copy_dense_yale(const YALE_STORAGE* rhs, int8_t l_dtype);
DENSE_STORAGE* scast_copy_dense_list(const LIST_STORAGE* rhs, int8_t l_dtype);
YALE_STORAGE* scast_copy_yale_dense(const DENSE_STORAGE* rhs, int8_t l_dtype);
YALE_STORAGE* scast_copy_yale_list(const LIST_STORAGE* rhs, int8_t l_dtype);
LIST_STORAGE* scast_copy_list_yale(const YALE_STORAGE* rhs, int8_t l_dtype);
LIST_STORAGE* scast_copy_list_dense(const DENSE_STORAGE* rhs, int8_t l_dtype);

/* nmatrix.c */
void cast_copy_value_single(void* to, const void* from, int8_t l_dtype, int8_t r_dtype);
int8_t nm_dtypestring_to_dtype(VALUE str);
int8_t nm_dtypesymbol_to_dtype(VALUE sym);
int8_t nm_stypestring_to_stype(VALUE str);
int8_t nm_stypesymbol_to_stype(VALUE sym);
int8_t nm_guess_dtype(VALUE v);
size_t* nm_interpret_shape_arg(VALUE arg, size_t* rank);
NMATRIX* nm_create(int8_t stype, void* storage);
void Init_nmatrix();

#endif
