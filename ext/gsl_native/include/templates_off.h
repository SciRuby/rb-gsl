#ifdef FUNCTION
#undef FUNCTION
#endif

#ifdef CONCAT4
#undef CONCAT4
#endif

#ifdef CONCAT4x
#undef CONCAT4x
#endif

#ifdef CONCAT3
#undef CONCAT3
#endif

#ifdef CONCAT3x
#undef CONCAT3x
#endif

#ifdef CONCAT2
#undef CONCAT2
#endif

#ifdef CONCAT2x
#undef CONCAT2x
#endif

#ifdef GSL_TYPE
#undef GSL_TYPE
#endif

#ifdef REAL_GSL_TYPE
#undef REAL_GSL_TYPE
#endif

#ifdef QUALIFIED_GSL_TYPE
#undef QUALIFIED_GSL_TYPE
#endif

#ifdef VIEW
#undef VIEW
#endif

#ifdef REAL_VIEW
#undef REAL_VIEW
#endif

#ifdef QUALIFIED_VIEW
#undef QUALIFIED_VIEW
#endif

#ifdef QUALIFIED_REAL_GSL_TYPE
#undef QUALIFIED_REAL_GSL_TYPE
#endif

#ifdef QUALIFIED_REAL_VIEW
#undef QUALIFIED_REAL_VIEW
#endif

#ifdef USES_LONGDOUBLE
#undef USES_LONGDOUBLE
#endif

#ifdef SHORT_REAL
#undef SHORT_REAL
#endif

#ifndef USE_QUALIFIER
#ifdef QUALIFIER
#undef QUALIFIER
#endif
#endif

#undef BASE
#undef BASE_EPSILON
#undef SHORT
#undef ATOMIC
#undef MULTIPLICITY
#undef IN_FORMAT
#undef OUT_FORMAT
#undef ATOMIC_IO
#undef ZERO
#undef ONE
#undef NAME
#undef STRING
#undef EXPAND

#ifdef RUBY_3
#undef memcpy
#define memcpy ruby_nonempty_memcpy
#endif