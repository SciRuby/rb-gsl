#ifndef ___RB_GSL_HISTOGRAM3D_H___
#define ___RB_GSL_HISTOGRAM3D_H___

typedef struct {
  size_t nx, ny, nz;
  double * xrange ;
  double * yrange ;
  double * zrange ;
  double * bin ;
} mygsl_histogram3d ;

typedef struct {
  gsl_histogram2d h;
} mygsl_histogram3d_view;


mygsl_histogram3d* mygsl_histogram3d_alloc(const size_t nx, const size_t ny,
             const size_t nz);
void mygsl_histogram3d_free (mygsl_histogram3d * h);
mygsl_histogram3d* mygsl_histogram3d_calloc_uniform(const size_t nx,
                const size_t ny,
                const size_t nz,
                const double xmin,
                const double xmax,
                const double ymin,
                const double ymax,
                const double zmin,
                const double zmax);
mygsl_histogram3d* mygsl_histogram3d_calloc(const size_t nx,
              const size_t ny,
              const size_t nz);
int mygsl_histogram3d_set_ranges_uniform (mygsl_histogram3d * h,
            double xmin, double xmax,
            double ymin, double ymax,
            double zmin, double zmax);
int mygsl_histogram3d_set_ranges (mygsl_histogram3d * h,
          const double xrange[], size_t xsize,
          const double yrange[], size_t ysize,
          const double zrange[], size_t zsize);
int mygsl_histogram3d_memcpy(mygsl_histogram3d * dest, const mygsl_histogram3d * src);
mygsl_histogram3d* mygsl_histogram3d_clone(const mygsl_histogram3d * src);
int mygsl_histogram3d_fread(FILE * stream, mygsl_histogram3d * h);
int mygsl_histogram3d_fwrite(FILE * stream, const mygsl_histogram3d * h);
int mygsl_histogram3d_increment(mygsl_histogram3d * h, double x, double y, double z);
int mygsl_histogram3d_accumulate (mygsl_histogram3d * h,
          double x, double y, double z, double weight);
int mygsl_histogram3d_increment2(mygsl_histogram3d * h,
         double x, double y, double z);
int mygsl_histogram3d_accumulate2(mygsl_histogram3d * h,
          double x, double y, double z, double weight);
double mygsl_histogram3d_get (const mygsl_histogram3d * h, const size_t i,
            const size_t j, const size_t k);
int mygsl_histogram3d_get_xrange(const mygsl_histogram3d * h, const size_t i,
         double *xlower, double *xupper);
int mygsl_histogram3d_get_yrange(const mygsl_histogram3d * h, const size_t j,
         double *ylower, double *yupper);
int mygsl_histogram3d_get_zrange(const mygsl_histogram3d * h, const size_t k,
         double *zlower, double *zupper);
int mygsl_histogram3d_find (const mygsl_histogram3d * h,
          const double x, const double y, const double z,
          size_t * i, size_t * j, size_t *k);
gsl_histogram2d* mygsl_histogram3d_xyproject(const mygsl_histogram3d * h3,
               size_t kstart, size_t kend);
gsl_histogram2d* mygsl_histogram3d_xzproject(const mygsl_histogram3d * h3,
               size_t jstart, size_t jend);
gsl_histogram2d* mygsl_histogram3d_yzproject(const mygsl_histogram3d * h3,
               size_t istart, size_t iend);
int mygsl_histogram3d_scale(mygsl_histogram3d * h, double scale);
int mygsl_histogram3d_shift(mygsl_histogram3d * h, double shift);
double mygsl_histogram3d_xmax(const mygsl_histogram3d * h);
double mygsl_histogram3d_xmin(const mygsl_histogram3d * h);
double mygsl_histogram3d_ymax(const mygsl_histogram3d * h);
double mygsl_histogram3d_ymin(const mygsl_histogram3d * h);
double mygsl_histogram3d_zmax(const mygsl_histogram3d * h);
double mygsl_histogram3d_zmin(const mygsl_histogram3d * h);
double mygsl_histogram3d_max_val(const mygsl_histogram3d * h);
void mygsl_histogram3d_max_bin(const mygsl_histogram3d * h,
             size_t *imax_out, size_t *jmax_out, size_t *kmax_out);
double mygsl_histogram3d_min_val(const mygsl_histogram3d * h);
void mygsl_histogram3d_min_bin(const mygsl_histogram3d * h,
             size_t *imin_out, size_t *jmin_out, size_t *kmin_out);
double mygsl_histogram3d_sum (const mygsl_histogram3d * h);
double mygsl_histogram3d_xmean (const mygsl_histogram3d * h);
double mygsl_histogram3d_ymean (const mygsl_histogram3d * h);
double mygsl_histogram3d_zmean (const mygsl_histogram3d * h);
double mygsl_histogram3d_xsigma(const mygsl_histogram3d * h);
double mygsl_histogram3d_ysigma(const mygsl_histogram3d * h);
double mygsl_histogram3d_zsigma(const mygsl_histogram3d * h);
void mygsl_histogram3d_reset(mygsl_histogram3d * h);
int mygsl_histogram3d_equal_bins_p(const mygsl_histogram3d * h1,
           const mygsl_histogram3d * h2);
int mygsl_histogram3d_add(mygsl_histogram3d * h1, const mygsl_histogram3d * h2);
int mygsl_histogram3d_sub(mygsl_histogram3d * h1, const mygsl_histogram3d * h2);
int mygsl_histogram3d_mul(mygsl_histogram3d * h1, const mygsl_histogram3d * h2);
int mygsl_histogram3d_div(mygsl_histogram3d * h1, const mygsl_histogram3d * h2);

#endif
