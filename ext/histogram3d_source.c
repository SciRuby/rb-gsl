/*
  histogram3d_source.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

/*
  This source code is created modifying gsl-1.5/histogram/ *2d.c
 */

#include "rb_gsl_histogram.h"

mygsl_histogram3d* mygsl_histogram3d_alloc(const size_t nx, const size_t ny,
					   const size_t nz)
{
  mygsl_histogram3d *h = NULL;
  if (nx == 0) GSL_ERROR_VAL ("histogram3d length nx must be positive integer",
			      GSL_EDOM, 0);
  if (ny == 0) GSL_ERROR_VAL ("histogram3d length ny must be positive integer",
			      GSL_EDOM, 0);
  if (nz == 0) GSL_ERROR_VAL ("histogram3d length nz must be positive integer",
			      GSL_EDOM, 0);
  h = (mygsl_histogram3d *) malloc(sizeof(mygsl_histogram3d));
  if (h == NULL) GSL_ERROR_VAL ("failed to allocate space for histogram3d struct",
				GSL_ENOMEM, 0);
  h->xrange = (double *) malloc ((nx + 1) * sizeof (double));
  if (h->xrange == 0) {
    free (h);         /* exception in constructor, avoid memory leak */
    GSL_ERROR_VAL ("failed to allocate space for histogram3d x ranges",
		   GSL_ENOMEM, 0);
  }
  h->yrange = (double *) malloc ((ny + 1) * sizeof (double));
  if (h->yrange == 0) {
    free (h->xrange);
    free (h);         /* exception in constructor, avoid memory leak */
    GSL_ERROR_VAL ("failed to allocate space for histogram3d y ranges",
		   GSL_ENOMEM, 0);
  }
  h->zrange = (double *) malloc ((nz + 1) * sizeof (double));
  if (h->zrange == 0) {
    free (h->xrange);
    free (h->yrange);
    free (h);         /* exception in constructor, avoid memory leak */
    GSL_ERROR_VAL ("failed to allocate space for histogram3d z ranges",
		   GSL_ENOMEM, 0);
  }
  h->bin = (double *) malloc (nx*ny*nz*sizeof (double));
  if (h->bin == 0) {
    free (h->xrange);
    free (h->yrange);
    free (h->zrange);
    free (h);         /* exception in constructor, avoid memory leak */
    GSL_ERROR_VAL ("failed to allocate space for histogram bins",
		   GSL_ENOMEM, 0);
  }
  h->nx = nx;
  h->ny = ny;
  h->nz = nz;
  return h;
}

mygsl_histogram3d* mygsl_histogram3d_calloc_uniform(const size_t nx, 
						    const size_t ny,
						    const size_t nz,
						    const double xmin,
						    const double xmax,
						    const double ymin,
						    const double ymax,
						    const double zmin,
						    const double zmax)
{
  mygsl_histogram3d *h;
  size_t i;
  h = mygsl_histogram3d_alloc(nx, ny, nz);
  for (i = 0; i < nx + 1; i++)
    h->xrange[i] = xmin + ((double) i / (double) nx) * (xmax - xmin);
  for (i = 0; i < ny + 1; i++)
    h->yrange[i] = ymin + ((double) i / (double) ny) * (ymax - ymin);
  for (i = 0; i < nz + 1; i++)
    h->zrange[i] = zmin + ((double) i / (double) nz) * (zmax - zmin);
  return h;
}

mygsl_histogram3d* mygsl_histogram3d_calloc(const size_t nx, 
					    const size_t ny,
					    const size_t nz)
{
  mygsl_histogram3d *h;
  size_t i;
  h = mygsl_histogram3d_alloc(nx, ny, nz);
  for (i = 0; i < nx + 1; i++)
    h->xrange[i] = i;
  for (i = 0; i < ny + 1; i++)
    h->yrange[i] = i;
  for (i = 0; i < nz + 1; i++)
    h->zrange[i] = i;
  for (i = 0; i < nx*ny*nz; i++)
    h->bin[i] = 0;
  return h;
}

int mygsl_histogram3d_set_ranges_uniform (mygsl_histogram3d * h, 
					  double xmin, double xmax,
					  double ymin, double ymax,
					  double zmin, double zmax)
{
  size_t i;
  const size_t nx = h->nx, ny = h->ny, nz = h->nz;

  if (xmin >= xmax) GSL_ERROR_VAL ("xmin must be less than xmax", GSL_EINVAL, 0);
  if (ymin >= ymax) GSL_ERROR_VAL ("ymin must be less than ymax", GSL_EINVAL, 0);
  if (zmin >= zmax) GSL_ERROR_VAL ("zmin must be less than zmax", GSL_EINVAL, 0);

  for (i = 0; i <= nx; i++)
    h->xrange[i] = xmin + ((double) i / (double) nx) * (xmax - xmin);
  for (i = 0; i <= ny; i++)
    h->yrange[i] = ymin + ((double) i / (double) ny) * (ymax - ymin);
  for (i = 0; i <= nz; i++)
    h->zrange[i] = zmin + ((double) i / (double) nz) * (zmax - zmin);
  for (i = 0; i < nx*ny*nz; i++) h->bin[i] = 0;
  return GSL_SUCCESS;
}

int mygsl_histogram3d_set_ranges(mygsl_histogram3d * h, 
				 const double xrange[], size_t xsize,
				 const double yrange[], size_t ysize,
				 const double zrange[], size_t zsize)
{
  size_t i;
  const size_t nx = h->nx, ny = h->ny, nz = h->nz;
  if (xsize != (nx + 1))
    GSL_ERROR_VAL ("size of xrange must match size of histogram", 
		   GSL_EINVAL, 0);
  if (ysize != (ny + 1))
    GSL_ERROR_VAL ("size of yrange must match size of histogram", 
		   GSL_EINVAL, 0);
  if (zsize != (nz + 1))
    GSL_ERROR_VAL ("size of yrange must match size of histogram", 
		   GSL_EINVAL, 0);
  memcpy(h->xrange, xrange, sizeof(double)*xsize);
  memcpy(h->yrange, yrange, sizeof(double)*ysize);
  memcpy(h->zrange, zrange, sizeof(double)*zsize);
  for (i = 0; i < nx*ny*nz; i++) h->bin[i] = 0;
  return GSL_SUCCESS;
}

void mygsl_histogram3d_free (mygsl_histogram3d * h)
{
  free (h->xrange);
  free (h->yrange);
  free (h->zrange);
  free (h->bin);
  free (h);
}

/*****/

int mygsl_histogram3d_memcpy(mygsl_histogram3d * dest, const mygsl_histogram3d * src)
{  
  size_t nx = src->nx;  
  size_t ny = src->ny;  
  size_t nz = src->nz;  
  if (dest->nx != src->nx || dest->ny != src->ny ||  dest->nz != src->nz) {
    GSL_ERROR ("histograms have different sizes, cannot copy",   
	       GSL_EINVAL);    
  }
  memcpy(dest->xrange, src->xrange, sizeof(double)*(nx+1));
  memcpy(dest->yrange, src->yrange, sizeof(double)*(ny+1));
  memcpy(dest->zrange, src->zrange, sizeof(double)*(nz+1));
  memcpy(dest->bin, src->bin, sizeof(double)*nx*ny*nz);
  return GSL_SUCCESS;
}

mygsl_histogram3d* mygsl_histogram3d_clone(const mygsl_histogram3d * src)
{
  mygsl_histogram3d *h;
  h = mygsl_histogram3d_alloc(src->nx, src->ny, src->nz);
  mygsl_histogram3d_memcpy(h, src);
  return h;
}

/*****/

int mygsl_histogram3d_fread(FILE * stream, mygsl_histogram3d * h)
{  
  int status = gsl_block_raw_fread (stream, h->xrange, h->nx + 1, 1);  
  if (status)    return status;  
  status = gsl_block_raw_fread (stream, h->yrange, h->ny + 1, 1);  
  if (status)    return status;  
  status = gsl_block_raw_fread (stream, h->zrange, h->nz + 1, 1);  
  if (status)    return status;  
  status = gsl_block_raw_fread (stream, h->bin, h->nx * h->ny * h->nz, 1);  
  return status;
}

int mygsl_histogram3d_fwrite(FILE * stream, const mygsl_histogram3d * h)
{  
  int status = gsl_block_raw_fwrite (stream, h->xrange, h->nx + 1, 1);  
  if (status)    return status;  
  status = gsl_block_raw_fwrite (stream, h->yrange, h->ny + 1, 1);  
  if (status)    return status;  
  status = gsl_block_raw_fwrite (stream, h->zrange, h->nz + 1, 1);  
  if (status)    return status;  
  status = gsl_block_raw_fwrite (stream, h->bin, h->nx * h->ny * h->nz, 1);  
  return status;
}

int mygsl_histogram3d_increment(mygsl_histogram3d * h, double x, double y, double z)
{  
  int status = mygsl_histogram3d_accumulate (h, x, y, z, 1.0);  
  return status;
}

int mygsl_find (const size_t n, const double range[], const double x, size_t * i);
int mygsl_find2d (const size_t nx, const double xrange[],
		  const size_t ny, const double yrange[],
		  const double x, const double y,
		  size_t * i, size_t * j);
int mygsl_find3d (const size_t nx, const double xrange[],
		  const size_t ny, const double yrange[],
		  const size_t nz, const double zrange[],
		  const double x, const double y, const double z,
		  size_t * i, size_t * j, size_t *k);

int mygsl_histogram3d_accumulate(mygsl_histogram3d * h, 
				 double x, double y, double z, double weight)
{  
  const size_t nx = h->nx;  
  const size_t ny = h->ny;  
  const size_t nz = h->nz;  
  size_t i = 0, j = 0, k = 0;  
  int status = mygsl_find3d (h->nx, h->xrange, h->ny, h->yrange, h->nz, h->zrange,
			     x, y, z, &i, &j, &k);  
  if (status) return GSL_EDOM; 
  if (i >= nx) GSL_ERROR ("index lies outside valid range of 0 .. nx - 1",   
			  GSL_ESANITY);    
  if (j >= ny) GSL_ERROR ("index lies outside valid range of 0 .. ny - 1",  
			  GSL_ESANITY);   
  if (k >= nz) GSL_ERROR ("index lies outside valid range of 0 .. nz - 1",  
			  GSL_ESANITY);   
  h->bin[i*ny*nz + j*nz + k] += weight;  
  return GSL_SUCCESS;
}

int mygsl_histogram3d_increment2(mygsl_histogram3d * h, 
				 double x, double y, double z)
{
  return mygsl_histogram3d_accumulate2(h, x, y, z, 1.0);
}

int mygsl_histogram3d_accumulate2(mygsl_histogram3d * h, 
				  double x, double y, double z, double weight)
{  
  const size_t nx = h->nx;  
  const size_t ny = h->ny;  
  const size_t nz = h->nz;  
  size_t i = 0, j = 0, k = 0;
  int status;
  if (x < h->xrange[0]) x = h->xrange[0] + 4*GSL_DBL_EPSILON;
  if (x > h->xrange[h->nx]) x = h->xrange[h->nx] - 4*GSL_DBL_EPSILON;
  if (y < h->yrange[0]) y = h->yrange[0] + 4*GSL_DBL_EPSILON;
  if (y > h->yrange[h->ny]) y = h->yrange[h->ny] - 4*GSL_DBL_EPSILON;
  if (z < h->zrange[0]) z = h->zrange[0] + 4*GSL_DBL_EPSILON;
  if (z > h->zrange[h->nz]) z = h->zrange[h->nz] - 4*GSL_DBL_EPSILON;
  status = mygsl_find3d (h->nx, h->xrange, h->ny, h->yrange, h->nz, h->zrange,
			 x, y, z, &i, &j, &k);  
  if (status) return GSL_EDOM; 
  if (i >= nx) GSL_ERROR ("index lies outside valid range of 0 .. nx - 1",   
			  GSL_ESANITY);    
  if (j >= ny) GSL_ERROR ("index lies outside valid range of 0 .. ny - 1",  
			  GSL_ESANITY);   
  if (k >= nz) GSL_ERROR ("index lies outside valid range of 0 .. nz - 1",  
			  GSL_ESANITY);   
  h->bin[i*ny*nz + j*nz + k] += weight;  
  return GSL_SUCCESS;
}

double mygsl_histogram3d_get(const mygsl_histogram3d * h, const size_t i, 
			     const size_t j, const size_t k)
{  
  const size_t nx = h->nx;  
  const size_t ny = h->ny;  
  const size_t nz = h->nz;  
  if (i >= nx) GSL_ERROR_VAL ("index i lies outside valid range of 0 .. nx - 1", 
			      GSL_EDOM, 0);   
  if (j >= ny) GSL_ERROR_VAL ("index j lies outside valid range of 0 .. ny - 1",    
			      GSL_EDOM, 0);    
  if (k >= nz) GSL_ERROR_VAL ("index k lies outside valid range of 0 .. nz - 1",    
			      GSL_EDOM, 0);    
  return h->bin[i*ny*nz + j*nz + k];
}

int mygsl_histogram3d_get_xrange(const mygsl_histogram3d * h, const size_t i,  
			       double *xlower, double *xupper)
{ 
  const size_t nx = h->nx;  
  if (i >= nx) GSL_ERROR ("index i lies outside valid range of 0 .. nx - 1", 
			  GSL_EDOM);    
  *xlower = h->xrange[i];  
  *xupper = h->xrange[i + 1];  
  return GSL_SUCCESS;
}

int mygsl_histogram3d_get_yrange(const mygsl_histogram3d * h, const size_t j,  
			       double *ylower, double *yupper)
{ 
  const size_t ny = h->ny;  
  if (j >= ny) GSL_ERROR ("index j lies outside valid range of 0 .. ny - 1", 
			  GSL_EDOM);    
  *ylower = h->yrange[j];  
  *yupper = h->yrange[j + 1];  
  return GSL_SUCCESS;
}

int mygsl_histogram3d_get_zrange(const mygsl_histogram3d * h, const size_t k,  
			       double *zlower, double *zupper)
{ 
  const size_t nz = h->nz;  
  if (k >= nz) GSL_ERROR ("index k lies outside valid range of 0 .. nz - 1", 
			  GSL_EDOM);    
  *zlower = h->zrange[k];  
  *zupper = h->zrange[k + 1];  
  return GSL_SUCCESS;
}

int mygsl_histogram3d_find (const mygsl_histogram3d * h,
			    const double x, const double y, const double z,
			    size_t * i, size_t * j, size_t *k)
{ 
  int status = mygsl_find(h->nx, h->xrange, x, i);  
  if (status) GSL_ERROR ("x not found in range of h", GSL_EDOM);    
  status = mygsl_find (h->ny, h->yrange, y, j);  
  if (status) GSL_ERROR ("y not found in range of h", GSL_EDOM);    
  status = mygsl_find (h->nz, h->zrange, z, k);  
  if (status) GSL_ERROR ("z not found in range of h", GSL_EDOM);    
  return GSL_SUCCESS;
}

gsl_histogram2d* mygsl_histogram3d_xyproject(const mygsl_histogram3d * h3,
					     size_t kstart, size_t kend)
{
  gsl_histogram2d *h2;
  double count;
  size_t i, j, k;
  h2 = gsl_histogram2d_calloc(h3->nx, h3->ny);
  gsl_histogram2d_set_ranges(h2, h3->xrange, h3->nx+1, h3->yrange, h3->ny+1);
  for (i = 0; i < h3->nx; i++) {
    for (j = 0; j < h3->ny; j++) {
      count = 0;
      for (k = kstart; k <= kend; k++) {
	if (k >= h3->nz) break;
	count += mygsl_histogram3d_get(h3, i, j, k);
      }
      h2->bin[i*h2->ny + j] = count;
    }
  }
  return h2;
}

gsl_histogram2d* mygsl_histogram3d_xzproject(const mygsl_histogram3d * h3,
					     size_t jstart, size_t jend)
{
  gsl_histogram2d *h2;
  double count;
  size_t i, j, k;
  h2 = gsl_histogram2d_calloc(h3->nx, h3->nz);
  gsl_histogram2d_set_ranges(h2, h3->xrange, h3->nx+1, h3->zrange, h3->nz+1);
  for (i = 0; i < h3->nx; i++) {
    for (k = 0; k < h3->nz; k++) {
      count = 0;
      for (j = jstart; j <= jend; j++) {
	if (j >= h3->ny) break;
	count += mygsl_histogram3d_get(h3, i, j, k);
      }
      h2->bin[i*h2->ny + k] = count;
    }
  }
  return h2;
}

gsl_histogram2d* mygsl_histogram3d_yzproject(const mygsl_histogram3d * h3,
					     size_t istart, size_t iend)
{
  gsl_histogram2d *h2;
  double count;
  size_t i, j, k;
  h2 = gsl_histogram2d_calloc(h3->ny, h3->nz);
  gsl_histogram2d_set_ranges(h2, h3->yrange, h3->ny+1, h3->zrange, h3->nz+1);
  for (j = 0; j < h3->ny; j++) {
    for (k = 0; k < h3->nz; k++) {
      count = 0;
      for (i = istart; i <= iend; i++) {
	if (i >= h3->nx) break;
	count += mygsl_histogram3d_get(h3, i, j, k);
      }
      h2->bin[j*h2->ny + k] = count;
    }
  }
  return h2;
}

int mygsl_histogram3d_scale(mygsl_histogram3d * h, double scale)
{
  size_t i, n;
  n = h->nx*h->ny*h->nz;
  for (i = 0; i < n; i++) h->bin[i] *= scale;
  return GSL_SUCCESS;
}

int mygsl_histogram3d_shift(mygsl_histogram3d * h, double shift)
{
  size_t i, n;
  n = h->nx*h->ny*h->nz;
  for (i = 0; i < n; i++) h->bin[i] += shift;
  return GSL_SUCCESS;
}

double mygsl_histogram3d_xmax(const mygsl_histogram3d * h)
{
  const int nx = h->nx;  
  return h->xrange[nx];
}
double mygsl_histogram3d_xmin(const mygsl_histogram3d * h)
{
  return h->xrange[0];
}

double mygsl_histogram3d_ymax(const mygsl_histogram3d * h)
{
  const int ny = h->ny;  
  return h->yrange[ny];
}
double mygsl_histogram3d_ymin(const mygsl_histogram3d * h)
{
  return h->yrange[0];
}

double mygsl_histogram3d_zmax(const mygsl_histogram3d * h)
{
  const int nz = h->nz;  
  return h->zrange[nz];
}

double mygsl_histogram3d_zmin(const mygsl_histogram3d * h)
{
  return h->zrange[0];
}

double mygsl_histogram3d_max_val(const mygsl_histogram3d * h)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  const size_t nz = h->nz;
  size_t i;
  double max = h->bin[0];
  for (i = 0; i < nx*ny*nz; i++) {
    if (h->bin[i] > max) max = h->bin[i];
  }
  return max;
}

void mygsl_histogram3d_max_bin(const mygsl_histogram3d * h, 
			       size_t *imax_out, size_t *jmax_out, size_t *kmax_out)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  const size_t nz = h->nz;
  size_t imax = 0, jmax = 0, kmax = 0, i, j, k;
  double max = h->bin[0], x;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
	x = h->bin[i * ny*nz + j*nz + k];
	if (x > max) {
	  max = x;
	  imax = i;
	  jmax = j;
	  kmax = k;
	}
      }
    }
  }
  *imax_out = imax;
  *jmax_out = jmax;
  *kmax_out = kmax;
}

double mygsl_histogram3d_min_val(const mygsl_histogram3d * h)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  const size_t nz = h->nz;
  size_t i;
  double min = h->bin[0];
  for (i = 0; i < nx*ny*nz; i++) {
    if (h->bin[i] < min) min = h->bin[i];
  }
  return min;
}

void mygsl_histogram3d_min_bin(const mygsl_histogram3d * h, 
			       size_t *imin_out, size_t *jmin_out, size_t *kmin_out)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  const size_t nz = h->nz;
  size_t imin = 0, jmin = 0, kmin = 0, i, j, k;
  double min = h->bin[0], x;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
	x = h->bin[i * ny*nz + j*nz + k];
	if (x < min) {
	  min = x;
	  imin = i;
	  jmin = j;
	  kmin = k;
	}
      }
    }
  }
  *imin_out = imin;
  *jmin_out = jmin;
  *kmin_out = kmin;
}

double mygsl_histogram3d_sum (const mygsl_histogram3d * h)
{
  const size_t n = h->nx*h->ny*h->nz;
  double sum = 0;
  size_t i = 0;
  while (i < n) sum += h->bin[i++];
  return sum;
}

double mygsl_histogram3d_xmean (const mygsl_histogram3d * h)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  const size_t nz = h->nz;
  size_t i, j, k;
  double wmean = 0;
  double W = 0;
  for (i = 0; i < nx; i++) {
    double xi = (h->xrange[i + 1] + h->xrange[i]) / 2.0;
    double wi = 0;
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
	double wijk = h->bin[i * ny *nz + j * nz + k];
	if (wijk > 0) wi += wijk;
      }
    }
    if (wi > 0) {
      W += wi;
      wmean += (xi - wmean) * (wi / W);
    }
  }
  return wmean;
}

double mygsl_histogram3d_ymean (const mygsl_histogram3d * h)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  const size_t nz = h->nz;
  size_t i, j, k;
  double wmean = 0;
  double W = 0;
  for (j = 0; j < ny; j++) {
    double yj = (h->yrange[j + 1] + h->yrange[j]) / 2.0;
    double wj = 0;
    for (i = 0; i < nx; i++) {
      for (k = 0; k < nz; k++) {
	double wijk = h->bin[i * ny *nz + j * nz + k];
	if (wijk > 0) wj += wijk;
      }
    }
    if (wj > 0) {
      W += wj;
      wmean += (yj - wmean) * (wj / W);
    }
  }
  return wmean;
}

double mygsl_histogram3d_zmean (const mygsl_histogram3d * h)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  const size_t nz = h->nz;
  size_t i, j, k;
  double wmean = 0;
  double W = 0;
  for (k = 0; k < nz; k++) {
    double zk = (h->zrange[k + 1] + h->zrange[k]) / 2.0;
    double wk = 0;
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
	double wijk = h->bin[i * ny *nz + j * nz + k];
	if (wijk > 0) wk += wijk;
      }
    }
    if (wk > 0) {
      W += wk;
      wmean += (zk - wmean) * (wk / W);
    }
  }
  return wmean;
}

double mygsl_histogram3d_xsigma(const mygsl_histogram3d * h)
{
  const double xmean = mygsl_histogram3d_xmean(h);
  const size_t nx = h->nx, ny = h->ny, nz = h->nz;
  size_t i, j, k;
  double wvariance = 0, W = 0;
  for (i = 0; i < nx; i++) {
    double xi = (h->xrange[i + 1] + h->xrange[i]) / 2 - xmean;
    double wi = 0;
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
	double wijk = h->bin[i * ny*nz + j*nz + k];
	if (wijk > 0) wi += wijk;
      }
    }
    if (wi > 0) {
      W += wi;
      wvariance += ((xi * xi) - wvariance) * (wi / W);
    }
  }
  return  sqrt(wvariance);
}

double mygsl_histogram3d_ysigma(const mygsl_histogram3d * h)
{
  const double ymean = mygsl_histogram3d_ymean(h);
  const size_t nx = h->nx, ny = h->ny, nz = h->nz;
  size_t i, j, k;
  double wvariance = 0, W = 0;
  for (j = 0; j < ny; j++) {
    double yj = (h->yrange[j + 1] + h->yrange[j]) / 2 - ymean;
    double wj = 0;
    for (i = 0; i < nx; i++) {
      for (k = 0; k < nz; k++) {
	double wjjk = h->bin[i * ny*nz + j*nz + k];
	if (wjjk > 0) wj += wjjk;
      }
    }
    if (wj > 0) {
      W += wj;
      wvariance += ((yj * yj) - wvariance) * (wj / W);
    }
  }
  return  sqrt(wvariance);
}

double mygsl_histogram3d_zsigma(const mygsl_histogram3d * h)
{
  const double zmean = mygsl_histogram3d_zmean(h);
  const size_t nx = h->nx, ny = h->ny, nz = h->nz;
  size_t i, j, k;
  double wvariance = 0, W = 0;
  for (k = 0; k < nz; k++) {
    double zk = (h->zrange[k + 1] + h->zrange[k]) / 2 - zmean;
    double wk = 0;
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
	double wijk = h->bin[i * ny*nz + j*nz + k];
	if (wijk > 0) wk += wijk;
      }
    }
    if (wk > 0) {
      W += wk;
      wvariance += ((zk * zk) - wvariance) * (wk / W);
    }
  }
  return  sqrt(wvariance);
}

void mygsl_histogram3d_reset(mygsl_histogram3d * h)
{
  size_t i;
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  const size_t nz = h->nz;

  for (i = 0; i < nx * ny * nz; i++)
    {
      h->bin[i] = 0;
    }
}

int mygsl_histogram3d_equal_bins_p(const mygsl_histogram3d * h1,
				   const mygsl_histogram3d * h2)
{
  size_t i;
  if ((h1->nx != h2->nx) || (h1->ny != h2->ny) || (h1->nz != h2->nz)) return 0;
  for (i = 0; i <= h1->nx; i++) 
    if (h1->xrange[i] != h2->xrange[i]) return 0;
  for (i = 0; i <= h1->ny; i++)
    if (h1->yrange[i] != h2->yrange[i]) return 0;
  for (i = 0; i <= h1->nz; i++)
    if (h1->zrange[i] != h2->zrange[i]) return 0;
  return 1;
}

int mygsl_histogram3d_add(mygsl_histogram3d * h1, const mygsl_histogram3d * h2)
{
  size_t i;
  if (!mygsl_histogram3d_equal_bins_p(h1, h2))
    GSL_ERROR ("histograms have different binning", GSL_EINVAL);
  for (i = 0; i < (h1->nx) * (h1->ny) * (h1->nz); i++)
    h1->bin[i] += h2->bin[i];
  return GSL_SUCCESS;
}

int mygsl_histogram3d_sub(mygsl_histogram3d * h1, const mygsl_histogram3d * h2)
{
  size_t i;
  if (!mygsl_histogram3d_equal_bins_p(h1, h2))
    GSL_ERROR ("histograms have different binning", GSL_EINVAL);
  for (i = 0; i < (h1->nx) * (h1->ny) * (h1->nz); i++)
    h1->bin[i] -= h2->bin[i];
  return GSL_SUCCESS;
}

int mygsl_histogram3d_mul(mygsl_histogram3d * h1, const mygsl_histogram3d * h2)
{
  size_t i;
  if (!mygsl_histogram3d_equal_bins_p(h1, h2))
    GSL_ERROR ("histograms have different binning", GSL_EINVAL);
  for (i = 0; i < (h1->nx) * (h1->ny) * (h1->nz); i++)
    h1->bin[i] *= h2->bin[i];
  return GSL_SUCCESS;
}

int mygsl_histogram3d_div(mygsl_histogram3d * h1, const mygsl_histogram3d * h2)
{
  size_t i;
  if (!mygsl_histogram3d_equal_bins_p(h1, h2))
    GSL_ERROR ("histograms have different binning", GSL_EINVAL);
  for (i = 0; i < (h1->nx) * (h1->ny) * (h1->nz); i++)
    h1->bin[i] /= h2->bin[i];
  return GSL_SUCCESS;
}
