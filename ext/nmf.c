/**
 * NMF: Non-Negative Matrix Factorization
 *
 * Written by Roman Shterenzon
 * (Slightly modified by Y.Tsunesada: just added "const" qualifiers etc.)
 */

#include <math.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h> /* for multiplication */

#define THRESH 0.000001
#define MAXITER 1000

#undef DEBUG

#define mm(a, b) gsl_matrix_mult(a, b)
//gsl_matrix * gsl_matrix_mult(gsl_matrix *a, gsl_matrix *b)
gsl_matrix * gsl_matrix_mult(const gsl_matrix *a, const gsl_matrix *b)
{
  gsl_matrix *c;
  c = gsl_matrix_alloc(a->size1, b->size2);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, a, b, 0.0, c);
  return c;
}

// pretty print
//void pp(gsl_matrix *m)
void pp(const gsl_matrix *m)
{
  int r, c;

  for(r=0; r < (int) m->size1; r++) {
    for(c=0; c < (int) m->size2; c++) {
      printf(" %.2f", gsl_matrix_get(m, r, c));
    }
    printf("\n");
  }
}

/* Returns a distance cost */
//double difcost(gsl_matrix *a, gsl_matrix *b)
double difcost(const gsl_matrix *a, const gsl_matrix *b)
{
  int i, j;
  double dif=0, d;

  for (i=0; i < (int) a->size1; i++)
  {
    for (j=0; j < (int) a->size2; j++)
    {
      d = gsl_matrix_get(a, i, j) - gsl_matrix_get(b, i, j);
      dif += d*d;
	//      dif += pow(gsl_matrix_get(a, i, j) - gsl_matrix_get(b, i, j), 2);
    }
  }
  return dif;
}

static void initmatrix(gsl_matrix *m, double min, double max)
{
  int i,j;
  double val;

  srand(time(NULL));

  for(i=0; i < (int) m->size1; i++)
  {
    for(j=0; j < (int) m->size2; j++)
    {
      val = min + (int) (max * (rand() / (RAND_MAX + min)));
      gsl_matrix_set(m, i, j, val);
    }
  }
}

static double update(gsl_matrix *v, gsl_matrix *w, gsl_matrix *h)
{
  double dist = 0;
  gsl_matrix *wt=NULL, *ht=NULL, *wh=NULL;
  gsl_matrix *w_h=NULL, *wt_w=NULL;
  gsl_matrix *wt_v = NULL;
  gsl_matrix *v_ht=NULL, *wt_w_h=NULL, *w_h_ht=NULL;

  wt = gsl_matrix_alloc(w->size2, w->size1);
  gsl_matrix_transpose_memcpy(wt, w);
  ht = gsl_matrix_alloc(h->size2, h->size1);
  gsl_matrix_transpose_memcpy(ht, h);

  // wt * v
  wt_v = mm(wt, v);

  // wt * w * h
  wt_w = mm(wt, w);
  wt_w_h = mm(wt_w, h);
  gsl_matrix_free(wt_w);

  // h = h.mul_elements(wt * v).div_elements(wt * w * h)
  gsl_matrix_mul_elements(h, wt_v);
  gsl_matrix_div_elements(h, wt_w_h);
  gsl_matrix_free(wt_v);
  gsl_matrix_free(wt_w_h);

  // v * ht
  v_ht = mm(v, ht);

  // w * h * ht
  w_h = mm(w, h);
  w_h_ht = mm(w_h, ht);
  gsl_matrix_free(w_h);

  // w = w.mul_elements(v * ht).div_elements(w * h * ht)
  gsl_matrix_mul_elements(w, v_ht);
  gsl_matrix_div_elements(w, w_h_ht);
  gsl_matrix_free(v_ht);
  gsl_matrix_free(w_h_ht);

  gsl_matrix_free(wt);
  gsl_matrix_free(ht);

  wh = mm(w, h);
  dist = difcost(v, wh);
  gsl_matrix_free(wh);

  // w and h were modified in place
  return dist;
}

/* The main thing - compute the nmf */
int gsl_matrix_nmf(gsl_matrix *v, int cols, gsl_matrix **w, gsl_matrix **h)
{
  double dist = 1;
  int iter = 1;
  double min, max;

#ifdef DEBUG
  printf("\nCols: %d\nv:\n", cols);
  pp(v);
#endif

  gsl_matrix_minmax(v, &min, &max);

#ifdef DEBUG
  printf("Min: %f, Max: %f\n", min, max);
#endif
  *w = gsl_matrix_alloc(v->size1, cols);
  initmatrix(*w, min, max/2); // the multiplicative rules tend to increase w

  *h = gsl_matrix_alloc(cols, v->size2);
  initmatrix(*h, min, max);

  while(dist >= THRESH && iter < MAXITER)
  {
    dist = update(v, *w, *h);
#ifdef DEBUG
    printf("Iteration: %d, distance: %f\n", iter, dist);
    printf("\nw:\n");
    pp(*w);
    printf("\nh:\n");
    pp(*h);
    printf("\n");
#endif
    iter++;
  }
#ifdef DEBUG
  printf("Ended\n");
#endif
  return GSL_SUCCESS;
}
