#include <stdlib.h>
#include <stdio.h>
#include "hsl_mi20d.h"

/* MI21 Fortran routines (no C interface available) */
/* As these are F77 style codes, we assume that the C binding merely appends
   an underscore to the name of the Fortran routine, and that all data types
   match and are pass by reference. THIS WILL NOT WORK ON ALL COMPILERS. */
void mi21id_(int *icntl, double *cntl, int *isave, double *rsave);
void mi21ad_(int *iact, const int *n, double *w, const int *ldw, int *locy,
      int *locz, double *resid, int *icntl, double *cntl, int *info, int *isave,
      double *rsave);

/* Generate example matrix */
void matrix_gen(const int n, int **ptr, int **col, double **val) {
   int i,nnz,p;

   nnz = n + 2*(n-1);
   *ptr = (int*) malloc((n+1)*sizeof(int));
   *col = (int*) malloc(nnz*sizeof(int));
   *val = (double*) malloc(nnz*sizeof(double));
   p = 0; /* pointer to next empty position */
   for(i=0; i<n; i++) {
     (*ptr)[i] = p;
     if (i==0) {   /* first row */
       (*col)[p] = i;    (*col)[p+1] = i+1;
       (*val)[p] = 2.0;  (*val)[p+1] = -1.0;
       p = p+2;
     } else if (i==n-1) { /* last row */
       (*col)[p] = i-1;  (*col)[p+1] = i;
       (*val)[p] = -1.0; (*val)[p+1] = 2.0;
       p = p+2;
     } else {
       (*col)[p] = i-1;  (*col)[p+1] = i;   (*col)[p+2] = i+1;
       (*val)[p] = -1.0; (*val)[p+1] = 2.0; (*val)[p+2] = -1.0;
       p = p+3;
     }
   }
   (*ptr)[n] = nnz;
}

/* Calculate b = Ax */
void spmv(int n, const int ptr[], const int col[], const double val[],
      double b[], const double x[]) {
   int i, j;

   for(i=0; i<n; i++) {
      b[i] = 0;
      for(j=ptr[i]; j<ptr[i+1]; j++) {
         b[i] += val[j] * x[col[j]];
      }
   }
}

int main(void) {
   const int n = 10; /* size of system to solve */

   /* matrix data */
   int *ptr, *col;
   double *val;

   /* derived types */
   void *keep;
   struct mi20_control control;
   struct mi20_info info;

   /* Arrays and scalars required by the CG code mi21 */
   double cntl[5], rsave[6];
   int icntl[8],isave[10],info21[4];
   double w[n*4];
   double resid;
   int locy, locz, iact, i;
   
   /* generate matrix A */
   matrix_gen(n, &ptr, &col, &val);
   
   /* Prepare to use the CG code mi21 with preconditioning */
   mi21id_(icntl, cntl, isave, rsave);
   icntl[3-1] = 1;

   /* set right hand side to vector of ones */
   for(i=0; i<n; i++) w[i] = 1;

   /* initalize control */
   mi20_default_control(&control);
   
   /* call mi20_setup */
   mi20_setup_csr(n, ptr, col, val, &keep, &control, &info);
   if (info.flag < 0) {
     printf("Error return from mi20_setup\n");
     return 1;
   }

   /* solver loop */
   iact = 0;
   while(1) {
      mi21ad_(&iact, &n, w, &n, &locy, &locz, &resid, icntl, cntl, info21,
         isave, rsave);
      if (iact == -1) {
         printf("Error in solver loop\n");
         break;
      } else if (iact == 1) {
         printf(" Convergence in %d iterations\n", info21[2-1]);
         printf(" 2-norm of residual = %e\n", resid);
         break;
      } else if (iact == 2) {
         spmv(n, ptr, col, val, &w[n*(locy-1)], &w[n*(locz-1)]);
      } else if (iact == 3) {
         mi20_precondition(&w[n*(locz-1)], &w[n*(locy-1)], &keep, &control, 
            &info);
         if (info.flag < 0) {
            printf("Error return from mi20_precondition\n");
            break;
         }
      }
   }

   /* deallocation */
   mi20_finalize(&keep, &control, &info);
   free(ptr); free(col); free(val);

   return 0; /* success */
}
