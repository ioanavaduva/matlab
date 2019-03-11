#include <stdlib.h>
#include <stdio.h>
#include "hsl_mi20d.h"

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

int main(void) {
   const int n = 10; /* size of system to solve */

   /* matrix data */
   int *ptr, *col;
   double *val;

   /* derived types */
   void *keep;
   struct mi20_control control;
   struct mi20_solve_control solve_control;
   struct mi20_info info;

   /* others */
   int i;
   double rhs[n];
   double sol[n];

   /* generate matrix A */
   matrix_gen(n, &ptr, &col, &val);
   
   /* set right hand side to vector of ones */
   for(i=0; i<n; i++) rhs[i] = 1;
   
   /* initalize controls */
   mi20_default_control(&control);
   mi20_default_solve_control(&solve_control);

   /* call mi20_setup */
   mi20_setup_csr(n, ptr, col, val, &keep, &control, &info);
   if (info.flag < 0) {
     printf("Error return from mi20_setup\n");
     return 1;
   }

   mi20_solve(rhs, sol, &keep, &control, &solve_control, 
	      &info);
   if (info.flag < 0) {
     printf("Error return from mi20_solve\n");
   } else {
     printf(" Convergence in %d iterations\n", info.iterations);
     printf(" 2-norm of residual = %e\n", info.residual);
   } 

   /* deallocation */
   mi20_finalize(&keep, &control, &info);
   free(ptr); free(col); free(val);

   return 0; /* sucess */
}
