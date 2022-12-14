/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, lab, kv;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *SOL, *EX_SOL, *X;
  double *AB;
  double *MB;
  
  double temp, relres;

  double opt_alpha;

  /* Size of the problem */
  NRHS=1;
  nbpoints=12;
  la=nbpoints-2;

  /* Dirichlet Boundary conditions */
  T0=5.0;
  T1=20.0;

  printf("--------- Poisson 1D (iter) ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  SOL=(double *) calloc(la, sizeof(double)); 
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  /* Setup the Poisson 1D problem */
  /* General Band Storage */
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS_iter.dat");
  write_vec(EX_SOL, &la, "EX_SOL_iter.dat");
  write_vec(X, &la, "X_grid_iter.dat");

  kv=0;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;
  
  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "EX_AB_iter.dat");
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "MY_AB_iter.dat");

  /********************************************/
  /* Solution (Richardson with optimal alpha) */

  /* Computation of optimum alpha */
  //opt_alpha = richardson_alpha_opt(&la);
  opt_alpha = 0.5;
  printf("Optimal alpha for simple Richardson iteration is : %lf\n",opt_alpha); 

  /* Solve */
  double tol=1e-3;
  int maxit=1000;
  double *resvec;
  int nbite=0;

  resvec=(double *) calloc(maxit, sizeof(double));

  /* Solve with Richardson alpha */
  richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
  printf("Nb iter for richardson alpha: %d\n", nbite);

  /* Richardson General Tridiag */

  /* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
  SOL=(double *) calloc(la, sizeof(double));
  kv = 1;
  ku = 1;
  kl = 1;
  MB = (double *) malloc(sizeof(double)*(lab)*la);
  extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
  write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "J_MB_iter.dat");

  /* Solve with General Richardson */
  richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
  
  write_vec(SOL, &la, "J_MB_SOL_iter.dat");
  write_vec(resvec, &nbite, "J_RESVEC_iter.dat");

  relres = make_relres(EX_SOL,SOL, &la);
  printf("\nNb iter for Jacobi : %d\n", nbite);
  printf("The relative forward error for Jacobi is relres = %e\n",relres);

  MB = (double *) malloc(sizeof(double)*(lab)*la);
  extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
  write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "GS_MB_iter.dat");

  /* Solve with General Richardson */
  SOL=(double *) calloc(la, sizeof(double));
  richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
  
  /* Write solution */
  write_vec(SOL, &la, "GS_MB_SOL_iter.dat");
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  relres = make_relres(EX_SOL,SOL, &la);

  printf("\nNb iter for Gauss-Seidel : %d\n", nbite);
  printf("The relative forward error for Gauss-Seidel is relres = %e\n",relres);

  /* Write convergence history */
  write_vec(resvec, &nbite, "GS_RESVEC_iter.dat");
  
  free(resvec);
  free(RHS);
  free(SOL);
  free(EX_SOL);
  free(X);
  free(AB);
  free(MB);
  printf("\n\n--------- End of iterative methods -----------\n");
}
