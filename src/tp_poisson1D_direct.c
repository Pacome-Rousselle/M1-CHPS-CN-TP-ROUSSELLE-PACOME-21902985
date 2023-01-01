/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *TEST, *EX_RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;

  double temp, relres;

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  EX_RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  TEST=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(EX_RHS, &la, "EX_RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
  
  printf("DGBMV\n");
  // TEST <- AB*EX_SOL
  cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,kl,ku,1.0,AB+1,lab,EX_SOL,1,0.0,TEST,1);
  write_vec(TEST, &la, "TEST.dat");

  /* Relative forward error */
  relres = make_relres(EX_RHS,TEST, relres);
  printf("\nThe relative forward error for dgbmv is relres = %e\n",relres);

  printf("\nSolution with LAPACK\n");

  /* LU Factorization */
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(TEST,&la,&T0,&T1);

  info=0;
  ipiv = (int *) calloc(la, sizeof(int));
  dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "EX_LU.dat");

  /* Solution (Triangular) */
  if (info==0){
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, TEST, &la, &info, la);
    write_vec(TEST, &la, "SOL_LU.dat");
    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  }else{
    printf("\n INFO = %d\n",info);
  }

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(TEST,&la,&T0,&T1);

  ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "MY_LU.dat");
  /* Solution (Triangular) */

  if (info==0){
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, TEST, &la, &info, la);
    write_vec(TEST, &la, "MY_SOL_LU.dat");
    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  }else{
    printf("\n INFO = %d\n",info);
  }

  /* Relative forward error */
  relres = make_relres(EX_SOL,TEST, relres);
  printf("\nThe relative forward error is relres = %e\n",relres);

  /* It can also be solved with dgbsv */
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);
  dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, EX_RHS, &la, &info);
  write_xy(EX_RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  // relres = make_relres(EX_RHS,TEST, relres);
  // printf("\nThe relative forward error is relres = %e\n",relres);

  free(EX_RHS);
  free(TEST);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
