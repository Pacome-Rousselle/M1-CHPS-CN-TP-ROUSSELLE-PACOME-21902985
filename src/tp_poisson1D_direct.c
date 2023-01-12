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
  double *MY_RHS, *EX_RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;

  double temp, relres;

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D (direct) ---------\n\n");
  EX_RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  MY_RHS=(double *) malloc(sizeof(double)*la);

  //Setting functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(EX_RHS, &la, "EX_RHS_direct.dat");
  write_vec(EX_SOL, &la, "EX_SOL_direct.dat");
  write_vec(X, &la, "X_grid_direct.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "MY_AB_direct.dat");
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "EX_AB_direct.dat");
  
  printf("DGBMV\n");
  // MY_RHS <- AB*EX_SOL
  cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,kl,ku,1.0,AB+1,lab,EX_SOL,1,0.0,MY_RHS,1);
  write_vec(MY_RHS, &la, "MY_RHS_direct.dat");

  /* Validation for dgbmv : relative forward error */
  relres = make_relres(EX_RHS,MY_RHS, &la);
  printf("\nThe relative forward error for dgbmv is relres = %e\n",relres);

  printf("\nSolution with LAPACK\n");

  /* LU Factorization */
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(MY_RHS,&la,&T0,&T1);

  info=0;
  ipiv = (int *) calloc(la, sizeof(int));
  dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "EX_LU_direct.dat");

  /* Solution (Triangular) */
  if (info==0){
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, MY_RHS, &la, &info, la);
    write_vec(MY_RHS, &la, "EX_SOL_LU_direct.dat");
    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  }else{
    printf("\n INFO = %d\n",info);
  }

  /* Validation of LU for tridiagonal matrix (dgbtrf/dgbtrs) */

  relres = make_relres(EX_SOL,MY_RHS, &la);
  printf("\nThe relative forward error for dgbtrf/dgbtrs is relres = %e\n",relres);

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(MY_RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

  ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "MY_LU_direct.dat");

  /* Solution (Triangular) */

  if (info==0){
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, MY_RHS, &la, &info, la);
    write_vec(MY_RHS, &la, "MY_SOL_LU_direct.dat");
    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  }else{
    printf("\n INFO = %d\n",info);
  }

  /* Validation of our LU method for tridiagonal matrix */
  relres = make_relres(EX_SOL,MY_RHS, &la);
  printf("\nThe relative forward error for dgbtrftridiag is relres = %e\n",relres);

  /* It can also be solved with dgbsv */
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

  dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, EX_RHS, &la, &info);
  write_xy(EX_RHS, X, &la, "SOL_direct.dat");

  /* Relative forward error for dgbsv*/
  relres = make_relres(EX_SOL,EX_RHS,&la);
  printf("\nThe relative forward error for dgbsv is relres = %e\n",relres);

  /* Complexity testing */
  clock_t top;
  
  FILE *direct = fopen("direct_graphd.dat", "w");
  if(direct == NULL){perror("direct_graphd.dat");}

  FILE *dgbtrf = fopen("dgbtrf_graphd.dat", "w");
  if(dgbtrf == NULL){perror("dgbtrf_graphd.dat");}

  FILE *errav_direct = fopen("ferror_graphd.dat", "w");
  if(errav_direct == NULL){perror("ferror_graphd.dat");}

  for(nbpoints = 100; nbpoints < 10000; nbpoints += 100)
  {
    // Realloc
    la=nbpoints-2;
    EX_RHS=(double *) realloc(EX_RHS, sizeof(double)*la);
    EX_SOL=(double *) realloc(EX_SOL, sizeof(double)*la);
    AB = (double *) realloc(AB, sizeof(double)*lab*la);
    X = (double *) realloc(X, sizeof(double)*la);
    set_grid_points_1D(X, &la);
    free(ipiv); ipiv = (int *) calloc(la, sizeof(int)); // realloc at 0

    // dgbtrftridiag
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);

    top = clock();
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    fprintf(dgbtrf, "%f ",((double)(clock() - top)/CLOCKS_PER_SEC));

    // dgbtrf
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);

    top = clock();
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    fprintf(dgbtrf, "%f\n",((double)(clock() - top)/CLOCKS_PER_SEC));

    // dgbtrftridiag + dgbtrs
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);

    top = clock();
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, EX_RHS, &la, &info, la);    
    fprintf(direct, "%f ",((double)(clock() - top)/CLOCKS_PER_SEC));

    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
    relres = make_relres(EX_SOL,EX_RHS,&la);
    fprintf(errav_direct, "%e ", relres);

    // dgbtrf + dgbtrs
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);

    top = clock();
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, EX_RHS, &la, &info, la);    
    fprintf(direct, "%f ",((double)(clock() - top)/CLOCKS_PER_SEC));

    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
    relres = make_relres(EX_SOL,EX_RHS,&la);
    fprintf(errav_direct, "%e ", relres);

    // dgbsv
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);

    top = clock();
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, EX_RHS, &la, &info);
    fprintf(direct, "%f\n",((double)(clock() - top)/CLOCKS_PER_SEC));

    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
    relres = make_relres(EX_SOL,EX_RHS,&la);
    fprintf(errav_direct, "%e\n", relres);
  }

  free(EX_RHS);
  free(MY_RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End of direct methods -----------\n");
}
