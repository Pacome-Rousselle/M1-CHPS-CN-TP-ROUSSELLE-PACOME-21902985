#include "lib_poisson1D.h"

int main(int argc, char **argv)
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
  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  printf("--------- Poisson 1D (iterative) ---------\n\n");
  EX_RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  MY_RHS=(double *) malloc(sizeof(double)*la);
  AB = (double *) malloc(sizeof(double)*lab*la);

  //Setting functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
  write_vec(EX_RHS, &la, "EX_RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  //Richardson

  //Jacobi
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  //Gauss-Seidel
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  
  free(EX_RHS);
  free(MY_RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End of iterative methods -----------\n");    
}