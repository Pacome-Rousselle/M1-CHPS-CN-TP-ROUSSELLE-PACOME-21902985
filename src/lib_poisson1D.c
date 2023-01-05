/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

/**********************************************/
/*          Direct methods section            */
/**********************************************/
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  // lab = nb colonnes
  // la = nb lignes
  
  for(int i = 0; i<(*la); i++)
  {
    // ld*j + i => (i,j)
    AB[i*(*lab)] = 0;
    AB[i*(*lab)+1] = -1;
    AB[i*(*lab)+2] = 2;
    AB[i*(*lab)+3] = -1;
  }
  AB[1] = 0; AB[(*lab)*(*la)-1] = 0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // Left condition T(0) = T0
  int i = 0;
  RHS[i] = *BC0;
  // Vector body between left and right conditions set at 0
  for(i = 1; i < (*la)-1; i++)
  {
    RHS[i] = 0;
  }
  // Right condition T(n-1) = T1
  RHS[i] = *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // Invariant extraction
  double diff = *BC1 - *BC0;
  for(int i = 0; i < *la; i++)
  {
    // T(x) = T0 + x(T1-T0)
    EX_SOL[i] = *BC0 + X[i]*diff;
  }
}  

void set_grid_points_1D(double* x, int* la){
  // X is equally segmented with steps of 1/la
  // X(i-1) = X(i) - 1/la
  // X(i+1) = X(i) + 1/la
  for(int i = 0; i < *la; i++)
  {
    x[i] = (i+1)/ (1.0*(*la)+1);
  }
}

double make_relres(double *analytic, double *experimental, double relres)
{
  double *cpy = malloc(sizeof(analytic));
  cblas_dcopy(1,analytic,1,cpy,1); //cpy is a copy of analytic
  // norm of x is the square root of the sum of cpy's arguments squared 
  // (||x|| = ddot(cpy))
  double normx = cblas_ddot(1,cpy,1,cpy,1); 
  sqrt(normx);

  // (x = -^x + x) 
  cblas_daxpy(1,-1.0,experimental,1,analytic,1);

  // Relative forward error :
  relres = (sqrt(cblas_ddot(1,analytic,1,analytic,1)))/normx;
  free(cpy);
  return relres;
}

int indexABCol(int i, int j, int *lab){
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info)
{
  double pivot = 0;

  for (int i = 1; i < (*la); i++)
  {
    pivot = AB[i*(*lab)+1]*AB[(i-1)*(*lab)+3]; // pivot = b(n-1)*c(n-1)
    pivot /= AB[(i-1)*(*lab)+2]; // pivot = b(n-1)c(n-1)/a(n-1)
    pivot *= -1.0; // pivot = -(b(n-1)c(n-1)/a(n-1))

    AB[i*(*lab)+2] += pivot; // dn = a(n-1) - b(n-1)*c(n-1)/a(n-1)
    AB[(i-1)*(*lab)+3] /= AB[(i-1)*(*lab)+2]; 
  }
  return *info;
}

/**********************************************/
/*         Iterative methods section          */
/**********************************************/
void eig_poisson1D(double* eigval, int *la){
}

double eigmax_poisson1D(int *la){
  return 0;
}

double eigmin_poisson1D(int *la){
  return 0;
}

double richardson_alpha_opt(int *la){
  return 0;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){

}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){

}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){

}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){

}
/**********************************************/
/*           File writing section             */
/**********************************************/
void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  


