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
  int ld;
  for(int i = 0; i<(*la); i++)
  {
    // ld*j + i => (i,j)
    ld = i*(*lab);
    if((*kv) > 0)
    {
      AB[ld + *kv] = 0.0;
    }
    AB[ld+*kv] = -1.0;
    AB[ld+*kv+1] = 2.0;
    AB[ld+*kv+2] = -1.0;
  }
  AB[0]=0.0;
  if((*kv) > 0){AB[1] = 0.0;}
  AB[(*lab)*(*la)-1] = 0.0;
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

double make_relres(double *analytic, double *experimental, int *la)
{
  // norm of x 
  double normx = cblas_dnrm2(*la,analytic,1);

  // (x = -^x + x) 
  cblas_daxpy(*la,-1.0,experimental,1,analytic,1);

  // Relative forward error : ||(x = -^x + x)||/||x||
  return (cblas_dnrm2(*la,analytic,1)/normx);
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
  double lambda = 1/((*la)+1); // h
  lambda *= M_PI_2; // h*pi/2
  lambda *= (*la); // n*h*pi/2
  lambda = sin(lambda); // sin(n*h*pi/2)
  lambda *= lambda; // sin²(n*h*pi/2)
  lambda *= 4; // 4sin²(n*h*pi/2)
  return lambda;
}

double eigmin_poisson1D(int *la){
  double lambda = 1/((*la)+1); // h
  lambda *= M_PI_2; // h*pi/2
  lambda = sin(lambda); // sin(h*pi/2)
  lambda *= lambda; // sin²(h*pi/2)
  lambda *= 4; // 4sin²(h*pi/2)
  return lambda;
}

double richardson_alpha_opt(int *la){
  double opt = eigmax_poisson1D(la) + eigmin_poisson1D(la);
  opt = 2/opt;
  return opt;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  double *tmp = malloc(sizeof(double)*(*la));
  cblas_dcopy(*la,RHS,1,tmp,1);

  //  Initialize r^0
  //  r^0 := b - Ax^0
  cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,tmp,1);

  //  ||tmp||
  double check = cblas_dnrm2(*la,tmp,1);
  *nbite=0;
  resvec[*nbite] = check;

  while(check > (*tol) && *nbite < *maxit)
  {
    cblas_dcopy(*la,RHS,1,tmp,1);
    //  tmp^k := b - Ax^k
    cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,tmp,1);

    //  x^(k+1) := x^k + α*tmp^k
    cblas_daxpy(*la,(*alpha_rich),tmp,1,X,1);

    //  ||tmp||
    check = cblas_dnrm2(*la,tmp,1);
    (*nbite)++;
    resvec[*nbite] = check;
  }
  free(tmp);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  int ld;
  for(int i = 0; i<(*la); i++)
  {
    // ld*j + i => (i,j)
    ld = i*(*lab);
    MB[ld+*kv] = AB[ld+*kv];
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  int ld;
  for(int i = 0; i<(*la); i++)
  {
    // ld*j + i => (i,j)
    ld = i*(*lab);
    MB[ld+*kv] = AB[ld+*kv];
    MB[ld+*kv+1] = AB[ld+*kv+1];
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite){
  //r
  double *tmp = malloc(sizeof(double)*(*la));
  cblas_dcopy(*la,RHS,1,tmp,1);

  //  Initialize r^0
  //  r^0 := b - Ax^0
  cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,tmp,1);

  int info; int NHRS = 1;
  int *ipiv = (int *)calloc(*la, sizeof(int));
  int ku_mb = (*ku)-1;
  dgbtrf_(la,la,kl,&ku_mb,MB,lab,ipiv,&info);

  if(info != 0)
  {
    printf("LU factorisation error : info = %d\n", info);
  }

  //  ||tmp||
  double check = cblas_dnrm2(*la,tmp,1);
  *nbite=0;
  resvec[*nbite] = check;
  
  while((check > (*tol)) && ((*nbite) < (*maxit)))
  {
    cblas_dcopy(*la,RHS,1,tmp,1);
    //  tmp^k := b - Ax^k
    cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,tmp,1);

    //  tmp^k := D^-1*tmp^k (= tmp^k/D)
    dgbtrs_("N",la,kl,&ku_mb,&NHRS,MB,lab,ipiv,tmp,la,&info,*la);

    //  x^(k+1) := x^k + tmp^k
    cblas_daxpy(*la,1.0,tmp,1,X,1);

    //  ||tmp||
    check = cblas_dnrm2(*la,tmp,1);
    (*nbite)++;
    resvec[*nbite] = check;
  }
  free(ipiv);
  free(tmp);
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


