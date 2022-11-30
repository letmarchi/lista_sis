#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double **LeMatriz(char nome[20], int *m, int *n){
  double **a;
  int i, j;
  FILE *arq;
  
  arq=fopen(nome, "r");
  fscanf(arq, "%d", m);
  fscanf(arq, "%d", n);
  a = (double **) malloc(*m*sizeof(double *));
  for (i=0;i<*m;i++) a[i]= (double *) malloc((*n)*sizeof(double));
  for (i=0;i<*m;i++){
    for (j=0;j<*n;j++){
      fscanf(arq, "%lf", &a[i][j]);
    }
  }
  return a;
}

void **ImprimeMatriz(double **M, int *m, int *n){
  int i, j;
  
  for (i=0;i<*m;i++){
    for (j=0;j<*n;j++){
      printf("%lf ", M[i][j]);
    }
  puts("");
  }
  
}

double *LeVetor(char *nome, int *m){
  FILE *fp = fopen(nome, "r");
  int i, j;
  double *vetor;  

  fscanf(fp, "%d", m);
   
  
  vetor = (double *)malloc( *m*sizeof(double *));
  
  for(i=0;i<*m;i++){
          fscanf(fp, "%lf", &vetor[i]);
  }

  return vetor;
}


void ImprimeVetor(double *v, int j){
  int i;
  for(i=0; i<j; i++) printf("v[%d] = %g\n", i, v[i]);
}

double NormaMatriz(double **v, int m, int n, int p){
  double max=0, aux;
  int i,j;

  if(p==0){
    for(i=0;i<m;i++){
      aux=0;
      for(j=0;j<n;j++){
        aux+=fabs(v[i][j]);
      }
      if(aux>max) max=aux;
    }
  }
  else if(p==1){
    for(i=0;i<n;i++){
      aux=0;
      for(j=0;j<m;j++){
        aux+=fabs(v[j][i]);
      }
      if(aux>max) max=aux;
    }
  }
  else{
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
          max+=pow(v[i][j],2);
        }
    }
    max=sqrt(max);
  }
  return max;
}

double NormaVetor(double *v, int m, int p){
  
  int i;
  double max=0;

  if(p==0){
    for(i=0;i<m;i++) {
      if(fabs(v[i])>max) max=fabs(v[i]);
    }
  }
  else{
    for(i=0;i<m;i++) {
      max+=pow(fabs(v[i]),p); 
     }
    max=pow(max,1.0/p); 
  }
  return max;
}

double *MultiVetMat(double **M, double *x, int m, int n){
  int i, j, k;
  double *r;

  r=malloc(m*sizeof(double *));

  for(i=0; i<m; i++){
    for(j=0; j<m; j++){
      r[i]+=M[i][j]*x[j];
    }
  }
  return r;
}

double MultVetor(double *x, double *y, int t){
  int i;
  double l=0;
  for(i=0;i<t;i++){
    l+=x[i]*y[i];
  }
  return l;
}

void TrocaLinha(double **k, double **l){
  double *temp; 

  temp = *k;
  *k = *l;
  *l = temp;
}

void **LUPivot(double **M, double ***L, double ***U, double **b, int m, int n ){
  int i, j, k, p;
  double lbd, pivot;

  *L = calloc(m, sizeof(double *));
  for(i=0; i<m; i++) (*L)[i] = calloc(m, sizeof(double));

  *U = malloc(m*sizeof(double *));
  for(i=0; i<m; i++) (*U)[i] = malloc(m*sizeof(double));

  *b = malloc(m*sizeof(double));
  for(i=0; i<m; i++) (*b)[i] = M[i][n-1];

  for(i=0; i<m; i++){
    for(j=0; j<m; j++){
      (*U)[i][j] = M[i][j];
    }
  }

  for(j=0; j<m-1; j++){
    pivot = (*U)[j][j];
    p=0;

      for(i=j+1; i<m; i++){
        if(fabs((*U)[i][j])>pivot){
          pivot=fabs((*U)[i][j]);
          p=i;
        }
      }

      if(p){
        TrocaLinha(&(*U)[p], &(*U)[j]);
        TrocaLinha(&(*L)[p], &(*L)[j]);
        lbd=(*b)[p];
        (*b)[p]=(*b)[j];
        (*b)[j]=lbd;
      }    

    for(i=j+1; i<m; i++){
      lbd=(*L)[i][j] = (*U)[i][j]/(*U)[j][j];
      for(k=j; k<m; k++) (*U)[i][k] -= lbd*(*U)[j][k];
    }
  }
  for(i=0; i<m; i++){
    (*L)[i][i] = 1;
  }
}

double *SubstituicaoDireta(double **M, int m, double *b){
  double *v, sum; 
  int i,j;

  v = malloc(m* sizeof(double));

  for(i=0; i<m; i++){
    sum = 0;
    for(j=i-1; j>=0; j--){
      sum+= M[i][j]*v[j];
    }
    v[i] = (b[i] - sum)/M[i][i];
  }
  return  v;
}

double *SubstituicaoReversa(double **M, int m, double *b){
  double *v,sum;
  int i,j;

  v = malloc(m*sizeof(double));

  for(i=m-1; i>=0; i--){
    sum=0;
    for(j=i+1; j<= m-1; j++) {
      sum+= M[i][j]*v[j];
    }
    v[i] = (b[i]-sum)/M[i][i];
  }
  return v;
}

double f1(double *x){
  return 4*(x[0])-(x[1])+(x[2])-(x[0]*x[3]);
}

double f2(double *x){
  return -(x[0])+3*(x[1])-2*(x[2])-(x[1]*x[3]);
}

double f3(double *x){
  return x[0]-2*x[1]+3*x[2]-x[2]*x[3];
}

double f4(double *x){
  return pow(x[0],2)+pow(x[1],2)+pow(x[2],2)-1;
}

void Jacobiana(double **J, double *x) {
  J[0][0]=4-x[3];
  J[0][1]=-1;
  J[0][2]=1;
  J[0][3]=-x[0];
  J[1][0]=-1;
  J[1][1]=3-x[3];
  J[1][2]=-2;
  J[1][3]=-x[1];
  J[2][0]=1;
  J[2][1]=-2;
  J[2][2]=3-x[3];
  J[2][3]=-x[2];
  J[3][0]=2*x[0];
  J[3][1]=2*x[1];
  J[3][2]=2*x[2];
  J[3][3]=0;
}

int main() {
int i, itera=0, m;
double *b, *x, *x2, *x3, **J, **L, **U, norma, tol=1e-8, (*f[])(double*)={f1,f2,f3,f4}; 
FILE *arq;
  
 //chute 1

  x = LeVetor("chute1.dat", &m);
  printf("Chute 1:");
  ImprimeVetor(x, m);

  J = (double **) malloc(m*sizeof(double *));
  for(i=0; i<m;i++) J[i] = (double *) malloc((m+1)*sizeof(double));

  do{
      Jacobiana(J, x);
      for(i=0; i<m; i++) J[i][m] = -f[i](x);

      LUPivot(J, &L, &U, &b, m, m+1);

      b = SubstituicaoDireta(L, m, b);
      b = SubstituicaoReversa(U, m, b);
      itera++;

      norma = NormaVetor(b, m, 0);
      for(i=0; i<m; i++) x[i] += b[i];
  
  }while ((norma>tol)&&(itera<20));

  puts("");
  printf("solucoes: ");
  ImprimeVetor(x, m);

  //chute 2

  x2 = LeVetor("chute2.dat", &m);
  printf("Chute 2:");
  ImprimeVetor(x2, m);

  free(J);

  do{
      Jacobiana(J, x2);
      for(i=0; i<m; i++) J[i][m] = -f[i](x2);

      LUPivot(J, &L, &U, &b, m, m+1);

      b = SubstituicaoDireta(L, m, b);
      b = SubstituicaoReversa(U, m, b);
      itera++;

      norma = NormaVetor(b, m, 0);
      for(i=0; i<m; i++) x2[i] += b[i];
  
  }while ((norma>tol)&&(itera<20));

  puts("");
  printf("solucoes: ");
  ImprimeVetor(x2, m);

  //chute 3

  x3 = LeVetor("chute3.dat", &m);
  printf("Chute 3:");
  ImprimeVetor(x3, m);

  J = (double **) malloc(m*sizeof(double *));
  for(i=0; i<m;i++) J[i] = (double *) malloc((m+1)*sizeof(double));

  do{
      Jacobiana(J, x3);
      for(i=0; i<m; i++) J[i][m] = -f[i](x3);

      LUPivot(J, &L, &U, &b, m, m+1);

      b = SubstituicaoDireta(L, m, b);
      b = SubstituicaoReversa(U, m, b);
      itera++;

      norma = NormaVetor(b, m, 0);
      for(i=0; i<m; i++) x3[i] += b[i];
  
  }while ((norma>tol)&&(itera<20));

  puts("");
  printf("solucoes: ");
  ImprimeVetor(x3, m);
  
  return 0;
}