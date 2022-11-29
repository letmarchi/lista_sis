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

double *Residuos(double **M, int m, int n, double *x){
  int i;
  double *r,*s;

  r=malloc(m*sizeof(double *));
  s=MultiVetMat(M, x, m, n);
  for(i=0;i<m;i++)
    r[i]=M[i][n-1]-s[i];

  return r;
}

int main() {
double **M, **C,*b ,*v,*d,*l1, dx, tolerance=1e-8;
int m, n, l, i, it=0, p=0;
FILE *arq;

M = LeMatriz("Matrix.dat",&m, &n);
v = LeVetor("vetor.dat", &l);

//teste

return 0;
}