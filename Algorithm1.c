#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){
  // Initial matrices
  double *A, *Ad, *Adu, *Adl,  *B, *L, *R;
  // Initial vectors
  double *x, *y;
  // Initial parameters
  double l,r;

  // Computed Matrices
  double *BtB;
  // Computer vectors
  double *b;


  /* creation of parameters for BLAS routines */
  char lower = 'L';
  char upper = 'U';
  char no = 'N';
  char yes ='Y';
  char trans = 'T';
  char left ='L';
  char right = 'R';
  int zero = 0;
  int one = 1;
  int two = 2;
  int info;

  int i; // iterator
  int n; // for the dimension/size problem
  int nm1; //n minus 1




  // Matrix A is tridiagonal and stored using 3 arrays
  // A = (double *) malloc( n*n * sizeof(double) );

  // TODO: Change the allocation size to store B in a Packed Storage. This is valid because B is symmetric
  B = (double *) malloc( n*n * sizeof(double) );
  L = (double *) malloc( n*n * sizeof(double) );
  R = (double *) malloc( n*(n-1) * sizeof(double) );

  // Allocate 3 arrays for the tridiagonal matrix A
  Ad = (double *) malloc( n*1 * sizeof(double) );       //A diagonal
  Adu = (double *) malloc( (n-1)*1 * sizeof(double) );  // A upper diagonal
  Adl = (double *) malloc( (n-1)*1 * sizeof(double) );  // A lower diagonal


  // Matrices to store computed results
  BtB = (double *) malloc( n*n * sizeof(double) );
  RA = (double *) malloc( (n-1)*n * sizeof(double) );
  P = (double *) malloc( n*(n-1) * sizeof(double) );
  // Vectors to store computed results
  b = (double *) malloc( n * sizeof(double) );
  
  // Seed for randomizer
  srand( time(NULL) );
  double random = (double)rand()/(double)RAND_MAX;

    printf("random =%f \n", random);


  // Initialize y
  for(i=0;i<n;i++){
    y[i] = (double)rand()/(double)RAND_MAX;
  }

  // Initialize Ad,Adu,Adl
  for(i=0;i<n-1;i++){
    Ad[i] = (double)rand()/(double)RAND_MAX;
    Adu[i] = (double)rand()/(double)RAND_MAX;
    Adl[i] = Adu[i];
  }
    Ad[n] =(double)rand()/(double)RAND_MAX;

  // TODO : Create copy of A(Ad,Adu,Adl)
  

  // 1.WORKING  solve A*y'=y
  dgtsv(&n, &one, Adl, Ad, Adu, y, &n, &info);
  // 2. compute BtB=B*B (B*B = B*transpose(B))
  dsymm(&left, &upper, &n, &n, &one, B, &n, B, &n, &zero, BtB, &n);
  // 3.
  dsbmv(&upper, &n, &two, &one, BtB, &n, y, &one, &zero, b, &one);
  // 4. compute RA=R*A
  dsymm(&right, &upper, &nm1, &n, &one, A, &n, R, &n, &zero, RA, &n);

  // TODO: for multiple iterations, insert loop beginning here

  // 5. compute v=R*x
  dgemv(&trans, &nm1, &n, &one, R, &n, x, &one, &zero, v, &one);

  // TODO: compute L

  // 6. compute P=transpose(RA)*L
  dgemm(&trans, &no, &n, &nm1, &nm1, &one, RA, &n, L, &nm1, &zero, P, &one);

  // TODO: Set Q = copy(BtB). Use Band Storage!!!
  // 7. compute Q = P*RA+Q
  dgemm(&no, &no, &n, &nm1, &n, &one, P, &n, RA, &nm1, &one, Q, &one);
  //8. solve Q*b'=b. If psbv leads to errors, use gbsv
  //dgbsv(&n, &two, &two, &one, QBandStorage, &n, ??ipiv?, b, &n, &info );
  dpbsv(&upper, &n, &two, &one, QBandStorage, &n, b, &n, &info );

  // TODO: set x=copy(b)?? to check

  // compute x=A*b
  dgbmv(&no, &n, &n, &one, &one, &one, A, &n, b, &one, &zero, x, &one);

  // Freeing all the allocated memory
  free(B);
  free(L);
  free(R);
  free(y);
  free(Ad);
  free(Adu);
  free(Adl);
  free(BtB);
  free(RA);
  free(P);
  free(b);

}
