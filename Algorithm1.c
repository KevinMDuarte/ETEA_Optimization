#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"
#include <string.h>

int main(int argc, char *argv[]){
 
  double s_initial,s_elapsed;	

  // Initial matrices
  double *A, *Ad, *Adu, *Adl,  *B, *L, *R;
  // Initial vectors
  double *x, *y, *ycpy;
  // Initial parameters
  double l,r;
  r = 0.45;
  l = 0.75;	 


  // Computed Matrices
  double *BtB, *RA, *P, *Q;
  // Computer vectors
  double *b, *v;


  /* creation of parameters for BLAS routines */
  char lower = 'L';
  char upper = 'U';
  char no = 'N';
  char yes ='Y';
  char trans = 'T';
  char left ='L';
  char right = 'R';
  int iZero = 0;
  int iOne = 1;
  double dZero = 0.;
  double dOne = 1.;
  int two = 2;
  int three = 3;
  int info;
  int *ipiv;

  int n = 10;
  int nm1 = n-1; //n minus 1

  // Iterators
  int i,j;


  // TODO: Change the allocation size to store B in a Packed Storage. This is valid because B is symmetric
  B = (double *) malloc( n*n * sizeof(double) );

  A = (double *) malloc( n*n * sizeof(double) );
  L = (double *) malloc( nm1*nm1 * sizeof(double) );
  R = (double *) malloc( n*nm1 * sizeof(double) );
  y = (double *) malloc( 1*n * sizeof(double) );
  ycpy = (double *) malloc( 1*n * sizeof(double) );
  x = (double *) malloc( 1*n * sizeof(double) );

  // Allocate 3 arrays for the tridiagonal matrix A
  Ad = (double *) malloc( 1*n * sizeof(double) );       //A diagonal
  Adu = (double *) malloc( 1*nm1 * sizeof(double) );  // A upper diagonal
  Adl = (double *) malloc( 1*nm1 * sizeof(double) );  // A lower diagonal


  //Matrices to store computed results
  BtB = (double *) malloc( n*n * sizeof(double) );
  RA = (double *) malloc( n*nm1 * sizeof(double) );
  P = (double *) malloc( nm1*n * sizeof(double) );
  Q = (double *) malloc( n*n * sizeof(double) );
  //Vectors to store computed results
  b = (double *) malloc( 1*n * sizeof(double) );
  v = (double *) malloc( 1*nm1 * sizeof(double) );
  ipiv = (int *) malloc( 1*n*sizeof(int) );


  srand( time(NULL) );
  double random = (double)rand()/(double)RAND_MAX;

	//printf("random =%f \n", random);

// Initialize ipiv
for(i=0;i<n;i++){
	ipiv[i] = 0;
}

// Initialize and set values of y
for(i=0;i<n;i++){
	//y[i] = (double)rand()/(double)RAND_MAX;
	y[i] = 0.5;
}

// Copy y into ycpy
memcpy(ycpy,y,n*sizeof(double));

// Set values of Ad,Adu,Adl
for(j=0;j<n-1;j++){
	Ad[j] = 2.;//(double)rand()/(double)RAND_MAX;
	Adu[j] = 3.;//(double)rand()/(double)RAND_MAX;
	Adl[j] = Adu[j];
}
Ad[n-1]= 2.;//(double)rand()/(double)RAND_MAX;

// Initialize B with 0
for(j=0;j<n;j++){
	for(i=0;i<n;i++){
		B[j*n+i] = 0;
	}
}

// Set values of B
B[0] = 4. ;
for(j=1;j<n;j++){
	B[j*n+j] = 4.;
	B[j*n+j-1] = 5.;
	B[(j-1)*n+j] = B[j*n+j-1];
}

// Initialize BtB with 0
for(j=0;j<n;j++){
	for(i=0;i<n;i++){
		BtB[j*n+i] = 0;
	}
}

// Initialize b with 0
for(j=0;j<n;j++){
	b[j]=0;
}

// Initialize A with 0
for(j=0;j<n;j++){
	for(i=0;i<n;i++){
		A[j*n+i] = 0;
	}
}

// Set values of A using Ad,Adu and Adl
A[0] = Ad[0];
for(j=1;j<n;j++){
	A[j*n+j] = Ad[j];
	A[j*n+j-1] = Adl[j-1];
	A[(j-1)*n+j] = Adu[j-1];
}

// Initialize R with 0
for(j=0;j<n;j++){
	for(i=0;i<nm1;i++){
		R[j*nm1+i] = 0;
	}
}
  
// Set values to R using parameter r
for(j=0;j<nm1;j++){
	R[j*nm1+j]= -r;
	R[(j+1)*nm1+j]= 1;
}

// Initialize RA with 0
for(i=0;i<nm1;i++){
	for(j=0;j<n;j++){
		RA[i*n+j] = 0;
	}
}

// Initialize x with 0
for(i=0;i<n;i++){
	x[i] = 0;
}

// Initialize v with 0
for(i=0;i<nm1;i++){
	v[i] = 0;
}

// Initialize L with 0
for(i=0;i<nm1;i++){
	for(j=0;j<nm1;j++){
		L[i*nm1+j] = 0;
	}
}

// Set values to L
for(i=0;i<nm1;i++){
	L[i*nm1+i] = l;
}

// Initialize P with 0
for(i=0;i<n;i++){
	for(j=0;j<nm1;j++){
		P[i*nm1+j] = 0;
	}
}
// TODO: Fix Initializations above to match the column-major definition
  
// Initialize Q with 0
for(j=0;j<n;j++){
	for(i=0;i<n;i++){
		Q[j*n+i] = 0;
	}
}

/*---------------------------------*/
/* Computational/Algorithm Section */
/*---------------------------------*/

/* begin of time measurement*/
	s_initial = dsecnd();

//printf("Before the computation \n");

  // 1.  solve A*y'=y. y gets overwritten with the solution
dgtsv(&n, &iOne, Adl, Ad, Adu, y, &n, &info);
//printf("After 1. \n");

  // 2. compute BtB=B*B (B*B = B*transpose(B) because of the symmetric property)
dsymm(&left, &upper, &n, &n, &dOne, B, &n, B, &n, &dZero, BtB, &n);
//printf("After 2. \n");

  // 3. copmute b=BtB*y
  // TODO: Storage method needs to be changed to be able to use sbmv!
  //dsbmv(&upper, &n, &two, &dOne, BtB, &n, y, &iOne, &dZero, y, &iOne);
dsymv(&upper, &n, &dOne, BtB, &n, y, &iOne, &dZero, b, &iOne);
//printf("After 3. \n");

  // 4. compute RA=R*A
dsymm(&right, &upper, &nm1, &n, &dOne, A, &n, R, &nm1, &dZero, RA, &nm1);
//dgemm(&no, &no, &nm1, &n, &n, &dOne, R, &nm1, A, &n, &dZero, RA, &nm1);
//printf("After 4. \n");
  

// For multiple iteration the code will loop from stepp 5. to step 10.

	
  // 4bis. copy BtB onto Q
memcpy(Q,BtB,n*n*sizeof(double));

  // 5. compute v=R*x
dgemv(&no, &nm1, &n, &dOne, R, &nm1, x, &iOne, &dZero, v, &iOne);
//printf("After 5. \n");

  // 6. compute P=transpose(RA)*L
//cblas_dgemm(CblasRowMajor,&trans, &no, &n, &nm1, &nm1, &dOne, RA, &nm1, L, &nm1, &dZero, P, &n);
//cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, nm1, nm1, 1. , RA, n, L, nm1, dZero, P, nm1);
dgemm(&trans, &no, &n, &nm1, &nm1, &dOne, RA, &nm1, L, &nm1, &dZero, P, &n);
//printf("After 6. \n");
 
  
  // TODO: Set Q = copy(BtB). Use Band Storage!!!
  // 7. compute Q = P*RA+Q
dgemm(&no, &no, &n, &n, &nm1, &dOne, P, &n, RA, &nm1, &dOne, Q, &n);
//printf("After 7. \n");

  // 8. solve Q*b'=b. 
  // TODO: Adapt to be able to use dgbsv(band) instead of dgesv(sym)
//dgbsv(&n, &three, &three, &iOne, Q, &n, ipiv, b, &n, &info );
dgesv( &n, &iOne, Q, &n, ipiv, b, &n, &info );

//printf("After 8. \n");

//LAPACKE_dgesv(LAPACK_COL_MAJOR, n, 1, Q, n, ipiv, b, n);

/*
int *iter;
double *xx;
xx = (double *) malloc( 1*n * sizeof(double) );
iter = (int *) malloc(1* sizeof(int) );

LAPACKE_dsgesv (LAPACK_COL_MAJOR, n, 1, Q, n, ipiv, b, n, xx , n, iter);

//PRINTS
printf("INFO = %d \n", info);
for(i=0;i<n;i++){
	printf("IPIV[%d] = %d ", i, ipiv[i]);
}
printf("\n");

printf("ITER = %d\n",iter);
for(i=0;i<n;i++){
	printf("xx[%d] = %d ", i, xx[i]);
}
*/

//dpbsv(&upper, &n, &three, &iOne, Q, &n, b, &n, &info); THIS ONLY WORKS FOR SPD
//dsysv(&upper, &n, &iOne, Q, &n, &ipiv, b, &n, );
  // TODO: set x=copy(b)?? to check

  // 9. Compute x=A*b
//dgbmv(&no, &n, &n, &iOne, &iOne, &dOne, A, &n, b, &iOne, &dZero, x, &iOne);
dsymv(&upper, &n, &dOne, A, &n, b, &iOne, &dZero, x, &iOne);

//printf("\nAfter computation \n");
//printf("\n");

/* end of time measurement*/
	s_elapsed = (dsecnd() - s_initial);


/*---------------------------------------------*/
/* Print every matrix and vector defined above */
/*---------------------------------------------*/

/*

// Print A
for(i=0;i<n;i++){
	for(j=0;j<n;j++){
		printf("A[%d] = %f ",j*n+i,A[j*n+i]);
	}
printf("\n");
}

printf("\n");
  
// Print B
for(i=0;i<n;i++){
	for(j=0;j<n;j++){
		printf("B[%d] = %f ",j*n+i,B[j*n+i]);
	}
printf("\n");
}

printf("\n");

// Print BtB
for(i=0;i<n;i++){
	for(j=0;j<n;j++){
		printf("BtB[%d] = %f ",j*n+i,BtB[j*n+i]);
	}
printf("\n");
}

printf("\n");
  
// Print L
for(i=0;i<nm1;i++){
	for(j=0;j<nm1;j++){
		printf("L[%d] = %f ",j*nm1+i,L[j*nm1+i]);
	}
printf("\n");
}

printf("\n");  
  
// Print P
for(i=0;i<n;i++){
	for(j=0;j<nm1;j++){
		printf("P[%d] = %f ",j*n+i,P[j*n+i]);
	}
printf("\n");
}

printf("\n");


// Print Q
for(i=0;i<n;i++){
	for(j=0;j<n;j++){
		printf("Q[%d] = %f ",j*n+i,Q[j*n+i]);
	}
printf("\n");
}

printf("\n");  

// Print R
for(i=0;i<nm1;i++){
	for(j=0;j<n;j++){
		printf("R[%d] = %f ",j*nm1+i,R[j*nm1+i]);
	}
printf("\n");
}

printf("\n");

// Print RA
for(i=0;i<nm1;i++){
	for(j=0;j<n;j++){
		printf("RA[%d] = %f ",j*nm1+i,RA[j*nm1+i]);
	}
printf("\n");
}

printf("\n");

// Print Ad
for(i=0;i<n;i++){
	printf("Ad[%d] = %f ",i,Ad[i]);
}

printf("\n");
printf("\n");    

// Print Adl
for(i=0;i<nm1;i++){
	printf("Adl[%d] = %f ",i,Adl[i]);
}

printf("\n");
printf("\n");    

// Print Adu
for(i=0;i<nm1;i++){
	printf("Adu[%d] = %f ",i,Adu[i]);
}

printf("\n");
printf("\n"); 
  
// Print b
for(i=0;i<n;i++){
	printf("b[%d] = %f ",i,b[i]);
}

printf("\n");
printf("\n");    
  
// Print v
for(i=0;i<nm1;i++){
	printf("v[%d] = %f ",i,v[i]);
}

printf("\n");
printf("\n");  

// Print x
for(i=0;i<n;i++){
	printf("x[%d] = %f ",i,x[i]);
}

printf("\n");
printf("\n");

// Print y
for(i=0;i<n;i++){
	printf("y[%d] = %f ",i,y[i]);
}

printf("\n");
printf("\n");

// Print ycpy
for(i=0;i<n;i++){
	printf("ycpy[%d] = %f ",i,ycpy[i]);
}
 
printf("\n");
printf("\n");

*/

printf("%.5f seconds \n\n", (s_elapsed));

/*---------------------------------------------------------------*/
/* Free up all the matrices and vectors defined for this problem */
/*---------------------------------------------------------------*/

free(A);	
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
free(x);
free(v);

/* END OF PROGRAM */

}
