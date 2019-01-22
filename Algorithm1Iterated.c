#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"
#include <string.h>

int main(int argc, char *argv[]){


#ifndef max
    #define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
    #define min(a,b) ((a) < (b) ? (a) : (b))
#endif
  
  // Time measurement variables
  double s_initial,s_elapsed;
  // Iterators
  int i,j,k;	
  int n;

  // Initial matrices
  double *A, *Ad, *Adu, *Adl, *Adcpy, *Aducpy, *Adlcpy, *B, *L, *R;
  // Initial vectors
  double *x, *y, *ycpy;
  // Initial parameters and iterator
  double l,r;
  int iterations;

  int threads = 32;
  n = 4000; 		// Problem Size
  iterations = 1;	// Number of iterations
  r = 0.45;		// Parameter for Matrix R
  l = 0.75;	 	// Paramter for Matrix L TODO: Check if randomizing l changes anything


  // Computed Matrices
  double *BtB, *P, *Q, *RA;
  double *ABand, *BBand, *BtBBand, *QBand;
  // Computer vectors
  double *b, *bcpy, *f, *v;


  /* creation of parameters for BLAS routines */
  char lower = 'L';
  char upper = 'U';
  char trans = 'T';
  char no = 'N';
  char yes ='Y';
  char left ='L';
  char right = 'R';
  double mdOne = -1.;
  int iZero = 0;
  double dZero = 0.;
  int iOne = 1;
  double dOne = 1.;
  int two = 2;
  int three = 3;
  int five = 5;
  int info;
  int *ipiv;

  int nm1 = n-1;

  

  // Allocate initial matrices and vectors	
  A = (double *) malloc( n*n * sizeof(double) );
  B = (double *) malloc( n*n * sizeof(double) );
  L = (double *) malloc( nm1*nm1 * sizeof(double) );
  R = (double *) malloc( n*nm1 * sizeof(double) );
  x = (double *) malloc( 1*n * sizeof(double) );
  y = (double *) malloc( 1*n * sizeof(double) );
  ycpy = (double *) malloc( 1*n * sizeof(double) );


  // Allocate 3 arrays for the tridiagonal matrix A
  Ad = (double *) malloc( 1*n * sizeof(double) );       //A diagonal
  Adu = (double *) malloc( 1*nm1 * sizeof(double) );  // A upper diagonal
  Adl = (double *) malloc( 1*nm1 * sizeof(double) );  // A lower diagonal
  Adcpy = (double *) malloc( 1*n * sizeof(double) );       //A diagonal
  Aducpy = (double *) malloc( 1*nm1 * sizeof(double) );	// A upper diagonal
  Adlcpy = (double *) malloc( 1*nm1 * sizeof(double) );	// A lower diagonal


  //Matrices to store computed results
  BtB = (double *) malloc( n*n * sizeof(double) );
  P = (double *) malloc( nm1*n * sizeof(double) );
  Q = (double *) malloc( n*n * sizeof(double) );
  RA = (double *) malloc( n*nm1 * sizeof(double) );

  ABand = (double *) malloc( n*3 * sizeof(double) );
  BBand = (double *) malloc( n*2 * sizeof(double) );
  BtBBand = (double *) malloc( n*3 * sizeof(double) );
  QBand = (double *) malloc( n*10 * sizeof(double) );


  //Vectors to store computed results
  b = (double *) malloc( 1*n * sizeof(double) );
  bcpy = (double *) malloc( 1*n * sizeof(double) );
  f = (double *) malloc( 1*n * sizeof(double) );
  v = (double *) malloc( 1*nm1 * sizeof(double) );


  ipiv = (int *) malloc( 1*n*sizeof(int) );


  srand( time(NULL) );
  double random = (double)rand()/(double)RAND_MAX;

	//printf("random =%f \n", random);

// Initialize ipiv
for(i=0;i<n;i++){
	ipiv[i] = 0;
}


// Set values of Ad,Adu,Adl
for(j=0;j<n-1;j++){
	Ad[j] = 2.;//(double)rand()/(double)RAND_MAX;
	Adu[j] = 3.;//(double)rand()/(double)RAND_MAX;
	Adl[j] = Adu[j];
}
Ad[n-1]= 2.;//(double)rand()/(double)RAND_MAX;

memcpy(Adcpy,Ad,n*1*sizeof(double));
memcpy(Aducpy,Adu,nm1*1*sizeof(double));
memcpy(Adlcpy,Adl,nm1*1*sizeof(double));


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

// Initialize ABand with 0
for(j=0;j<n;j++){
	for(i=0;i<3;i++){
		ABand[i*3+j] = 0;
	}
}

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

// Initialize BtBBand with 0
for(j=0;j<n;j++){
	for(i=0;i<3;i++){
		BtBBand[i*3+j] = 0;
	}
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

// Initialize QBand with 0
for(j=0;j<n;j++){
	for(i=0;i<10;i++){
		QBand[i*10+j] = 0;
	}
}

// Initialize R with 0
for(j=0;j<n;j++){
	for(i=0;i<nm1;i++){
		R[j*nm1+i] = 0;
	}
}
  
// Set values to R using parameter r9912
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

// Initialize b with 0
for(j=0;j<n;j++){
	b[j]=0;
}

// Initialize f with 0
for(i=0;i<n;i++){
	f[i] = 0;
}

// Initialize v with 0
for(i=0;i<nm1;i++){
	v[i] = 0;
}

// Initialize x with 0
for(i=0;i<n;i++){
	x[i] = 0;
}

// Initialize and set values of y
for(i=0;i<n;i++){
	//y[i] = (double)rand()/(double)RAND_MAX;
	y[i] = 0.5;
}

// Copy y into ycpy
memcpy(ycpy,y,n*sizeof(double));


/*-------------------------*/
/* Pre-Computation Section */
/*-------------------------*/


 //Converts upper triangle of full rank matrix A into band storage matrix ABand
for (i = nm1; i >=0; i--) {
    for (j = 0; j < min(2, (i+1)); j++) {
	ABand[i*3+(2-j)-1] = A[i*n+(i-j)];
    }
}
//Converts lower triangle of full rank matrix A into band storage matrix ABand
for(i=0;i<nm1;i++){
	ABand[3*i+2]=ABand[3*(i+1)];
}



 //Converts upper triangle of full rank matrix B into band storage matrix BBand
for (i = nm1; i >=0; i--) {
    for (j = 0; j < min(2, (i+1)); j++) {
	BBand[i*2+(1-j)] = B[i*n+(i-j)];
    }
}


/*---------------------------------*/
/* Computational/Algorithm Section */
/*---------------------------------*/
mkl_set_num_threads(threads);

//begin of time measurement
s_initial = dsecnd();

//printf("Before the computation \n");

  // 1.  solve A*y'=y. y gets overwritten with the solution
dgtsv(&n, &iOne, Adl, Ad, Adu, y, &n, &info);
//printf("After 1. \n");

  // 2. compute BtB=B*B (B*B = B*transpose(B) because of the symmetric property)
dsymm(&left, &upper, &n, &n, &dOne, B, &n, B, &n, &dZero, BtB, &n);
//dsyrk(&upper, &no, &n, &n, &dOne, B, &n, &dOne, BtB, &n);
//printf("After 2. \n");




//Store BtB in BtBBand using the Band Storage Method
//BtBBand[i*3 + j] = BtB[i*n+(n-j)];

/*
//Converts upper triangle of a full rank matrix into a band storage matrix
for (i = nm1; i >=0; i--) {
//printf("Outside: i=%d\n", i);
    for (j = 0; j < min(3, (i+1)); j++) {
	//printf("BtB[%d] = %f\n",i*n+(i-j),BtB[i*n+(i-j)]);
	BtBBand[i*3+(3-j)-1] = BtB[i*n+(i-j)];
	//printf("Inside: j=%d\n",j);
	//printf("i*3+(3-j)-1=%d ",i*3+(3-j)-1);
	//printf("BtB[%d]=%f \n",i*n+(i-j),BtB[i*n+(i-j)]);
    }
}
*/

  // 3. copmute b=BtB*y and create a copy of b (bcpy)
  // TODO: Storage method needs to be changed to be able to use sbmv!
//dsbmv(&upper, &n, &two, &dOne, BtBBand, &three, y, &iOne, &dZero, b, &iOne);
dsymv(&upper, &n, &dOne, BtB, &n, y, &iOne, &dZero, b, &iOne);
memcpy(bcpy,b,n*sizeof(double));

//printf("After 3. \n");

  // 4. compute RA=R*A
dsymm(&right, &upper, &nm1, &n, &dOne, A, &n, R, &nm1, &dZero, RA, &nm1);
//dgemm(&no, &no, &nm1, &n, &n, &dOne, R, &nm1, A, &n, &dZero, RA, &nm1);
//printf("After 4. \n");
  

// For multiple iteration the code will loop from stepp 5. to step 10.
for(k=0;k<iterations;k++){

	  // 4bis. copy BtB onto Q
	  //       copy bcpy onto b
	memcpy(Q,BtB,n*n*sizeof(double));
	memcpy(b,bcpy,n*sizeof(double));

	  // 5. compute v=R*x
	dgemv(&no, &nm1, &n, &dOne, R, &nm1, x, &iOne, &dZero, v, &iOne);
	//printf("After 5. \n");


	  // 6. compute P=transpose(RA)*L
	dgemm(&trans, &no, &n, &nm1, &nm1, &dOne, RA, &nm1, L, &nm1, &dZero, P, &n);
	//printf("After 6. \n");

	 
	  // TODO: Set Q = copy(BtB). Use Band Storage!!!
	  // 7. compute Q = P*RA+Q
	dgemm(&no, &no, &n, &n, &nm1, &dOne, P, &n, RA, &nm1, &dOne, Q, &n);
	//printf("After 7. \n");


/*
	 //Converts upper triangle of a full rank matrix into a band storage matrix
	for (i = nm1; i >=0; i--) {
	    for (j = 0; j < min(4, (i+1)); j++) {
		QBand[i*10+(7-j)-1] = Q[i*n+(i-j)];
	    }
	}
	// Converts the lower tringular part into the band storage matrix
	for(j=0;j<n;j++){
		//printf("Outside: %d\n",j);
		for(i=0;i<min(3,n-j-1);i++){
			//printf("Inside: %d\n",i);
			QBand[j*10+7+i]=QBand[(j+1+i)*10+(5-i)];
			//printf("j*7+5+i = %d\t (j+i+1)*7+(3-i)=%d\n",j*7+5+i,(j+i+1)*7+(3-i));	
		}
	}


	int ldab = 10;
*/
	// 8. solve Q*b'=b. 
	//dgbsv(&n, &three, &three, &iOne, QBand, &ldab, ipiv, b, &n, &info );
	dgesv( &n, &iOne, Q, &n, ipiv, b, &n, &info );
	//dpbsv(&upper, &n, &three, &iOne, Q, &n, b, &n, &info); THIS ONLY WORKS FOR SPD
	//dsysv(&upper, &n, &iOne, Q, &n, &ipiv, b, &n, &info);

	//printf("After 8. \n");

	  // 9. Compute x=A*b
	//dgbmv(&no, &n, &n, &iOne, &iOne, &dOne, ABand, &three, b, &iOne, &dZero, x, &iOne);
	dsymv(&upper, &n, &dOne, A, &n, b, &iOne, &dZero, x, &iOne);
}

//printf("\nAfter computation \n");
//printf("\n");

//10. compute ycpy = ycpy-x
daxpy(&n, &mdOne, x, &iOne, ycpy, &iOne); 

//11.  solve A*ycpy'=ycpy. y gets overwritten with the solution
dgtsv(&n, &iOne, Adlcpy, Adcpy, Aducpy, ycpy, &n, &info);


memcpy(f,ycpy,n*1*sizeof(double));


//12. TODO: SBMV  NOT WORKING YET, CHECK 
//dsbmv(&upper, &n, &iOne, &mdOne, BBand, &two, ycpy, &iOne, &dOne, f, &iOne);
dsymv(&upper, &n, &mdOne, B, &n, ycpy, &iOne, &dOne, f, &iOne);

//end of time measurement
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

// Print ABand
for(i=0;i<3;i++){
	for(j=0;j<n;j++){
		printf("ABand[%d] = %f ",j*3+i,ABand[j*3+i]);
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
  


// Print Adcpy
for(i=0;i<n;i++){
	printf("Adcpy[%d] = %f ",i,Adcpy[i]);
}

printf("\n");
printf("\n");    

// Print Adlcpy
for(i=0;i<nm1;i++){
	printf("Adlcpy[%d] = %f ",i,Adlcpy[i]);
}

printf("\n");
printf("\n");    

// Print Aducpy
for(i=0;i<nm1;i++){
	printf("Aducpy[%d] = %f ",i,Aducpy[i]);
}

printf("\n");
printf("\n");
  
// Print B
for(i=0;i<n;i++){
	for(j=0;j<n;j++){
		printf("B[%d] = %f ",j*n+i,B[j*n+i]);
	}
printf("\n");
}

printf("\n");

// Print BBand
for(i=0;i<2;i++){
	for(j=0;j<n;j++){
		printf("BBand[%d] = %f ",j*2+i,BBand[j*2+i]);
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
 
// Print BtBBand
for(i=0;i<3;i++){
	for(j=0;j<n;j++){
		printf("BtBB[%d] = %f ",j*3+i,BtBBand[j*3+i]);
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

// Print QBand
for(i=0;i<10;i++){
	for(j=0;j<n;j++){
		printf("QB[%d] = %f ",j*10+i,QBand[j*10+i]);
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


// Print b
for(i=0;i<n;i++){
	printf("b[%d] = %f ",i,b[i]);
}

printf("\n");
printf("\n");    
 

// Print f
for(i=0;i<n;i++){
	printf("f[%d] = %f ",i,f[i]);
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

printf("Computation time: %.5f seconds \nSize n: %d \nIterations : %d", (s_elapsed), n, iterations);

/*---------------------------------------------------------------*/
/* Free up all the matrices and vectors defined for this problem */
/*---------------------------------------------------------------*/

free(A);
free(Ad);
free(Adu);
free(Adl);
free(Adcpy);
free(Aducpy);
free(Adlcpy);
free(ABand);
free(B);
free(BBand);
free(BtB);
free(BtBBand);
free(L);
free(Q);
free(QBand);
free(R);
free(RA);
free(P);
free(b);
free(bcpy);
free(f);
free(v);
free(x);
free(y);
free(ycpy);



/* END OF PROGRAM */

}
