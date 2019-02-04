#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"
#include <string.h>

int main(int argc, char *argv[]){

/****************************/
/*Help Functions Definitions*/
/****************************/
#ifndef max
    #define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
    #define min(a,b) ((a) < (b) ? (a) : (b))
#endif
  


void sbmmcompact(double bd,double bdu,double *BtBc){
	BtBc[5*0+2]=pow(bd,2)+pow(bdu,2);
	BtBc[5*1+2]=pow(bd,2)+2*pow(bdu,2);

	BtBc[5*0+0]=pow(bdu,2);
	BtBc[5*1+0]=pow(bdu,2);
	BtBc[5*1+4]=pow(bdu,2);

	BtBc[5*0+1]=2*bd*bdu;
	BtBc[5*1+1]=2*bd*bdu;
	BtBc[5*1+3]=2*bd*bdu;
}

void computeRAcompact(double r, double ad, double adu, double *RAc){
	RAc[0] = adu;
	RAc[1] = -r*adu+ad;
   	RAc[2] = -r*ad+adu;
    	RAc[3] = -r*adu;
}


void computeLRAcompact(double *L, double *RAc,int n, int nm1, double *LRAc){
int i;

for(i=0;i<(n-2);i++){
	LRAc[(i+1)*4+1] = L[i*nm1+i]*RAc[1];
	LRAc[i*4+2] = L[i*nm1+i]*RAc[2];
	LRAc[(i+2)*4+0] = L[i*nm1+i]*RAc[0];
        LRAc[i*4+3] = L[(i+1)*nm1+(i+1)]*RAc[3];
} 
LRAc[nm1*4+1] = L[(i-1)*nm1+(i-1)]*RAc[1];
LRAc[(n-2)*4+2] = L[(i-1)*nm1+(i-1)]*RAc[2];

}

void sbmvcompact(double *Bc, double *y, int n, int nm1, double *b){
int i,k;

for(i=1;i<nm1;i++){
	for(k=max(0,2-i);k<=min(4,n+1-i);k++){
		b[i] += Bc[5*1+k]*y[i-2+k];	
		//printf("i: %d \t k : %d \t 5*1+k : %d \t i-2+k : %d \n",i,k,5*1+k,i-2+k);
	}
}

for(k=0;k<3;k++){
	b[0] += Bc[5*0+2-k]*y[k];
	b[nm1] += Bc[5*0+k]*y[n-3+k];
}


}


void sbcompactstorage(int ad, int adu, double *Ac){

Ac[5*0+2]=ad;
Ac[5*1+2]=ad;

Ac[5*0+1]=adu;
Ac[5*1+1]=adu;
Ac[5*1+3]=adu;

}

void addmBand(double *Bc, double *Ac, int n, int nm1, double *Rc){

int i,j,k;
memcpy(Rc,Ac,n*n*sizeof(double));
//*R = *A; //TODO: can the sum be applied to A right away or are the original values of A required later on? 

for(i=0;i<n;i++){
	for(k=max(0,3-i);k<=min(4,n-i);k++){
		Rc[i*10+4+k] = Ac[i*10+4+k]+Bc[5*1+k];
		//printf("i = %d \t k = %d \n",i,k);
	}
}

for(j=0;j<3;j++){
	Rc[(2-j)*10+4+j] = Ac[(2-j)*10+4+j]+Bc[0*5+j];
	Rc[(n-3+j)*10+8-j]= Ac[(n-3+j)*10+8-j]+Bc[0*5+j];
	//printf("(n-3+j)*n+nm1-1 = %d \n",(n-3+j)*n+nm1-1);
}

}


void computeQBandcompact(double *LRAc, int n, int nm1, double *Q){

int i,j,k,count = 4,band=0;

//Upper triangular
for(i=0;i<n;i++){
	for(j=i;j<=min(nm1,i+3);j++){
		for(k=1;k<=count;k++){
			//printf("i=%d    j=%d    k=%d \n",i,j,k);
			Q[j*10+(6-band)] += LRAc[i*4+(3-count+k)]*LRAc[j*4+k-1];
		}
	count--;
	band++;
	}
count = 4;
band = 0;
}


//Lower triangular

count = 3;
for(i=nm1;i>=0;i--){
	for(j=i-1;j>=max(0,i-3);j--){
		for(k=1;k<=count;k++){
			//printf("i=%d    j=%d    k=%d    count=%d    3-count+k=%d  j*n+i=%d   LRAc[%d]=%f \n",i,j,k,count,3-count+k,j*n+i,(i-1)*4+(3-count+k),LRAc[(i-1)*4+(3-count+k)]);
			Q[j*10+7+band] += LRAc[j*4+(3-count+k)]*LRAc[i*4+k-1];
		}
	count--;
	band++;
	}
count = 3;
band = 0;
}

}

/***********/
/*Variables*/
/***********/

  // Time measurement variables
  double s_initial,s_elapsed;
  // Iterators
  int i, j, k, iterations;	
  int n, nm1;
  // Initial parameters and iterator setter
  double l,r;

  // Initial matrices
  double *A, *Ad, *Adu, *Adl, *Adcpy, *Aducpy, *Adlcpy, *B, *L, *R;
  // Initial vectors
  double *x, *y, *ycpy;


/*-----INPUT PARAMETERS-----*/
  int threads = 1;
  n = 10; 		// Problem Size
  nm1 = n-1;
  iterations = 1;	// Number of iterations
  r = 0.45;		// Parameter for Matrix R
  l = 0.75;	 	// Paramter for Matrix L TODO: Check if randomizing l changes anything
 double ad = 2,adu = 3;
 double bd = 4,bdu = 5;


  // Computed Matrices
  double *Ac, *Bc, *BtB, *BtBc, *LRAc, *P, *Q, *RA, *RAc;
  double *ABand, *BBand, *BtBBand, *QBand;
  // Computer vectors
  double *b, *bcpy, *dyx, *f, *v;


  /* BLAS/LAPACK PARAMETERS */
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


  // Allocate initial matrices and vectors	
  A = (double *) mkl_malloc( n*n * sizeof(double) ,64 );
  B = (double *) mkl_malloc( n*n * sizeof(double) ,64 );
  L = (double *) mkl_malloc( nm1*nm1 * sizeof(double) ,64 );
  R = (double *) mkl_malloc( n*nm1 * sizeof(double) ,64 );
  x = (double *) mkl_malloc( 1*n * sizeof(double) ,64 );
  y = (double *) mkl_malloc( 1*n * sizeof(double) ,64 );
  ycpy = (double *) mkl_malloc( 1*n * sizeof(double) ,64 );


  // Allocate 3 arrays for the tridiagonal matrix A
  Ad = (double *) mkl_malloc( 1*n * sizeof(double) ,64 );       //A diagonal
  Adu = (double *) mkl_malloc( 1*nm1 * sizeof(double) ,64 );  // A upper diagonal
  Adl = (double *) mkl_malloc( 1*nm1 * sizeof(double) ,64 );  // A lower diagonal
  Adcpy = (double *) mkl_malloc( 1*n * sizeof(double) ,64 );       //A diagonal
  Aducpy = (double *) mkl_malloc( 1*nm1 * sizeof(double) ,64 );	// A upper diagonal
  Adlcpy = (double *) mkl_malloc( 1*nm1 * sizeof(double) ,64 ); // A lower diagonal


  //Matrices to store computed results
  BtB = (double *) mkl_malloc( n*n * sizeof(double) ,64 );
  P = (double *) mkl_malloc( nm1*n * sizeof(double) ,64 );
  Q = (double *) mkl_malloc( n*n * sizeof(double) ,64 );
  RA = (double *) mkl_malloc( n*nm1 * sizeof(double) ,64 );

  ABand = (double *) mkl_malloc( n*3 * sizeof(double) ,64 );
  BBand = (double *) mkl_malloc( n*2 * sizeof(double) ,64 );
  BtBBand = (double *) mkl_malloc( n*3 * sizeof(double) ,64 );
  QBand = (double *) mkl_malloc( n*10 * sizeof(double) ,64 );

  BtBc = (double *) mkl_malloc( 2*5 * sizeof(double) ,64 );
  RAc = (double *) mkl_malloc( 1*4 * sizeof(double) ,64 );
  LRAc = (double *) mkl_malloc( n*4 * sizeof(double) ,64 );
  Ac = (double *) mkl_malloc( 2*5 * sizeof(double) ,64 );
  Bc = (double *) mkl_malloc( 2*5 * sizeof(double) ,64 );

  //Vectors to store computed results
  b = (double *) mkl_malloc( 1*n * sizeof(double) ,64 );
  bcpy = (double *) mkl_malloc( 1*n * sizeof(double) ,64 );
  f = (double *) mkl_malloc( 1*n * sizeof(double) ,64 );
  v = (double *) mkl_malloc( 1*nm1 * sizeof(double) ,64 );
  dyx = (double *) mkl_malloc( n*1 * sizeof(double) ,64 );

  ipiv = (int *) mkl_malloc( 1*n*sizeof(int) ,64 );



/**********************************************************************/
  // RANDOMIZER
  srand( time(NULL) );
  double random = (double)rand()/(double)RAND_MAX;
  //printf("random =%f \n", random);

/*--------Variables Initialisation-------*/


// Set values of Ad,Adu,Adl
for(j=0;j<n-1;j++){
	Ad[j] = ad;//(double)rand()/(double)RAND_MAX;
	Adu[j] = adu;//(double)rand()/(double)RAND_MAX;
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

// Initialize Ac with 0
for(j=0;j<2;j++){
	for(i=0;i<5;i++){
		Ac[j*4+i] = 0;
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

// Initialize Bc with 0
for(j=0;j<2;j++){
	for(i=0;i<5;i++){
		Bc[j*4+i] = 0;
	}
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

// Initialize BtBc with 0
for(j=0;j<5*2;j++){
	BtBc[j] = 0.;
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

// Initialize LRAc with 0
for(j=0;j<n;j++){
	for(i=0;i<4;i++){
		LRAc[j*4+i] = 0;
	}
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

// Initialize RAc with 0
for(j=0;j<4;j++){
	RAc[j] = 0.;
}

// Initialize b with 0
for(j=0;j<n;j++){
	b[j]=0;
}

// Initialize dyx with 0
for(i=0;i<n;i++){
	dyx[i] = 0;
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

// Initialize ipiv
for(i=0;i<n;i++){
	ipiv[i] = 0;
}




/***************************/
/* Pre-Computation Section */
/***************************/

sbcompactstorage(ad,adu,Ac);
sbcompactstorage(bd,bdu,Bc);


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

  // 1.  solve A*y'=y. y gets overwritten with the solution
dgtsv(&n, &iOne, Adl, Ad, Adu, y, &n, &info);


  // 2. compute BtB=B*B (B*B = B*transpose(B) because of the symmetric property)
//dsymm(&left, &upper, &n, &n, &dOne, B, &n, B, &n, &dZero, BtB, &n);
//dsyrk(&upper, &no, &n, &n, &dOne, B, &n, &dOne, BtB, &n);
sbmmcompact(bd,bdu,BtBc);

  // 3. copmute b=BtB*y and create a copy of b (bcpy)
  // TODO: Storage method needs to be changed to be able to use sbmv!
//dsbmv(&upper, &n, &two, &dOne, BtBBand, &three, y, &iOne, &dZero, b, &iOne);
//dsymv(&upper, &n, &dOne, BtB, &n, y, &iOne, &dZero, b, &iOne);
sbmvcompact(BtBc,y,n,nm1,b);
  
memcpy(bcpy,b,n*sizeof(double));

  // 4. compute RA=R*A
//dsymm(&right, &upper, &nm1, &n, &dOne, A, &n, R, &nm1, &dZero, RA, &nm1);
//dgemm(&no, &no, &nm1, &n, &n, &dOne, R, &nm1, A, &n, &dZero, RA, &nm1);
computeRAcompact(r,ad,adu,RAc);


// For multiple iteration the code will loop from stepp 4bis. to step 10.
for(k=0;k<iterations;k++){

	  // 4bis. copy BtB onto Q
	  //       copy bcpy onto b
	//memcpy(Q,BtB,n*n*sizeof(double));
	memcpy(b,bcpy,n*sizeof(double));

	  // 5. compute v=R*x
	//dgemv(&no, &nm1, &n, &dOne, R, &nm1, x, &iOne, &dZero, v, &iOne);
	for(i=0;i<nm1;i++){
		v[i]=-r*x[i]+x[i+1];
	}	

	  // 6. compute P=transpose(RA)*L
	//dgemm(&trans, &no, &n, &nm1, &nm1, &dOne, RA, &nm1, L, &nm1, &dZero, P, &n);
	computeLRAcompact(L, RAc, n, nm1, LRAc);

	 
	  // TODO: gemm + seperate addition to avoid memcopy ?
	  // 7. compute Q = P*RA+Q
	//dgemm(&no, &no, &n, &n, &nm1, &dOne, P, &n, RA, &nm1, &dOne, Q, &n);
	//computeQcompact(LRAc, n, nm1, Q);
	computeQBandcompact(LRAc, n, nm1, QBand);

	// 7bis. compute Q = Q+BtB
	addmBand(BtBc,QBand,n,nm1,QBand);

	int ldab = 10;
	// 8. solve Q*b'=b. 
	//dgesv( &n, &iOne, Q, &n, ipiv, b, &n, &info );
	//dsysv(&upper, &n, &iOne, Q, &n, &ipiv, b, &n, );

	//dgbsv(&n, &three, &three, &iOne, QBand, &ldab, ipiv, b, &n, &info );
	dgbtrf(&ldab, &n, &three, &three, QBand, &ldab, ipiv, &info);
	dgbtrs(&no, &n, &three, &three, &iOne, QBand, &ldab, ipiv, b, &n,&info);


	  // 9. Compute x=A*b
	//dgbmv(&no, &n, &n, &iOne, &iOne, &dOne, ABand, &three, b, &iOne, &dZero, x, &iOne);
	//dsymv(&upper, &n, &dOne, A, &n, b, &iOne, &dZero, x, &iOne);
	sbmvcompact(Ac, b, n, nm1, x);
}

//10. compute yinit = yinit-x and store a copy in dyx
daxpy(&n, &mdOne, x, &iOne, ycpy, &iOne); 

memcpy(dyx,ycpy,n*sizeof(double));


//11.  solve A*ycpy'=ycpy. y gets overwritten with the solution
dgtsv(&n, &iOne, Adlcpy, Adcpy, Aducpy, ycpy, &n, &info);

  //12. compute f = Bc*ycpy
//TODO: SBMV  NOT WORKING YET, CHECK 
//dsbmv(&upper, &n, &iOne, &mdOne, BBand, &two, ycpy, &iOne, &dOne, f, &iOne);
//dsymv(&upper, &n, &mdOne, B, &n, ycpy, &iOne, &dOne, f, &iOne);
sbmvcompact(Bc, ycpy, n, nm1, f);

  //13. compute f = dyx - f
for(i=0;i<n;i++){
	f[i]=dyx[i]-f[i];
}

// Print f
for(i=0;i<n;i++){
	printf("f[%d] = %f ",i,f[i]);
}

printf("\n");
printf("\n");

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

/*
//Print Ac
for(i=0;i<5;i++){
	for(j=0;j<2;j++){
		printf("Ac[%d] = %f ",j*5+i,Ac[j*5+i]);
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

//Print Bc
for(i=0;i<5;i++){
	for(j=0;j<2;j++){
		printf("Bc[%d] = %f ",j*5+i,Bc[j*5+i]);
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

// Print BtBc
for(i=0;i<5;i++){
	for(j=0;j<2;j++){
		printf("BtBc[%d] = %f ",j*5+i,BtBc[j*5+i]);
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

//Print LRAc
for(i=0;i<4;i++){
	for(j=0;j<n;j++){
		printf("LRAc[%d] = %f ",j*4+i,LRAc[j*4+i]);
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

// Print RAc
for(i=0;i<4;i++){
		printf("RAc[%d] = %f ",i,RAc[i]);
printf("\n");
}
printf("\n");

// Print b
for(i=0;i<n;i++){
	printf("b[%d] = %f ",i,b[i]);
}

printf("\n");
printf("\n");    
 
// Print dyx
for(i=0;i<n;i++){
	printf("dyx[%d] = %f ",i,dyx[i]);
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

printf("Computation time: %.5f seconds \nSize n: %d \nIterations: %d \n", (s_elapsed), n, iterations);

/*---------------------------------------------------------------*/
/* Free up all the matrices and vectors defined for this problem */
/*---------------------------------------------------------------*/

mkl_free(A);
mkl_free(Ac);
mkl_free(Ad);
mkl_free(Adu);
mkl_free(Adl);
mkl_free(Adcpy);
mkl_free(Aducpy);
mkl_free(Adlcpy);
mkl_free(ABand);
mkl_free(B);
mkl_free(BBand);
mkl_free(Bc);
mkl_free(BtB);
mkl_free(BtBBand);
mkl_free(BtBc);
mkl_free(L);
mkl_free(LRAc);
mkl_free(Q);
mkl_free(QBand);
mkl_free(R);
mkl_free(RA);
mkl_free(RAc);
mkl_free(P);
mkl_free(b);
mkl_free(bcpy);
mkl_free(dyx);
mkl_free(f);
mkl_free(v);
mkl_free(x);
mkl_free(y);
mkl_free(ycpy);

}
