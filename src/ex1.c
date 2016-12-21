#include <slepceps.h>
#include <petsctime.h>
#include <petscvec.h>
#include "stimsyr.h"
#include "read2.h"

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscReal      norm;
  PetscScalar    kr,ki,value[700];
  Vec            xr,xi;
  PetscInt       i,Istart,Iend,col[700],maxit,its,nconv,countcol;
  PetscInt		 nev, ncv, mpd;
  PetscLogDouble t1,t2,tt1,tt2;
//PetscBool      FirstBlock=PETSC_FALSE,LastBlock=PETSC_FALSE;
  PetscErrorCode ierr;
//PetscScalar    eigr;
//PetscScalar    eigi;
  int            mpiid;
  int            natomax=700;
  int			 ideter[natomax];

  char const* const fileName = argv[1];
  FILE* file = fopen(fileName, "r");
  Data getdata;
  PetscInt		 nlocal;
  
  Data_new(file, &getdata);
  nlocal = getdata.n/getdata.npar;
//printf("n=%ld\t nnz=%ld\t npar=%ld\t ntrou=%ld\t isz=%ld\n",getdata.n,getdata.nnz,getdata.npar,getdata.ntrou,getdata.isz);


  PetscScalar	 valxr[nlocal];
  PetscInt   	 indxr[nlocal];
//Vec            Vr,Vi;
  char           filename[PETSC_MAX_PATH_LEN]="FIL666";
  PetscViewer    viewer;
  PetscBool      ishermitian;
  PetscInt       kk,ll,iii2;
  long int       iii;
  long int		 ii;
  long int       tcountcol2,tcol[700],tcountcol[getdata.nnz];
  double         val[700];

  SlepcInitialize(&argc,&argv,(char*)0,NULL);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1-D t-J Eigenproblem, n=%D\n\n",getdata.n);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,getdata.n,getdata.n,getdata.nnz*getdata.npar,NULL,getdata.nnz*getdata.npar,NULL,&A);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,getdata.nnz*getdata.npar,NULL,getdata.nnz*getdata.npar,NULL);CHKERRQ(ierr);
//ierr = MatSetFromOptions(A);CHKERRQ(ierr);
//ierr = MatSetUp(A);CHKERRQ(ierr);

  MPI_Comm_rank(MPI_COMM_WORLD,&mpiid);
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  ierr = PetscTime(&tt1);CHKERRQ(ierr);

  for (i=Istart; i<Iend; i+=getdata.nnz) {
      tcountcol2=0;
      for(kk=0;kk<getdata.nnz;kk++){
          tcountcol[kk]=0;
      }
      iii=i+1;
      if(i%getdata.npar == 0 && mpiid==0){
        ierr = PetscTime(&t1);CHKERRQ(ierr);
      }
      unit_l1_(
            getdata.l1,
            getdata.l2,
            getdata.ktyp,
            &iii,
			&getdata.nnz,
            getdata.xjjxy,
            getdata.xjjz ,
            getdata.xtt ,
            tcountcol,
            &getdata.ntrou,
            &getdata.isz,
            tcol,
            val);
     if(i%getdata.npar == 0 && mpiid==0){
       ierr = PetscTime(&t2);CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_WORLD," i: %d\n mpiid: %d\ntime: %f\n",i,mpiid,t2-t1);CHKERRQ(ierr);
     }
      for(ll=0;ll<getdata.nnz;ll++){

      for(kk=0;kk<tcountcol[ll]+1;kk++){
          value[kk] = val[kk+tcountcol2];
          col[kk] = tcol[kk+tcountcol2]-1;
//        PetscPrintf(PETSC_COMM_WORLD,"value = %f col = %d\n",value[kk],col[kk]);
      }
      for(kk=tcountcol2+tcountcol[ll]+1;kk<700;kk++){
          value[kk] = 0.0;
          col[kk] = 0;
      }
      tcountcol2=tcountcol2 + tcountcol[ll]+1;
      countcol=tcountcol[ll]+1;
      if(i%getdata.npar == 0 && mpiid==0){
        ierr = PetscTime(&t1);CHKERRQ(ierr);
      }
      iii2=i+ll;
    ierr = MatSetValues(A,1,&iii2,countcol,col,value,INSERT_VALUES);CHKERRQ(ierr);
      if(i%getdata.npar == 0 && mpiid==0){
        ierr = PetscTime(&t2);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD," processor \ntime: %f\n",t2-t1);CHKERRQ(ierr);
      }
      }
  }
  ierr = PetscTime(&tt2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Time used to build the matrix: %f\n",tt2-tt1);CHKERRQ(ierr);


  ierr = PetscTime(&tt1);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = PetscTime(&tt2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Time used to assemble the matrix: %f\n",tt2-tt1);CHKERRQ(ierr);
//ierr = MatGetVecs(A,NULL,&xr);CHKERRQ(ierr);
//ierr = MatGetVecs(A,NULL,&xi);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&xr);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&xi);CHKERRQ(ierr);

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);

  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
  tol = 1.e-9;
  maxit = 10000000;
  ierr = EPSSetTolerances(eps,tol,maxit);CHKERRQ(ierr);
  nev  = 4;
  ncv  = 10;
  mpd  = 10;
  nev  = getdata.nroots;
  ierr = EPSSetDimensions(eps,nev,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);

  ierr = PetscTime(&t1);CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = PetscTime(&t2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Time used: %f\n",t2-t1);CHKERRQ(ierr);
  ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
//ierr = EPSPrintSolution(eps,NULL);CHKERRQ(ierr);

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n");CHKERRQ(ierr);

    for (i=0;i<nev;i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
      /*
         Compute the relative error associated to each eigenpair
      */
      ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif
      if (im!=0.0) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," %14f%+14fi %12g\n",(double)re,(double)im,(double)error);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   %18f       %12g\n",(double)re,(double)error);CHKERRQ(ierr);
      }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  }

  /*
     Save eigenvectors, if requested
  */
//PetscOptionsGetString(NULL,NULL,"-evecs",filename,PETSC_MAX_PATH_LEN,&evecs);
  EPSGetConverged(eps,&nconv);
  if (nconv>0) {
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
  	PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  	PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_SYMMODU);
    EPSIsHermitian(eps,&ishermitian);
    for (i=0;i<nconv;i++) {
      EPSGetEigenvector(eps,i,xr,xi);
      VecView(xr,viewer);
#if !defined(PETSC_USE_COMPLEX)
      if (!ishermitian) { VecView(xi,viewer); }
#endif
    }
    PetscViewerDestroy(&viewer);
  }

  /*
   * now analyzing the eigenvector
   */

  if (nconv>0) {

      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      ierr = EPSGetEigenpair(eps,1,&kr,&ki,xr,xi);CHKERRQ(ierr);

	  ierr = VecNorm(xr, NORM_2, &norm);CHKERRQ(ierr);
	  PetscPrintf(PETSC_COMM_WORLD," Norm = %18f \n", (double)norm);
  	  for (i=Istart; i<Iend; i+=1) {
		  indxr[i-Istart] = i;
	  }
	  ierr = VecGetValues(xr, nlocal, indxr, valxr);CHKERRQ(ierr);
  	  for (ii=Istart; ii<Iend; ii+=1) {
//  PetscPrintf(PETSC_COMM_WORLD," Element # = %d Value = %18f \n", i, valxr[i-Istart]); 
		  printf("\n Element # = %ld Value = %18f \n", ii, valxr[ii-Istart]); 
		  iii = ii+1;
		  getdet_(&iii, ideter);
		  for(kk=0; kk<8; kk++){
			  printf("%d ", ideter[kk]);
		  }
	  }

    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  }

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&xr);CHKERRQ(ierr);
  ierr = VecDestroy(&xi);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return 0;
}
