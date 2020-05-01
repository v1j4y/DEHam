#include <petsctime.h>
#include <petscvec.h>

#include "read2.h"
#include "stimsyr.h"
#include "get_s2.h"
#include "get_ntot.h"

#undef __FUNCT__
#define __FUNCT__ "main"

void solvequad(double *a, double *b, double *c, double *res){
    *res = -*b/(2.0*(*a)) + sqrt((*b)*(*b) - 4.0*(*a)*(*c))/(2.0*(*a));
}

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscReal      norm=0.0;
  PetscReal      norm2=0.0;
  PetscReal      norm3=0.0;
  PetscReal      norm4=0.0;
  PetscReal      normfin=0.0;
  PetscReal      normfin2=0.0;
  PetscReal      normfin3=0.0;
  PetscReal      normfin4=0.0;
  const int            natomax=900;
  PetscScalar    kr,ki,value[natomax];
  Vec            xr,xi;
  PetscInt       i,Istart,Iend,col[natomax],maxit,its,nconv,countcol;
  PetscInt		 nev, ncv, mpd;
  PetscLogDouble t1,t2,tt1,tt2;
  PetscErrorCode ierr;
  int            mpiid;
  SlepcInitialize(&argc,&argv,(char*)0,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpiid);

  char const* const fileName = argv[1];
  FILE* file = fopen(fileName, "r");
  Data getdata;
//PetscInt		 nlocal;
  
  /* gather the input data */
  if(mpiid==0)printf("Reading data\n");
  Data_new(file, &getdata);
  getdata.n = get_ntot(getdata.FAM1, getdata.natom, getdata.isz, getdata.ntrou, getdata.fix_trou1, getdata.fix_trou2);
  if(mpiid==0)printf("Done reading data\n");

//nlocal = getdata.n/getdata.npar;

  PetscScalar	 *valxr;
//PetscInt	   	 indxr[nlocal];
  char           filename[PETSC_MAX_PATH_LEN]="FIL666";
  PetscViewer    viewer;
  PetscBool      ishermitian;
  PetscInt       kk,ll,mm,nn,iii2,iiii;
  PetscInt       ii;
  long int       iii;
  long int       tcountcol2,tcol[natomax],tcountcol[getdata.nnz];
  double         val[natomax];
  PetscReal      xymat=0.0;
  PetscReal      xymat2=0.0;
  PetscReal      xymat3=0.0;
  PetscReal      xymat4=0.0;
  PetscReal      xymatfin=0.0;
  PetscReal      xymatfin2=0.0;
  PetscReal      xymatfin3=0.0;
  PetscReal      xymatfin4=0.0;
  PetscReal      weight3fin = 0.0;
  PetscReal      XS = 0.0;
  PetscReal      XS2 = 0.0;
  PetscReal      XS3 = 0.0;
  PetscReal      XS4 = 0.0;
  PetscReal      W3  = 0.0;
  PetscReal      weight3 = 0.0;
  PetscReal      trace1rdm=0.0;
  PetscReal      trace1rdmfin=0.0;
  PetscReal      trace2rdm=0.0;
  PetscReal      trace2rdmfin=0.0;
  IS                from, to; /* index sets that define the scatter */
//PetscInt          idx_to[nlocal], idx_from[nlocal];
  PetscScalar       *values;
  int 		   	 ndim=(getdata.natom/2)*((getdata.natom/2)-1)/2;
  double            a, b, c;
  double            gamma_p = 0.0, gamma_m = 0.0;
  double            gamma_pfin = 0.0, gamma_mfin = 0.0;
  double            nel, s2dens;
  double            nelfin, s2densfin;
//double            densmat2[getdata.natom][getdata.natom][getdata.natom][getdata.natom];
//memset(densmat2, 0, sizeof(densmat2));

  if(mpiid==0)printf("Initializing Slepc vars\n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1-D t-J Eigenproblem, n=%D\n\n",getdata.n);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,getdata.n,getdata.n,10*getdata.natom,NULL,10*getdata.natom,NULL,&A);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,10*getdata.natom,NULL,10*getdata.natom,NULL);CHKERRQ(ierr);
  if(mpiid==0)printf("Done Initializing Slepc\n");

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  ierr = PetscTime(&tt1);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," start: %d end: %d \n",Istart, Iend);CHKERRQ(ierr);

  for (i=Istart; i<Iend; i+=getdata.nnz) {
      tcountcol2=0;
      for(kk=0;kk<getdata.nnz;kk++){
          tcountcol[kk]=0;
      }
      iii=i+1;
      unit_l1_(
            getdata.l1,
            getdata.l2,
            getdata.ktyp,
            &iii,
            &getdata.nnz,
            getdata.xjjxy,
            getdata.xjjz ,
            getdata.xtt ,
            getdata.E ,
            tcountcol,
            &getdata.ntrou,
            &getdata.isz,
            &getdata.fix_trou1,
            &getdata.fix_trou2,
            &getdata.FAM1,
            &mpiid,
            tcol,
            val);
//    if(i%getdata.npar == 0 && mpiid==0){
//      ierr = PetscPrintf(PETSC_COMM_WORLD," i: %d \n",i);CHKERRQ(ierr);
//    }
      for(ll=0;ll<getdata.nnz;ll++){

//      printf("%d) ll=%d countcol=%d\n",i,ll,tcountcol[ll]+1);
        for(kk=0;kk<tcountcol[ll]+1;kk++){
            value[kk] = val[kk+tcountcol2];
            col[kk] = tcol[kk+tcountcol2]-1;
//          printf("%d) kk=%d col=%d val=%1.4f\n",i,kk,col[kk],value[kk]);
        }
        for(kk=tcountcol2+tcountcol[ll]+1;kk<natomax;kk++){
            value[kk] = 0.0;
            col[kk] = 0;
        }
        tcountcol2=tcountcol2 + tcountcol[ll]+1;
        countcol=tcountcol[ll]+1;
        iii2=i+ll;
        ierr = MatSetValues(A,1,&iii2,countcol,col,value,INSERT_VALUES);CHKERRQ(ierr);
      }
  }
  ierr = PetscTime(&tt2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Time used to build the matrix: %f\n",tt2-tt1);CHKERRQ(ierr);
  printf("time = %f mpiid = %d \n",tt2-tt1, mpiid);


  ierr = PetscTime(&tt1);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = PetscTime(&tt2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Time used to assemble the matrix: %f\n",tt2-tt1);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&xr);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&xi);CHKERRQ(ierr);

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);

  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
  tol = 1.e-9;
  maxit = 10000000;
  ierr = EPSSetTolerances(eps,tol,maxit);CHKERRQ(ierr);
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
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenvalues: %D\n",nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
//ierr = EPSPrintSolution(eps,NULL);CHKERRQ(ierr);

  /*
     Save eigenvectors, if  == ested
  */
  EPSGetConverged(eps,&nconv);
  if (getdata.print_wf) {
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
  	PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  	PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_SYMMODU);
    EPSIsHermitian(eps,&ishermitian);
    for (i=0;i<nev;i++) {
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

    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "               k          ||Ax-kx||/||kx||         <S>\n"
         "       ----------------- ----------------- ------------------\n");CHKERRQ(ierr);

   for(i=0;i<nev;i++){

      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      Vec               vec2;
      VecScatter        scatter; /* scatter context */
      ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
	  xymat = 0.0;
	  xymat2 = 0.0;
	  xymat3 = 0.0;
	  xymat4 = 0.0;
	  weight3 = 0.0;
	  norm = 0.0;
	  norm2 = 0.0;
	  norm3 = 0.0;
	  norm4 = 0.0;


//        ierr = PetscTime(&tt1);CHKERRQ(ierr);
//        ierr = VecGetArray(xr, &valxr);CHKERRQ(ierr);
          VecScatterCreateToAll(xr,&scatter,&vec2);
          VecScatterBegin(scatter,xr,vec2,INSERT_VALUES,SCATTER_FORWARD);
          VecScatterEnd(scatter,xr,vec2,INSERT_VALUES,SCATTER_FORWARD);
          ierr = VecGetArray(vec2,&values);CHKERRQ(ierr);
          get_s2(xr, &Istart, &Iend, values, &getdata.natom, &norm, &norm2, &norm3, &norm4, &xymat, &xymat2, &xymat3, &xymat4, &weight3,
                  &getdata.s21a1, &getdata.s21a2, &getdata.s21b1, &getdata.s21b2, &getdata.s22a1, &getdata.s22a2,
                  &getdata.s22b1, &getdata.s22b2,  &getdata.s23a1, &getdata.s23a2,
                  &getdata.s23b1, &getdata.s23b2, &getdata.postrou, natomax);
//        get_s2_cyclic(xr, &Istart, &Iend, values, &getdata.natom, &norm, &norm2, &norm3, &norm4, &xymat, &xymat2, &xymat3, &xymat4,
//                &getdata.s21a1, &getdata.s21a2, &getdata.s21b1, &getdata.s21b2, &getdata.s22a1, &getdata.s22a2,
//                &getdata.s22b1, &getdata.s22b2,  &getdata.s23a1, &getdata.s23a2,
//                &getdata.s23b1, &getdata.s23b2, &getdata.postrou, natomax);
//        get_1rdm(values, &Istart, &Iend, &getdata.natom, &trace1rdm, natomax);
//        get_2rdm(values, &Istart, &Iend, &getdata.natom, &trace2rdm, densmat2, natomax);
//        analyse_(valxr, (Iend-Istart), &Istart, &Iend, &xymat, &norm);
          VecRestoreArray(vec2,&values);
          ierr = VecRestoreArray(xr, &valxr);CHKERRQ(ierr);
          MPI_Reduce(&xymat, &xymatfin, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
          MPI_Reduce(&xymat2, &xymatfin2, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
          MPI_Reduce(&xymat3, &xymatfin3, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
          MPI_Reduce(&xymat4, &xymatfin4, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
          MPI_Reduce(&weight3, &weight3fin, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
          MPI_Reduce(&norm, &normfin, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
          MPI_Reduce(&norm2, &normfin2, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
          MPI_Reduce(&norm3, &normfin3, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
          MPI_Reduce(&norm4, &normfin4, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//        MPI_Reduce(&trace1rdm, &trace1rdmfin, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//          printf("done calc densmat\n");
//          for(ll=0;ll<getdata.natom/2;ll++){
//              for(kk=0;kk<getdata.natom/2;kk++){
//                  gamma_p = gamma_p + 0.5*(densmat2[ll][kk][kk][ll] + densmat2[ll][kk][ll][kk]);
//                  gamma_m = gamma_m + 0.5*(densmat2[ll][kk][kk][ll] - densmat2[ll][kk][ll][kk]);
//              }
//          }
//          MPI_Reduce(&trace2rdm, &trace2rdmfin, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//          MPI_Reduce(&gamma_p, &gamma_pfin, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//          MPI_Reduce(&gamma_m, &gamma_mfin, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//          if(mpiid==0){
//            for(kk=0;kk<getdata.natom;kk++){
//                for(ll=0;ll<getdata.natom;ll++){
//                    for(mm=0;mm<getdata.natom;mm++){
//                        for(nn=0;nn<getdata.natom;nn++){
////                          printf("%d\t%d\t%d\t%d\t%18f\n",kk,ll,mm,nn,densmat2[kk][ll][mm][nn]);
//                        }
//                    }
//                }
//            }
//            /* calc nel */
//            a=1.0;
//            b=-1.0;
//            c=-2.0*(gamma_mfin + gamma_pfin);
//            printf("\n  gp= %18f gm= %18f a=%18f b=%18f c=%18f\n", gamma_pfin, gamma_mfin, a, b, c);
//            nel = -b/(2.0*(a)) + sqrt((b)*(b) - 4.0*(a)*(c))/(2.0*(a));
////          solvequad(&a, &b, &c, &nel);
//
//            /* calc s^2 */
//            a=1.0;
//            b=1.0;
//            c=-1.0*((gamma_mfin - gamma_pfin) - nel*(nel - 4.0)/4.0);
//            s2dens = -b/(2.0*(a)) + sqrt((b)*(b) - 4.0*(a)*(c))/(2.0*(a));
////          solvequad(&a, &b, &c, &s2dens);
//            printf("\n mpiid = %d  # trace = %18f nel = %18f s2dens = %18f\n", mpiid, trace2rdmfin, nel, s2dens);
//          }

          if(!mpiid){
            XS=(1.0/2.0)*(-1.0+sqrt(1.0+(4.0*xymatfin/normfin)));
//          XS2=(1.0/2.0)*(-1.0+sqrt(1.0+(4.0*xymatfin2/normfin2)));
//          XS3=(1.0/2.0)*(-1.0+sqrt(1.0+(4.0*xymatfin3/normfin3)));
            XS2=(1.0/2.0)*(-1.0+sqrt(1.0+(4.0*xymatfin2)));
            XS3=(1.0/2.0)*(-1.0+sqrt(1.0+(4.0*xymatfin3)));
            XS4=(1.0/2.0)*(-1.0+sqrt(1.0+(4.0*xymatfin4/normfin4)));
            XS4=(1.0/2.0)*(-1.0+sqrt(1.0+(4.0*xymatfin4/normfin4)));
            W3=weight3fin/normfin2;
//          W3=weight3fin;
          }

  xymatfin = 0.0;
  normfin = 0.0;

      /*
       * Compute the relative error associated to each eigenpair
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
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   %18f     %12g %18f %18f %18f %18f\n",(double)re,(double)error,(double)XS,(double)XS2,(double)XS3, (double)W3);CHKERRQ(ierr);
      }
      VecScatterDestroy(&scatter);
      VecDestroy(&vec2); 
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
