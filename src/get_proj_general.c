#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "get_proj_general.h"


/*
 * 
 *-------------------------------------------
 * Calculate the projection for general size
 *-------------------------------------------
 * Input
 * =====
 * valxr    = The full vector
 * Istart   = Local starting id
 * Iend     = Local ending id
 * Output
 * =====
 * projvec
 */
void get_proj_general(PetscScalar *valxr, 
                      PetscInt *Istart, 
                      PetscInt *Iend, 
                      int *natom, 
                      int iroot, 
                      double *projvec, 
                      double *projvec2, 
                      const int natomax,
                      int sze,
                      int MS2,
                      int nholes,
                      double XS,
                      int colm,
                      double *projmatrix){

  int			 ideter[natomax];
  int			 ideter2[natomax];
  int 		   	 kk,kko,kok,kkio;
  int 		   	 mmo,mom,mmio;
  int 		   	 p, q;
  long int       ii;
  PetscInt      iiii;
  long int       iii;
  long int       iaa2, iaa;
  long int       nrow=-1, ncol=-1;
  double rest=0.0, norm=0.0;
  double normproj=0.0;
  int sumMs=0;
  int ntrouGauch = 0;
  int ntrouDroit = 0;
  int idx, idxprojm;
  int checkIonic;
  int pow3 = (int)pow(3,nholes);

  // Create a new dictionary
  struct Dictionary *d = dictionary_new();
  struct Dictionary *dadr = dictionary_new();
  
  idx = prepare_dictionary(d, sze, MS2);
  idx = prepare_dictionary_for_adressing(dadr, sze, MS2);

  int dim;
  dim = pow3*idx;
  double projvecmat[dim];
  for(ii=0;ii<dim;++ii) projvecmat[ii] = 0.0;
  for(ii=0;ii<idx;++ii) projvec[(iroot)*idx + ii] = 0.0;
  double msloc;
  //double normvec[dim];
  int ms[nholes];
  int idh[nholes];

  if(fabs(XS-0.0) < 1E-10){
    idxprojm = 0;
  }
  else if(fabs(XS-0.5) < 1E-10){
    idxprojm = 0;
  }
  else if(fabs(XS-1.0) < 1E-10){
    idxprojm = 1;
  }
  else if(fabs(XS-1.5) < 1E-10){
    idxprojm = 1;
  }
  else if(fabs(XS-2.0) < 1E-10){
    idxprojm = 2;
  }
  else if(fabs(XS-2.5) < 1E-10){
    idxprojm = 2;
  }
  else if(fabs(XS-3.0) < 1E-10){
    idxprojm = 3;
  }
  else if(fabs(XS-3.5) < 1E-10){
    idxprojm = 3;
  }
  else if(fabs(XS-4.0) < 1E-10){
    idxprojm = 4;
  }
  else if(fabs(XS-4.5) < 1E-10){
    idxprojm = 4;
  }
  else if(fabs(XS-5.0) < 1E-10){
    idxprojm = 5;
  }
  else if(fabs(XS-5.5) < 1E-10){
    idxprojm = 5;
  }
  else if(fabs(XS-6.0) < 1E-10){
    idxprojm = 6;
  }
  else if(fabs(XS-6.5) < 1E-10){
    idxprojm = 6;
  }
  else if(fabs(XS-7.0) < 1E-10){
    idxprojm = 7;
  }
  else if(fabs(XS-7.5) < 1E-10){
    idxprojm = 7;
  }
  else if(fabs(XS-8.0) < 1E-10){
    idxprojm = 8;
  }
  else if(fabs(XS-8.5) < 1E-10){
    idxprojm = 8;
  }

  int sizefac = (int)pow(3,nholes);
  double fachole[sizefac];
  //printf(" idx=%d dim=%d sizefac=%d sze=%d MS2=%d\n",idx,dim,sizefac,sze,MS2);
  prepareHueckelFactors(nholes, fachole, sizefac);

  for(ii=*Istart;ii<*Iend;ii++) {
    iii = ii + 1;
    iiii = ii;
    getdet_(&iii, ideter);
    //normproj += valxr[iiii]*valxr[iiii];
    
    checkIonic = 0;
    // Set ms to 0
    for(kk=0;kk<nholes;++kk)
      ms[kk]=0;
    
    // Set holes to -1
    for(kk=0;kk<nholes;++kk)
      idh[kk]=-1;
   
    // Calculate MS
    for(kk=0;kk<nholes;++kk){
      // Find position of hole
      for(kko=kk*3+0;kko<(kk+1)*3;++kko){
        if(ideter[kko]==3){
          idh[kk] = kko - kk*3;
        }
      }
      msloc = 0;
      for(kko=kk*3+0;kko<(kk+1)*3;++kko){
        if(ideter[kko]==1){
          msloc += 1;
        }
      }
      for(kko=(nholes*3*2)-1-kk*3;kko>(nholes*3*2)-1-(kk+1)*3;--kko){
        if(ideter[kko]==1){
          msloc += 1;
        }
      }
      ms[kk] = msloc;
    }
    
    for(kk=0;kk<nholes;++kk){
      if(idh[kk] == -1){
        checkIonic = 1;
      }
    }


    if(checkIonic == 0) {
      //if(ms[0]==5 || ms[1]==5){
      //  normproj += valxr[iiii]/sqrt(9);
      //}

      // Prepare address in base 6
      int addbase10,addbase6;
      addbase6 = 0;
      addbase10 = 0;
      for(kk=0;kk<nholes;++kk){
        addbase10 += ms[kk]*(int)pow(6,kk);
      }
      addbase6 = base10ToBase6(addbase10);

      //int idhole = idh[0] * 1 + idh[1] * 3 + idh[2] * 9;
      int idhole=0;
      for(kk=0;kk<nholes;++kk){
        idhole += idh[kk]*((int)pow(3,kk));
      }
      int idspin = charPtrToInt(dictionary_get(dadr,intToCharPtr(addbase6)));
      double fhole  = fachole[idhole];
      double fspin  = 1.0/charPtrToInt(dictionary_get(d, intToCharPtr(addbase6)));
      projvec[(iroot)*idx + idspin] += valxr[iiii]*fhole*sqrt(fspin);
      projvec2[(iroot)*idx + idspin] += valxr[iiii]*fhole*sqrt(fspin)*projmatrix[idxprojm * colm + idspin];
    }
  }
  //printf("idxprojm=%d XS=%6.4f \n",idxprojm,XS);
  //for(int i=0;i<6;++i){
  //  printf(" (%8.4f %8.4f ) ",projvec[(iroot)*idx + i],projmatrix[idxprojm * colm + i]);
  //}
  //printf("\n");

  // Free dictionary
  dictionary_free(d);
  dictionary_free(dadr);

} /** END **/

