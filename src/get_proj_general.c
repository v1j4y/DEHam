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
                      const int natomax,
                      int sze,
                      int MS2,
                      int nholes){

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
  int idx;
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
      //projvecmat[idhole + idspin] += valxr[iiii]*fhole*fspin;
      //projvecmat[idh[2]*3*3*idx + idh[1]*3*idx + idh[0]*idx + idspin] += valxr[iiii]*fhole*fspin;
      //printf("ii=%d val=%10.15f idx=%d\n",ii,valxr[iiii],idhole*idx + idspin);
      //projvecmat[idhole * idx + idspin] += valxr[iiii]*fhole*sqrt(fspin);
      //printf(" %3.5f %3.5f (%3.5f): %s | %3.5f | %d %d\n",fhole,sqrt(fspin),fhole*sqrt(fspin)*fhole*sqrt(fspin),intToCharPtr(addbase6),valxr[iiii],idspin,idhole);
      projvec[(iroot)*idx + idspin] += valxr[iiii]*fhole*sqrt(fspin);
      //if(addbase10==30){
      //  //normproj += valxr[iiii]*fhole;
      //  normproj += projvecmat[idhole * idx + idspin];
      //  for(kko=0;kko<6;++kko){
      //    printf(" %d ",ideter[kko]);
      //  }
      //  printf("\n");
      //  for(kko=11;kko>5;--kko){
      //    printf(" %d ",ideter[kko]);
      //  }
      //  printf("\n");
      //}
      //if(abs(valxr[iiii]) > 0.05) printf("- -- %3.5f\n",valxr[iiii]);
    }
  }
  //printf("normproj=%3.5f\n",normproj);

  // Sum over all hole positions
  //for(kko=0;kko<idx;++kko){
  //  for(ii=0;ii<nholes;++ii){
  //    idh[0]=ii;
  //    for(iii=0;iii<nholes;++iii){
  //      idh[1]=iii;
  //      //for(iiii=0;iiii<nholes;++iiii){
  //        //idh[2]=iiii;
  //        //int idhole = idh[0] * 1 + idh[1] * 3 + idh[2] * 9;
  //        int idhole = idh[0] * 1 + idh[1] * 3;
  //        //printf("%d %d %d idhole=%d val=%10.15f\n",idh[0],idh[1],idh[2],idhole,projvecmat[idhole*idx + kko]);
  //        //projvec[(iroot)*idx + kko] += projvecmat[idhole*idx + kko]; 
  //        if(kko==0){
  //          normproj += projvecmat[idhole * idx + kko];
  //        }
  //      //}
  //    }
  //  }
  //}
  //int id=0;
  //int p3=(int)pow(3,nholes);
  //for(kko=0;kko<idx;++kko){
  //  for(ii=0;ii<(int)pow(3,nholes);++ii){
  //    iii=ii;
  //    id=0;
  //    p3=(int)pow(3,nholes);
  //    for(kk=0;kk<3;++kk){
  //        id += p3 * (iii%3);
  //        p3 /= 3;
  //        iii=iii/3;
  //    }
  //    int idhole = id;
  //    projvec[(iroot)*idx + kko] += projvecmat[idhole*idx + kko]; 
  //  }
  //}

  //for(kko=0;kko<idx;++kko){
  //  printf("%d) %10.15f\n",kko,projvec[(iroot)*idx + kko]);
  //}

  // Free dictionary
  dictionary_free(d);
  dictionary_free(dadr);

} /** END **/

