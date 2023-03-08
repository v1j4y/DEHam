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
  double sq2;
  int sumMs=0;
  int ntrouGauch = 0;
  int ntrouDroit = 0;
  int idx;

  // Create a new dictionary
  struct Dictionary *d = dictionary_new();
  struct Dictionary *dadr = dictionary_new();
  
  idx = prepare_dictionary(d, sze, MS2);
  idx = prepare_dictionary_for_adressing(dadr, sze, MS2);

  int dim;
  dim = (int)pow(3,nholes)*idx;
  double projvecmat[dim];
  for(ii=0;ii<dim;++ii) projvecmat[ii] = 0.0;
  printf(" Dimension = %d\n",dim);
  double msloc;
  //double normvec[dim];
  double fachole[27];
  int ms[nholes];
  int idh[nholes];

  sq2 = sqrt(2.0);

  fachole[0] = 1.0/8.0;
  fachole[1] = sq2/8.0;
  fachole[2] = 1.0/8.0;
  fachole[3] = sq2/8.0;
  fachole[4] = 2.0/8.0;
  fachole[5] = sq2/8.0;
  fachole[6] = 1.0/8.0;
  fachole[7] = sq2/8.0;
  fachole[8] = 1.0/8.0;

  fachole[ 9] = sq2/8.0;
  fachole[10] = 2.0/8.0;
  fachole[11] = sq2/8.0;
  fachole[12] = 2.0/8.0;
  fachole[13] = 2.0*sq2/8.0;
  fachole[14] = 2.0/8.0;
  fachole[15] = sq2/8.0;
  fachole[16] = 2.0/8.0;
  fachole[17] = sq2/8.0;

  fachole[18] = 1.0/8.0;
  fachole[19] = sq2/8.0;
  fachole[20] = 1.0/8.0;
  fachole[21] = sq2/8.0;
  fachole[22] = 2.0/8.0;
  fachole[23] = sq2/8.0;
  fachole[24] = 1.0/8.0;
  fachole[25] = sq2/8.0;
  fachole[26] = 1.0/8.0;

  for(ii=*Istart;ii<*Iend;ii++) {
    iii = ii + 1;
    iiii = ii;
    getdet_(&iii, ideter);
    
    // Set ms to 0
    for(kk=0;kk<nholes;++kk)
      ms[kk]=0;
   
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

    //// Prepare address in base 6
    //int addbase10,addbase6;
    //addbase6 = 0;
    //addbase10 = 0;
    //for(kk=0;kk<nholes;++kk){
    //  addbase10 += ms[kk]*(int)pow(6,kk);
    //}
    //addbase6 = base10ToBase6(addbase10);
    ////for(kk=0;kk<natom;++kk)
    ////  printf(" %d ",ideter[kk]);
    //printf("\n");
    //for(kk=0;kk<nholes;++kk)
    //  printf(" %d ",ms[kk]);
    //printf("\n");
    //printf("idx=%d addbase6 = %d nstate=%d \n",charPtrToInt(dictionary_get(dadr,intToCharPtr(addbase6))),addbase6, charPtrToInt(dictionary_get(d, intToCharPtr(addbase6))));

    ////int idhole = idh[0] * 1 + idh[1] * 3 + idh[2] * 9;
    //int idspin = charPtrToInt(dictionary_get(dadr,intToCharPtr(addbase6)));
    ////double fhole  = fachole[idhole];
    //double fspin  = 1.0/charPtrToInt(dictionary_get(d, intToCharPtr(addbase6)));
    ////projvecmat[idhole + idspin] += valxr[iiii]*fhole*fspin;
    ////projvecmat[idh[2]*3*3*idx + idh[1]*3*idx + idh[0]*idx + idspin] += valxr[iiii]*fhole*fspin;
    ////projvecmat[idhole * (int)(pow(3,nholes)) + idspin] += valxr[iiii]*fhole*fspin;

  }

  // Free dictionary
  dictionary_free(d);
  dictionary_free(dadr);

} /** END **/

