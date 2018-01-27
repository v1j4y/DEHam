#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "get_dmat.h"


/*
 *---------------------------------------
 * One particle density matrix
 *----------------------------------------
 *
 * The One particle density matrix
 * Input
 * =====
 * valxr    = The full vector
 * Istart   = Local starting id
 * Iend     = Local ending id
 * natom    = number of sites
 * Output
 * =====
 * trace    = trace
 */
void get_1rdm(PetscScalar *valxr, PetscInt *Istart, PetscInt *Iend, int *natom, PetscReal *trace1rdm){

  const int            natomax=700;
  long int			 ideter[natomax];
  long int			 ideter2[natomax];
  int 		   	 kko,kok,kkio;
  long int       ii;
  PetscInt      iiii;
  long int       iii;
  long int       iaa2, iaa;
  int 		   	 ndim=(*natom)*(*natom)/8-(*natom)/2;
  double            densmat[ndim][ndim];
  memset(densmat, 0, sizeof(densmat[0][0]) * ndim * ndim);

            for(kko=0;kko<(*natom/2);kko++){
              for(kok=0;kok<(*natom/2);kok++){

          	  for(ii=*Istart;ii<*Iend;ii++) {
                      iii = ii + 1;
                      iiii = ii;
                      getdet_(&iii, ideter);
                      for(kkio=0;kkio<=*natom-1;kkio++){
                        ideter2[kkio]=ideter[kkio];
                      }
                      if(ideter[kko] != 3){
                        if(ideter[kok] == 3){
                          ideter2[kok]=ideter[kko];
                          ideter2[kko]=3;
                          adr_(ideter2, &iaa2);
                          densmat[kko][kok]=densmat[kko][kok]+valxr[iiii]*valxr[iaa2];
                        }
                      }
                      if(kko == kok && ideter[kko] != 3){
                        densmat[kko][kko]=densmat[kko][kko]+valxr[iiii]*valxr[iiii];
                      }

                  }
                  if(kko == kok){
                      *trace1rdm+=densmat[kko][kko];
                  }

              }
            }
} /** END **/

/*
 * 
 *----------------------------------------
 * two particle density matrix
 *----------------------------------------
 * Input
 * =====
 * valxr    = The full vector
 * Istart   = Local starting id
 * Iend     = Local ending id
 * Output
 * =====
 * trace    = trace
 */
void get_2rdm(PetscScalar *valxr, PetscInt *Istart, PetscInt *Iend, int *natom, PetscReal *trace2rdm, double densmat2[*natom][*natom][*natom][*natom]){

  const int            natomax=700;
  long int			 ideter[natomax];
  long int			 ideter2[natomax];
  int 		   	 kko,kok,kkio;
  int 		   	 mmo,mom,mmio;
  long int       ii;
  PetscInt      iiii;
  long int       iii;
  long int       iaa2, iaa;
  long int       nrow=-1, ncol=-1;
//int 		   	 ndim=(*natom/2)*((*natom/2)-1)/2;
//double            densmat2[ndim][ndim];
//memset(densmat2, 0, sizeof(densmat2[0][0]) * ndim * ndim);

            for(kko=0;kko<(*natom/2);kko++){
              for(kok=0;kok<(*natom/2);kok++){

               nrow=nrow+1;
               ncol=-1;
               for(mmo=0;mmo<(*natom/2);mmo++){
                 for(mom=0;mom<(*natom/2);mom++){

                   ncol=ncol+1;

          	   for(ii=*Istart;ii<*Iend;ii++) {
                       iii = ii + 1;
                       iiii = ii;
                       getdet_(&iii, ideter);
                       for(kkio=0;kkio<=*natom-1;kkio++){
                         ideter2[kkio]=ideter[kkio];
                       }
                       if(ideter[kko] == 3 && ideter[kok] == 3 && kko != kok && mmo != mom){
                          ideter2[kko]=ideter[mmo];
                          ideter2[mmo]=3;
                          ideter2[kok]=ideter[mom];
                          ideter2[mom]=3;
                          adr_(ideter2, &iaa2);
                          densmat2[kko][kok][mmo][mom]=densmat2[kko][kok][mmo][mom]+valxr[iiii]*valxr[iaa2];
                       }


                     if(kko == mmo && kok == mom && ideter[kko]==3 && ideter[kok]==3 && kko != kok){
                       densmat2[kko][kok][mmo][mom]=densmat2[kko][kok][mmo][mom]+valxr[iiii]*valxr[iiii];
                     }

                   }
                          printf("%d\t%d\t%d\t%d\t%18f\n",kko,kok,mmo,mom,densmat2[kko][kok][mmo][mom]);


                    if(kko == mmo && kok == mom)*trace2rdm+=densmat2[kko][kok][mmo][mom];
                 }
               }


              }
            }

} /** END **/

