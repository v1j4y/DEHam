#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "get_dmat.h"


/*
 * 
 *----------------------------------------
 * Calculate the projection for 6_2h
 *----------------------------------------
 * Input
 * =====
 * valxr    = The full vector
 * Istart   = Local starting id
 * Iend     = Local ending id
 * Output
 * =====
 * projvec
 */
void get_proj_9_3h(PetscScalar *valxr, PetscInt *Istart, PetscInt *Iend, int *natom, int iroot, double *projvec, const int natomax){

  int			 ideter[natomax];
  int			 ideter1[natomax];
  int 		   	 kko,kok,kkio;
  int 		   	 mmo,mom,mmio;
  int 		   	 p, q, r;
  long int       ii;
  PetscInt      iiii;
  long int       iii;
  long int       iaa2, iaa;
  long int       nrow=-1, ncol=-1;
  double ms1=0.0, ms2=0.0 ,rest=0.0, norm=0.0;
  double sq2;
  double normvec[21];
  int sumMs=0;
  int ntrouGauch = 0;
  int ntrouMilieu = 0;
  int ntrouDroit = 0;
  double projvecmat[3*3*3*21];
  for(ii=0;ii<3*3*3*21;++ii) projvecmat[ii] = 0.0;
  sq2 = sqrt(2.0);

          	   for(ii=*Istart;ii<*Iend;ii++) {
                       iii = ii + 1;
                       iiii = ii;
                       getdet_(&iii, ideter);

                       ms1 = 0.0;
                       ms2 = 0.0;
                       p = -1;
                       q = -1;
                       // Find if Ms=5/2
                       sumMs = 0;
                       ntrouGauch = 0;
                       ntrouMilieu = 0;
                       ntrouDroit = 0;
                       for(kko=0;kko<3;++kko){
                         if(ideter[kko]==3){
                           ntrouGauch += 1;
                           p = kko;
                         }
                       }
                       for(kko=3;kko<6;++kko){
                         if(ideter[kko]==3){
                           ntrouMilieu += 1;
                           q = kko - 3;
                         }
                       }
                       for(kko=6;kko<9;++kko){
                         if(ideter[kko]==3){
                           ntrouDroit += 1;
                           r = kko - 3;
                         }
                       }
                       if(ntrouGauch == 1 && ntrouMilieu == 1 && ntrouDroit == 1){
                         sumMs = 1;
                       }
                       else{
                         sumMs = 0;
                       }
                       // First box ms
                       for(kko=0;kko<=2;++kko){
                         if(sumMs==1){
                         if(ideter[kko]==1){
                           ms1 = ms1 + 0.5;
                         }
                         else if(ideter[kko]==2){
                           ms1 = ms1 - 0.5;
                         }
                         }
                       }
                       for(kok= 11;kok>=9;--kok){
                         if(sumMs==1){
                         if(ideter[kok]==1){
                           ms1 = ms1 + 0.5;
                         }
                         else if(ideter[kok]==2){
                           ms1 = ms1 - 0.5;
                         }
                         }
                       }
                       // Second box ms
                       for(kko=0;kko<=2;++kko){
                         if(sumMs==1){
                         if(ideter[kko]==1){
                           ms2 = ms2 + 0.5;
                         }
                         else if(ideter[kko]==2){
                           ms2 = ms2 - 0.5;
                         }
                         }
                       }
                       for(kok= 11;kok>=9;--kok){
                         if(sumMs==1){
                         if(ideter[kok]==1){
                           ms2 = ms2 + 0.5;
                         }
                         else if(ideter[kok]==2){
                           ms2 = ms2 - 0.5;
                         }
                         }
                       }
                       if(ideter[1] ==3 && ideter[4] == 3){
                         norm += valxr[iiii]*valxr[iiii];
                       }

                       if(fabs(ms-2.5) < 1e-10){
                         projvecmat[6*3*p + 6*q + 0] = valxr[iiii];
                         normvec[0] += 1.0;
                       }
                       else if(fabs(ms+2.5) < 1e-10){
                         projvecmat[6*3*p + 6*q + 1] = valxr[iiii];
                         normvec[1] += 1.0;
                       }
                       else if(fabs(ms-1.5) < 1e-10){
                         projvecmat[6*3*p + 6*q + 2] += valxr[iiii];
                         normvec[2] += 1.0;
                       }
                       else if(fabs(ms+1.5) < 1e-10){
                         projvecmat[6*3*p + 6*q + 3] += valxr[iiii];
                         normvec[3] += 1.0;
                       }
                       else if(fabs(ms-0.5) < 1e-10){
                         projvecmat[6*3*p + 6*q + 4] += valxr[iiii];
                         normvec[4] += 1.0;
                       }
                       else if(fabs(ms+0.5) < 1e-10){
                         projvecmat[6*3*p + 6*q + 5] += valxr[iiii];
                         normvec[5] += 1.0;
                       }
                       else{
                         rest+=valxr[iiii];
                       }

                       ms = 0.0;

                   }

               for(p=0;p<3;++p){
                 for(q=0;q<3;++q) {
                   projvecmat[p*3*6 + q*6 + 0] = projvecmat[p*3*6 + q*6 + 0]/sqrt(1.0);
                   projvecmat[p*3*6 + q*6 + 1] = projvecmat[p*3*6 + q*6 + 1]/sqrt(1.0);
                   projvecmat[p*3*6 + q*6 + 2] = projvecmat[p*3*6 + q*6 + 2]/sqrt(25.0);
                   projvecmat[p*3*6 + q*6 + 3] = projvecmat[p*3*6 + q*6 + 3]/sqrt(25.0);
                   projvecmat[p*3*6 + q*6 + 4] = projvecmat[p*3*6 + q*6 + 4]/sqrt(100.0);
                   projvecmat[p*3*6 + q*6 + 5] = projvecmat[p*3*6 + q*6 + 5]/sqrt(100.0);
                 }
               }
                   //projvec[(iroot)*6 + 0] = projvec[(iroot)*6 + 0]/sqrt(1.0);
                   //projvec[(iroot)*6 + 1] = projvec[(iroot)*6 + 1]/sqrt(1.0);
                   //projvec[(iroot)*6 + 2] = projvec[(iroot)*6 + 2]/sqrt(25.0);
                   //projvec[(iroot)*6 + 3] = projvec[(iroot)*6 + 3]/sqrt(25.0);
                   //projvec[(iroot)*6 + 4] = projvec[(iroot)*6 + 4]/sqrt(100.0);
                   //projvec[(iroot)*6 + 5] = projvec[(iroot)*6 + 5]/sqrt(100.0);
               //for(ii=0;ii<6;++ii){
               //  projvec[(iroot)*6 + ii] =    \
               //                               projvecmat[1*3*6 + 1*6 + ii];
               //}
               for(ii=0;ii<6;++ii){
                 projvecmat[0*3*6 + 0*6 + ii] *= 1.0/4.0; // (1,1) = 1.0/4.0
                 projvecmat[0*3*6 + 1*6 + ii] *= sq2/4.0; // (1,2) = 1.4/4.0
                 projvecmat[0*3*6 + 2*6 + ii] *= 1.0/4.0; // (1,3) = 1.0/4.0

                 projvecmat[1*3*6 + 0*6 + ii] *= sq2/4.0; // (2,1) = 1.4/4.0
                 projvecmat[1*3*6 + 1*6 + ii] *= 2.0/4.0; // (2,2) = 2.0/4.0
                 projvecmat[1*3*6 + 2*6 + ii] *= sq2/4.0; // (2,3) = 1.4/4.0

                 projvecmat[2*3*6 + 0*6 + ii] *= 1.0/4.0; // (3,1) = 1.0/4.0
                 projvecmat[2*3*6 + 1*6 + ii] *= sq2/4.0; // (3,2) = 1.4/4.0
                 projvecmat[2*3*6 + 2*6 + ii] *= 1.0/4.0; // (3,3) = 1.0/4.0
               }

               for(ii=0;ii<6;++ii){
                 projvec[(iroot)*6 + ii] =    \
                                              projvecmat[0*3*6 + 0*6 + ii] \  
                                            + projvecmat[0*3*6 + 1*6 + ii] \ 
                                            + projvecmat[0*3*6 + 2*6 + ii] \ 
                                            + projvecmat[1*3*6 + 0*6 + ii] \ 
                                            + projvecmat[1*3*6 + 1*6 + ii] \ 
                                            + projvecmat[1*3*6 + 2*6 + ii] \ 
                                            + projvecmat[2*3*6 + 0*6 + ii] \ 
                                            + projvecmat[2*3*6 + 1*6 + ii] \ 
                                            + projvecmat[2*3*6 + 2*6 + ii];
               }
               //projvec[(iroot)*6 + 0] = projvec[(iroot)*6 *6+ 0]*6/sqrt(9.0);
               //projvec[(iroot)*6 + 1] = projvec[(iroot)*6 + 1]/sqrt(9.0);
               //projvec[(iroot)*6 + 2] = projvec[(iroot)*6 + 2]/sqrt(25.0*9);
               //projvec[(iroot)*6 + 3] = projvec[(iroot)*6 + 3]/sqrt(25.0*9);
               //projvec[(iroot)*6 + 4] = projvec[(iroot)*6 + 4]/sqrt(100.0*9);
               //projvec[(iroot)*6 + 5] = projvec[(iroot)*6 + 5]/sqrt(100.0*9);
               //printf(" norm = %4.4f projnorm = %4.4f\n",norm,projvec[5]*projvec[5] +projvec[4]*projvec[4] +projvec[3]*projvec[3] +projvec[2]*projvec[2] +projvec[1]*projvec[1] +projvec[0]*projvec[0]);
               //printf(" rest=%4.4f norm1=%2.2f norm2=%2.2f norm3=%2.2f norm4=%2.2f norm5=%2.2f norm6=%2.2f\n",rest,normvec[0],normvec[1],normvec[2],normvec[3],normvec[4],normvec[5]);
               //printf(" proj ( 5/2,-5/2) : 1,4 %2.3f 1,5 %2.3f 1,6 %2.3f 2,4 %2.3f 2,5 %2.3f 2,6 %2.3f 3,4 %2.3f 3,5 %2.3f 3,6 %2.3f \n",projvecseparate[0], projvecseparate[1], projvecseparate[2], projvecseparate[3], projvecseparate[4], projvecseparate[5], projvecseparate[6], projvecseparate[7], projvecseparate[8]);
               //printf(" proj (-5/2, 5/2) : 1,4 %2.3f 1,5 %2.3f 1,6 %2.3f 2,4 %2.3f 2,5 %2.3f 2,6 %2.3f 3,4 %2.3f 3,5 %2.3f 3,6 %2.3f \n",projvecseparate[0+9], projvecseparate[1+9], projvecseparate[2+9], projvecseparate[3+9], projvecseparate[4+9], projvecseparate[5+9], projvecseparate[6+9], projvecseparate[7+9], projvecseparate[8+9]);
               //printf(" proj ( 3/2,-3/2) : 1,4 %2.3f 1,5 %2.3f 1,6 %2.3f 2,4 %2.3f 2,5 %2.3f 2,6 %2.3f 3,4 %2.3f 3,5 %2.3f 3,6 %2.3f \n",projvecseparate[0+18], projvecseparate[1+18], projvecseparate[2+18], projvecseparate[3+18], projvecseparate[4+18], projvecseparate[5+18], projvecseparate[6+18], projvecseparate[7+18], projvecseparate[8+18]);
               //printf(" proj (-3/2, 3/2) : 1,4 %2.3f 1,5 %2.3f 1,6 %2.3f 2,4 %2.3f 2,5 %2.3f 2,6 %2.3f 3,4 %2.3f 3,5 %2.3f 3,6 %2.3f \n",projvecseparate[0+27], projvecseparate[1+27], projvecseparate[2+27], projvecseparate[3+27], projvecseparate[4+27], projvecseparate[5+27], projvecseparate[6+27], projvecseparate[7+27], projvecseparate[8+27]);
               //printf(" proj ( 1/2,-1/2) : 1,4 %2.3f 1,5 %2.3f 1,6 %2.3f 2,4 %2.3f 2,5 %2.3f 2,6 %2.3f 3,4 %2.3f 3,5 %2.3f 3,6 %2.3f \n",projvecseparate[0+36], projvecseparate[1+36], projvecseparate[2+36], projvecseparate[3+36], projvecseparate[4+36], projvecseparate[5+36], projvecseparate[6+36], projvecseparate[7+36], projvecseparate[8+36]);
               //printf(" proj (-1/2, 1/2) : 1,4 %2.3f 1,5 %2.3f 1,6 %2.3f 2,4 %2.3f 2,5 %2.3f 2,6 %2.3f 3,4 %2.3f 3,5 %2.3f 3,6 %2.3f \n",projvecseparate[0+45], projvecseparate[1+45], projvecseparate[2+45], projvecseparate[3+45], projvecseparate[4+45], projvecseparate[5+45], projvecseparate[6+45], projvecseparate[7+45], projvecseparate[8+45]);
} /** END **/

