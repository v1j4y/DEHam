#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "get_s2.h"
#include "get_val_iaa2.h"

/* 
 * This function simply calculates the S^2 value of the wavefunction
 * Input
 * =====
 * Vr       = The full vector
 * Istart   = Local starting id of the vector
 * Iend     = Local vector ending id
 * valxr    = Local vector values
 * natom    = number of orbitals
 * Output
 * ======
 * norm     = norm of the vector 
 * xymat    = the S^2 value
 */

void get_s2_mov(Vec xr, PetscInt *Istart, PetscInt *Iend, PetscScalar *valxr, int *natom, PetscReal *norm, PetscReal *norm2, PetscReal *norm3, PetscReal *norm4, PetscReal *xymat, PetscReal *xymat2, PetscReal *xymat3, PetscReal *xymat4, PetscReal *weight3,
            int *s21a1, int *s21a2, int *s21b1,  int *s21b2,  int *s22a1,  int *s22a2,  int *s22b1, int *s22b2, int *s23a1,  int *s23a2,  int *s23b1, int *s23b2, int *postrou){
  const int            natomax=700;
  long int       iaa2, iaa;
  long int       iii;
  long int			 ideter[natomax];
  long int			 ideter2[natomax];
  int 		   	 kko,kok,kkio;
  long int       ii;
  double         xmat=0.0;
  double         xmat2=0.0;
  double         xmat3=0.0;
  double         xmat4=0.0;
  double         getvaliaa2;
  PetscLogDouble t1,t2,tt1,tt2;
  PetscErrorCode ierr;
  PetscInt      iiii;
  int            ntrouboit1=0;
  int            ntrouboit2=0;
  int            ntrouboit3=0;
  int            okboit1=0;
  int            okboit2=0;
  int            okboit3=0;
  int            mpiid;
  int            pos1=0;
  int            pos2=0;
  int            pos3=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpiid);
//if(!mpiid){printf("istart= %d ind = %d\n",*Istart,*Iend);}
//ierr = PetscTime(&tt1);CHKERRQ(ierr);
  	  for(ii=*Istart;ii<*Iend;ii++) {
              iii = ii + 1;
//            iiii = ii-*Istart;
              iiii = ii;
              xmat = 0.0;
              xmat2 = 0.0;
              xmat3 = 0.0;
              xmat4 = 0.0;
              ntrouboit1 = 0;
              ntrouboit2 = 0;
              ntrouboit3 = 0;
              okboit1 = 0;
              okboit2 = 0;
              okboit3 = 0;
              pos1 = 0;
              pos2 = 0;
              pos3 = 0;
              getdet_(&iii, ideter);
              *norm=*norm+valxr[iiii]*valxr[iiii];
              for(kko=*s21a1;kko<=*s21a2;kko++){
                  if(ideter[kko]==3){
                      ntrouboit1++;
                      pos1=kko;
                  }
              }
              for(kko=*s22a1;kko<=*s22a2;kko++){
                  if(ideter[kko]==3){
                      ntrouboit2++;
                      pos2=kko;
                  }
              }
              for(kko=*s23a1;kko<=*s23a2;kko++){
                  if(ideter[kko]==3){
                      ntrouboit3++;
                      pos3=kko;
                  }
              }
              if(ntrouboit1==1 && *s21a1 <= pos1 && pos1 <= *s21a2)okboit1=1;
              if(ntrouboit2==1 && pos2 == *postrou)okboit2=1;
              if(ntrouboit3==1 && pos3 == *postrou)okboit3=1;
              if(okboit1){
                *norm2=*norm2+valxr[iiii]*valxr[iiii];
              }
              if(okboit2){
                *norm3=*norm3+valxr[iiii]*valxr[iiii];
              }
              if(okboit3){
                *norm4=*norm4+valxr[iiii]*valxr[iiii];
              }
              /*
               * calculate the weight of ms=5/2
               *
               * loop over the determinants to see if we have a S=5/2
               */
              int countw = 0;
              for(kko=*s21a1;kko<=*s21a2;kko++){
                  if(ideter[kko] == 2) countw=1;
              }
              for(kok=*s21b1;kok<=*s21b2;kok++){
                  if(ideter[kok] == 2) countw=1;
              }
              if(countw==0 && okboit1){
                  *weight3 += (valxr[iiii]*valxr[iiii]);
              }
              for(kko=0;kko<=(*natom/2)-1;kko++){
                  for(kok=kko;kok<=(*natom/2)-1;kok++){
                    if(kok == kko && ideter[kok] != 3){
                      xmat=xmat+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                      if(okboit1){
                        if( kko >=*s21a1 && kko <=*s21a2){
                          if( kok >=*s21a1 && kok <=*s21a2){
                                xmat2=xmat2+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22a1 && kko <=*s22a2){
                          if( kok >=*s22a1 && kok <=*s22a2){
                                xmat3=xmat3+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23a1 && kko <=*s23a2){
                          if( kok >=*s23a1 && kok <=*s23a2){
                                xmat4=xmat4+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                    }
                    else{
                      if(ideter[kko] == 1 && ideter[kok] == 1){
                        xmat=xmat+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                        if(okboit1){
                          if( kko >=*s21a1 && kko <=*s21a2){
                            if( kok >=*s21a1 && kok <=*s21a2){
                                  xmat2=xmat2+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                        if(okboit2){
                          if( kko >=*s22a1 && kko <=*s22a2){
                            if( kok >=*s22a1 && kok <=*s22a2){
                                  xmat3=xmat3+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                        if(okboit3){
                          if( kko >=*s23a1 && kko <=*s23a2){
                            if( kok >=*s23a1 && kok <=*s23a2){
                                  xmat4=xmat4+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                      }
                      if(ideter[kko] == 2 && ideter[kok] == 2){
                        xmat=xmat+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                        if(okboit1){
                          if( kko >=*s21a1 && kko <=*s21a2){
                            if( kok >=*s21a1 && kok <=*s21a2){
                                  xmat2=xmat2+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                        if(okboit2){
                          if( kko >=*s22a1 && kko <=*s22a2){
                            if( kok >=*s22a1 && kok <=*s22a2){
                                  xmat3=xmat3+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                        if(okboit3){
                          if( kko >=*s23a1 && kko <=*s23a2){
                            if( kok >=*s23a1 && kok <=*s23a2){
                                  xmat4=xmat4+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                      }
                      if(ideter[kko] == 1 && ideter[kok] == 2){
                        xmat=xmat-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                        if(okboit1){
                          if( kko >=*s21a1 && kko <=*s21a2){
                            if( kok >=*s21a1 && kok <=*s21a2){
                                  xmat2=xmat2-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                        if(okboit2){
                          if( kko >=*s22a1 && kko <=*s22a2){
                            if( kok >=*s22a1 && kok <=*s22a2){
                                  xmat3=xmat3-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                        if(okboit3){
                          if( kko >=*s23a1 && kko <=*s23a2){
                            if( kok >=*s23a1 && kok <=*s23a2){
                                  xmat4=xmat4-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                        for(kkio=0;kkio<=*natom-1;kkio++){
                          ideter2[kkio]=ideter[kkio];
                        }
                        ideter2[kko]=2;
                        ideter2[kok]=1;
                        adr_(ideter2, &iaa2);
                        iaa2 = iaa2 - 1;
                        xmat=xmat+valxr[iiii]*valxr[iaa2];
                        if(okboit1){
                          if( kko >=*s21a1 && kko <=*s21a2){
                            if( kok >=*s21a1 && kok <=*s21a2){
                                  xmat2=xmat2+(valxr[iiii]*valxr[iaa2]);
                            }
                          }
                        }
                        if(okboit2){
                          if( kko >=*s22a1 && kko <=*s22a2){
                            if( kok >=*s22a1 && kok <=*s22a2){
                                  xmat3=xmat3+(valxr[iiii]*valxr[iaa2]);
                            }
                          }
                        }
                        if(okboit3){
                          if( kko >=*s23a1 && kko <=*s23a2){
                            if( kok >=*s23a1 && kok <=*s23a2){
                                  xmat4=xmat4+(valxr[iiii]*valxr[iaa2]);
                            }
                          }
                        }
                      }
                      if(ideter[kko] == 2 && ideter[kok] == 1){
                        xmat=xmat-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                        if(okboit1){
                          if( kko >=*s21a1 && kko <=*s21a2){
                            if( kok >=*s21a1 && kok <=*s21a2){
                                  xmat2=xmat2-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                        if(okboit2){
                          if( kko >=*s22a1 && kko <=*s22a2){
                            if( kok >=*s22a1 && kok <=*s22a2){
                                  xmat3=xmat3-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                        if(okboit3){
                          if( kko >=*s23a1 && kko <=*s23a2){
                            if( kok >=*s23a1 && kok <=*s23a2){
                                  xmat4=xmat4-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                            }
                          }
                        }
                        for(kkio=0;kkio<=*natom-1;kkio++){
                          ideter2[kkio]=ideter[kkio];
                        }
                        ideter2[kko]=1;
                        ideter2[kok]=2;
                        adr_(ideter2, &iaa2);
                        iaa2 = iaa2 - 1;
                        xmat=xmat+valxr[iiii]*valxr[iaa2];
                        if(okboit1){
                          if( kko >=*s21a1 && kko <=*s21a2){
                            if( kok >=*s21a1 && kok <=*s21a2){
                                  xmat2=xmat2+(valxr[iiii]*valxr[iaa2]);
                            }
                          }
                        }
                        if(okboit2){
                          if( kko >=*s22a1 && kko <=*s22a2){
                            if( kok >=*s22a1 && kok <=*s22a2){
                                  xmat3=xmat3+(valxr[iiii]*valxr[iaa2]);
                            }
                          }
                        }
                        if(okboit3){
                          if( kko >=*s23a1 && kko <=*s23a2){
                            if( kok >=*s23a1 && kok <=*s23a2){
                                  xmat4=xmat4+(valxr[iiii]*valxr[iaa2]);
                            }
                          }
                        }
                       }
                    }
                }
              }
              for(kko=(*natom/2);kko<=*natom-1;kko++){
                  for(kok=kko;kok<=*natom-1;kok++){
                    if(kok == kko && ideter[kok] != 3){
                      xmat=xmat+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                      if(okboit1){
                        if( kko >=*s21b1 && kko <=*s21b2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22b1 && kko <=*s22b2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23b1 && kko <=*s23b2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                    }
                    else{
                      if(ideter[kko] == 1 && ideter[kok] == 1){
                        xmat=xmat+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                      if(okboit1){
                        if( kko >=*s21b1 && kko <=*s21b2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22b1 && kko <=*s22b2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23b1 && kko <=*s23b2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      }
                      if(ideter[kko] == 2 && ideter[kok] == 2){
                        xmat=xmat+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                      if(okboit1){
                        if( kko >=*s21b1 && kko <=*s21b2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22b1 && kko <=*s22b2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23b1 && kko <=*s23b2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      }
                      if(ideter[kko] == 1 && ideter[kok] == 2){
                        xmat=xmat-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                      if(okboit1){
                        if( kko >=*s21b1 && kko <=*s21b2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22b1 && kko <=*s22b2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23b1 && kko <=*s23b2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                        for(kkio=0;kkio<=*natom-1;kkio++){
                          ideter2[kkio]=ideter[kkio];
                        }
                        ideter2[kko]=2;
                        ideter2[kok]=1;
                        adr_(ideter2, &iaa2);
                        iaa2 = iaa2 - 1;
                        xmat=xmat+valxr[iiii]*valxr[iaa2];
                      if(okboit1){
                        if( kko >=*s21b1 && kko <=*s21b2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22b1 && kko <=*s22b2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23b1 && kko <=*s23b2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                      }
                      if(ideter[kko] == 2 && ideter[kok] == 1){
                        xmat=xmat-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                      if(okboit1){
                        if( kko >=*s21b1 && kko <=*s21b2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22b1 && kko <=*s22b2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23b1 && kko <=*s23b2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                        for(kkio=0;kkio<=*natom-1;kkio++){
                          ideter2[kkio]=ideter[kkio];
                        }
                        ideter2[kko]=1;
                        ideter2[kok]=2;
                        adr_(ideter2, &iaa2);
                        iaa2 = iaa2 - 1;
                        xmat=xmat+valxr[iiii]*valxr[iaa2];
                      if(okboit1){
                        if( kko >=*s21b1 && kko <=*s21b2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22b1 && kko <=*s22b2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23b1 && kko <=*s23b2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                       }
                    }
                }
              }
              for(kko=0;kko<=(*natom/2)-1;kko++){
                  for(kok=(*natom/2);kok<=*natom-1;kok++){
                    if(kok == kko && ideter[kok] != 3){
                      xmat=xmat+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                      if(okboit1){
                        if( kko >=*s21a1 && kko <=*s21a2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22a1 && kko <=*s22a2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23a1 && kko <=*s23a2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4+(3.0/4.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                    }
                    else{
                      if(ideter[kko] == 1 && ideter[kok] == 1){
                        xmat=xmat+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                      if(okboit1){
                        if( kko >=*s21a1 && kko <=*s21a2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22a1 && kko <=*s22a2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23a1 && kko <=*s23a2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      }
                      if(ideter[kko] == 2 && ideter[kok] == 2){
                        xmat=xmat+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                      if(okboit1){
                        if( kko >=*s21a1 && kko <=*s21a2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22a1 && kko <=*s22a2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23a1 && kko <=*s23a2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4+(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      }
                      if(ideter[kko] == 1 && ideter[kok] == 2){
                        xmat=xmat-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                      if(okboit1){
                        if( kko >=*s21a1 && kko <=*s21a2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22a1 && kko <=*s22a2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23a1 && kko <=*s23a2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                        for(kkio=0;kkio<=*natom-1;kkio++){
                          ideter2[kkio]=ideter[kkio];
                        }
                        ideter2[kko]=2;
                        ideter2[kok]=1;
                        adr_(ideter2, &iaa2);
                        iaa2 = iaa2 - 1;
                        xmat=xmat+valxr[iiii]*valxr[iaa2];
                      if(okboit1){
                        if( kko >=*s21a1 && kko <=*s21a2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22a1 && kko <=*s22a2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23a1 && kko <=*s23a2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                      }
                      if(ideter[kko] == 2 && ideter[kok] == 1){
                        xmat=xmat-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                      if(okboit1){
                        if( kko >=*s21a1 && kko <=*s21a2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22a1 && kko <=*s22a2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23a1 && kko <=*s23a2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4-(1.0/2.0)*(valxr[iiii]*valxr[iiii]);
                          }
                        }
                      }
                        for(kkio=0;kkio<=*natom-1;kkio++){
                          ideter2[kkio]=ideter[kkio];
                        }
                        ideter2[kko]=1;
                        ideter2[kok]=2;
                        adr_(ideter2, &iaa2);
                        iaa2 = iaa2 - 1;
//                      if(!mpiid){if(iaa2 > *Iend || iaa2 < *Istart)printf("out iaa2 = %d\n",iaa2);}
                        xmat=xmat+valxr[iiii]*valxr[iaa2];
                      if(okboit1){
                        if( kko >=*s21a1 && kko <=*s21a2){
                          if( kok >=*s21b1 && kok <=*s21b2){
                                xmat2=xmat2+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                      if(okboit2){
                        if( kko >=*s22a1 && kko <=*s22a2){
                          if( kok >=*s22b1 && kok <=*s22b2){
                                xmat3=xmat3+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                      if(okboit3){
                        if( kko >=*s23a1 && kko <=*s23a2){
                          if( kok >=*s23b1 && kok <=*s23b2){
                                xmat4=xmat4+(valxr[iiii]*valxr[iaa2]);
                          }
                        }
                      }
                       }
                    }
                  }
                }
                *xymat=*xymat+xmat;
                *xymat2=*xymat2+xmat2;
                *xymat3=*xymat3+xmat3;
                *xymat4=*xymat4+xmat4;
//              if(mpiid==3)printf(" ii = %d norm = %18f %18f 3 = %18f 4 = %18f\n", ii, *norm2, *norm3, *xymat2, *xymat3);
          }

  ierr = PetscTime(&tt2);CHKERRQ(ierr);
//printf(" norm = %18f weight = %18f weight/N = %18f tmpwe = %18f\n", *norm2, *weight3, *weight3/(*norm2),tmpwe);
//printf(" norm = %18f %18f xymat = %18f %18f\n", *norm2, *norm3, *xymat2, *xymat3);
//ierr = PetscPrintf(PETSC_COMM_WORLD," Time used for the s2 loop: %f\n",tt2-tt1);CHKERRQ(ierr);
}
