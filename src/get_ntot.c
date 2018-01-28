#include "get_ntot.h"

int get_ntot(_Bool FAM1, int natom, long int isz, long int ntrou, long int fix_trou1, long int fix_trou2){
  int            tnt1, tnt2;
  int            natom2;
  if(FAM1){
      if(fix_trou1 == fix_trou2){
          natom2 = natom/2;
      }
      else{
          natom2 = fix_trou2 - fix_trou1;
      }
  }
  else{
      natom2 = natom;
  }

  tnt1 =         (int)ceil(exp(lgamma((double)(natom2+1)) - (lgamma((double)(natom2-ntrou+1)) + lgamma((double)(ntrou+1)))));
  printf("%10.5f | tnt1=%d\n",exp(lgamma((double)(natom2+1)) - (lgamma((double)(natom2-ntrou+1)) + lgamma((double)(ntrou+1)))),tnt1);
  int            nalpha, nbeta;

    if((((natom-ntrou) + 2*isz) % 2) == 0){
        nalpha=(natom-ntrou+2*isz)/2;
        nbeta=(natom -ntrou-2*isz)/2;
        if(((natom-ntrou)/2) == isz){
            nbeta=0;
        }
    }
    else{
        nalpha=(natom-ntrou+2*isz+1)/2;
        nbeta=(natom -ntrou-2*isz-1)/2;
        if(((natom-ntrou+1)/2) == isz){
            nbeta=0;
        }
    }

  tnt2 =         (int)ceil(exp(lgamma((double)(natom-ntrou+1)) - (lgamma((double)(nalpha+1)) + lgamma((double)(nbeta+1)))));
  printf("natom2=%d fix_trou1=%d fix_trou2=%d nalpha=%d nbeta=%d | | %d %d ntot=%d\n",natom2, fix_trou1, fix_trou2, nalpha, nbeta, tnt1, tnt2, tnt1*tnt2);
  return tnt1*tnt2;
}
