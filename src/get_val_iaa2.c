#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "get_val_iaa2.h"

/*
 * This function gets vector from a different processor
 * xr           the full vector
 * iaa2         adresse to get value from
 * getvaliaa2   the value
 */

void get_val_iaa2(Vec xr,long int *iaa2,double *getvaliaa2){
    Vec x; /* initial vector, destination vector */
    VecScatter scatter; /* scatter context */
    IS from, to; /* index sets that define the scatter */
    PetscScalar *values;
    int idx_from[] = {*iaa2}, idx_to[] = {0}; 
    VecCreateSeq(   PETSC_COMM_SELF,1,&x);
    ISCreateGeneral(PETSC_COMM_SELF,1,idx_from,PETSC_COPY_VALUES,&from); 
    ISCreateGeneral(PETSC_COMM_SELF,1,idx_to,  PETSC_COPY_VALUES,&to); 
    printf("in get_val");
    VecScatterCreate(xr,from,x,to,&scatter);
    VecScatterBegin(scatter,xr,x,INSERT_VALUES,SCATTER_FORWARD); 
    VecScatterEnd(scatter,xr,x,  INSERT_VALUES,SCATTER_FORWARD); 
    VecGetArray(x,&values);
    *getvaliaa2 = values[0];
    ISDestroy(&from);
    ISDestroy(&to);
    VecScatterDestroy(&scatter);
}
