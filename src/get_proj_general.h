#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "functions.h"

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
                      double *projmatrix);
