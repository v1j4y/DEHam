#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

void get_proj_9_3h(PetscScalar *valxr, PetscInt *Istart, PetscInt *Iend, int *natom, int iroot, double *projvec, const int natomax);
