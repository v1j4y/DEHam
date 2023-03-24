#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

void get_proj(PetscScalar *valxr, PetscInt *Istart, PetscInt *Iend, int *natom, int iroot, double *projvec, const int natomax);
