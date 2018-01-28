#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

void get_1rdm(PetscScalar *, PetscInt *, PetscInt *, int *,  PetscReal *, const int natomax);
void get_2rdm(PetscScalar *valxr, PetscInt *Istart, PetscInt *Iend, int *natom, PetscReal *trace2rdm, double densmat2[*natom][*natom][*natom][*natom], const int natomax);
