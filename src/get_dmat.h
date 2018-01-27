#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

void get_1rdm(PetscScalar *, PetscInt *, PetscInt *, int *,  PetscReal *);
void get_2rdm(PetscScalar *, PetscInt *, PetscInt *, int *,  PetscReal *, double ****);
