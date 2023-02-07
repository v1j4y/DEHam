#include <stdio.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

void get_s2(Vec xr, PetscInt *Istart, PetscInt *Iend, PetscScalar *valxr, int *natom, 
        PetscReal *norm, PetscReal *norm2, PetscReal *norm3, PetscReal *norm4, 
        PetscReal *xymat, PetscReal *xymat2, PetscReal *xymat3, PetscReal *xymat4, PetscReal *weight3,
            int *s21a1, int *s21a2, int *s21b1,  int *s21b2,  int *s22a1,  int *s22a2,  
            int *s22b1, int *s22b2, int *s23a1,  int *s23a2,  int *s23b1, int *s23b2, int *postrou1, int *postrou2, int *postrou3, const int natomax);

void get_s2_mov(Vec, PetscInt *, PetscInt *, PetscScalar *, int *, PetscReal *, PetscReal *,PetscReal *, PetscReal *, PetscReal *, PetscReal *,  PetscReal *, PetscReal *, PetscReal *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *, const int natomax);

void get_s2_cyclic(Vec, PetscInt *, PetscInt *, PetscScalar *, int *, PetscReal *, PetscReal *,PetscReal *, PetscReal *, PetscReal *, PetscReal *,  PetscReal *, PetscReal *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *,
        int *, const int natomax);
