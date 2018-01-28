#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include <petscsys.h>
#include <slepceps.h>

_Bool to_bool(const char* str);

typedef struct {
	PetscInt n;
	long int nnz,npar;
	long int ntrou,isz;
        _Bool FAM1;
	long int l1[700];
	long int l2[700];
	long int ktyp[700];
	double	xjjz[700];
	double	xjjxy[700]; 
	double	xtt[700]; 
	long int nroots;
	int natom;
	int s21a1;
	int s21a2;
	int s21b1;
	int s21b2;
	int s22a1;
	int s22a2;
	int s22b1;
	int s22b2;
	int s23a1;
	int s23a2;
	int s23b1;
	int s23b2;
	long int postrou;
	long int fix_trou1;
	long int fix_trou2;

} Data ;

void Data_new(FILE* , Data* );
