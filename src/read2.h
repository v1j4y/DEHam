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
	long int l1[900];
	long int l2[900];
	long int ktyp[900];
	double	xjjz[900];
	double	xjjxy[900]; 
	double	xtt[900]; 
	double  E[900]; 
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
	int postrou1;
	int postrou2;
	int postrou3;
	long int fix_trou1;
	long int fix_trou2;
    long int print_wf;

} Data ;

void Data_new(FILE* , Data* );
