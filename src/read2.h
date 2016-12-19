#include <stdio.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

typedef struct {
	PetscInt n;
	long int nnz,npar;
	long int ntrou,isz;
	long int l1[700];
	long int l2[700];
	long int ktyp[700];
	double	xjjz[700];
	double	xjjxy[700]; 
	double	xtt[700]; 

} Data ;

void Data_new(FILE* , Data* );
