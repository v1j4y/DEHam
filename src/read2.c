#include "read2.h"

/* Compare two strings s1 and s2, assuming s1 is terminated
 * by \n or a NULL, and s2 is terminated by a NULL. A match
 * returns 1, a non-match returns 0.
 */
int
strcmpst1nl (const char * s1, const char * s2)
{
  char s1c;
  do
    {
      s1c = *s1;
      if (s1c == '\n')
          s1c = 0;
      if (s1c != *s2)
          return 0;
      s1++;
      s2++;
    } while (s1c); /* already checked *s2 is equal */
  return 1;
}

void Data_new(FILE* file, Data* dat) { 

//ata* dat = (Data*)malloc(sizeof(Data));
    char line[256];
	char *token;
	const char *delim=",";
	size_t arrayIdx=0;
	int count=0;

	while (fgets(line, sizeof(line), file)) {
	
	    /* note that fgets doesn't strip the terminating \n, checking its
	       presence would allow to handle lines longer than sizeof(line) */
	  if (count != 32){
			count++;
			switch(count){
				case 1:
					dat->natom=atol(line);
					break;
				case 2:
					dat->nnz=atol(line);
					break;
				case 3:
					dat->npar=atol(line);
					break;
				case 4:
					dat->ntrou=atol(line);
					break;
				case 5:
					dat->isz=atol(line);
					break;
				case 6:
					//dat->FAM1 = to_bool(line);
          dat->FAM1 = 0;
          //dat->FAM1 = line && strcmp(line,"true")==0;
          dat->FAM1 = strcmpst1nl(line,"true");
					break;
				case 7:
					arrayIdx=0;
					for (token = strtok(line, delim); token != NULL; token = strtok(NULL, delim))
					{
					  double val;
					  char *unconverted;
					  /**
					   * Convert the next token to a float value
					   */
					  val = strtol(token, &unconverted, 10);
					  if (!isspace(*unconverted) && *unconverted != 0)
					  {
					    /**
					     * Bad input string.  Again, we just bail.
					     */
					    fprintf(stderr, "\"%s\" is not a valid floating-point number\n", token);
					    break;
					  }
					  else
					  {
					    dat->l1[arrayIdx++] = val;
					  }
					}
					break;
				case 8:
					arrayIdx=0;
					for (token = strtok(line, delim); token != NULL; token = strtok(NULL, delim))
					{
					  double val;
					  char *unconverted;
					  /**
					   * Convert the next token to a float value
					   */
					  val = strtol(token, &unconverted, 10);
					  if (!isspace(*unconverted) && *unconverted != 0)
					  {
					    /**
					     * Bad input string.  Again, we just bail.
					     */
					    fprintf(stderr, "\"%s\" is not a valid floating-point number\n", token);
					    break;
					  }
					  else
					  {
					    dat->l2[arrayIdx++] = val;
					  }
					}
					break;
				case 9:
					arrayIdx=0;
					for (token = strtok(line, delim); token != NULL; token = strtok(NULL, delim))
					{
					  double val;
					  char *unconverted;
					  /**
					   * Convert the next token to a float value
					   */
					  val = strtol(token, &unconverted, 10);
					  if (!isspace(*unconverted) && *unconverted != 0)
					  {
					    /**
					     * Bad input string.  Again, we just bail.
					     */
					    fprintf(stderr, "\"%s\" is not a valid floating-point number\n", token);
					    break;
					  }
					  else
					  {
					    dat->ktyp[arrayIdx++] = val;
					  }
					}
					break;
				case 10:
					arrayIdx=0;
					for (token = strtok(line, delim); token != NULL; token = strtok(NULL, delim))
					{
					  double val;
					  char *unconverted;
					  /**
					   * Convert the next token to a float value
					   */
					  val = strtod(token, &unconverted);
					  if (!isspace(*unconverted) && *unconverted != 0)
					  {
					    /**
					     * Bad input string.  Again, we just bail.
					     */
					    fprintf(stderr, "\"%s\" is not a valid floating-point number\n", token);
					    break;
					  }
					  else
					  {
					    dat->xjjz[arrayIdx++] = val;
					  }
					}
					break;
				case 11:
					arrayIdx=0;
					for (token = strtok(line, delim); token != NULL; token = strtok(NULL, delim))
					{
					  double val;
					  char *unconverted;
					  /**
					   * Convert the next token to a float value
					   */
					  val = strtod(token, &unconverted);
					  if (!isspace(*unconverted) && *unconverted != 0)
					  {
					    /**
					     * Bad input string.  Again, we just bail.
					     */
					    fprintf(stderr, "\"%s\" is not a valid floating-point number\n", token);
					    break;
					  }
					  else
					  {
					    dat->xjjxy[arrayIdx++] = val;
					  }
					}
					break;
				case 12:
					arrayIdx=0;
					for (token = strtok(line, delim); token != NULL; token = strtok(NULL, delim))
					{
					  double val;
					  char *unconverted;
					  /**
					   * Convert the next token to a float value
					   */
					  val = strtod(token, &unconverted);
					  if (!isspace(*unconverted) && *unconverted != 0)
					  {
					    /**
					     * Bad input string.  Again, we just bail.
					     */
					    fprintf(stderr, "\"%s\" is not a valid floating-point number\n", token);
					    break;
					  }
					  else
					  {
					    dat->xtt[arrayIdx++] = val;
					  }
					}
					break;
				case 13:
					arrayIdx=0;
					for (token = strtok(line, delim); token != NULL; token = strtok(NULL, delim))
					{
					  double val;
					  char *unconverted;
					  /**
					   * Convert the next token to a float value
					   */
					  val = strtod(token, &unconverted);
					  if (!isspace(*unconverted) && *unconverted != 0)
					  {
					    /**
					     * Bad input string.  Again, we just bail.
					     */
					    fprintf(stderr, "\"%s\" is not a valid floating-point number\n", token);
					    break;
					  }
					  else
					  {
					    dat->E[arrayIdx++] = val;
					  }
					  dat->xrep = val;
					}
					break;
                case 14:
					dat->nroots=atol(line);
					break;
                case 15:
					dat->s21a1=atol(line);
					break;
                case 16:
					dat->s21a2=atol(line);
					break;
                case 17:
					dat->s21b1=atol(line);
					break;
                case 18:
					dat->s21b2=atol(line);
					break;
                case 19:
					dat->s22a1=atol(line);
					break;
                case 20:
					dat->s22a2=atol(line);
					break;
                case 21:
					dat->s22b1=atol(line);
					break;
                case 22:
					dat->s22b2=atol(line);
					break;
                case 23:
					dat->s23a1=atol(line);
					break;
                case 24:
					dat->s23a2=atol(line);
					break;
                case 25:
					dat->s23b1=atol(line);
					break;
                case 26:
					dat->s23b2=atol(line);
					break;
                case 27:
					dat->postrou1=atol(line);
					break;
                case 28:
					dat->postrou2=atol(line);
					break;
                case 29:
					dat->postrou3=atol(line);
					break;
                case 30:
					dat->fix_trou1=atol(line);
					break;
                case 31:
					dat->fix_trou2=atol(line);
					break;
				case 32:
					dat->print_wf = atol(line);
					break;
                default:
                    printf("Done reading file\n");
                    break;
		} /* end of switch */

	  } /* end of the input file */

	} /* end of loop */

//return dat;
}

_Bool to_bool(const char* str) {
    PetscBool strflg=PETSC_FALSE;
    PetscStrcmp("true\n",str, &strflg);
    if(strflg != PETSC_TRUE) PetscStrcmp("True\n",str, &strflg);
    if(strflg != PETSC_TRUE) PetscStrcmp("TRUE\n",str, &strflg);
    if(strflg != PETSC_TRUE) PetscStrcmp("true",str, &strflg);
    if(strflg != PETSC_TRUE) PetscStrcmp("True",str, &strflg);
    if(strflg != PETSC_TRUE) PetscStrcmp("TRUE",str, &strflg);
    return strflg;
}

/*
int main(int argc, char* argv[])
{
    char const* const fileName = argv[1];
    FILE* file = fopen(fileName, "r");

	Data *getdata=Data_new(file);

    fclose(file);

	printf("n=%d,nnz=%d,npar=%d\n",getdata->n,getdata->nnz,getdata->npar);
	for(int i=0;i<7;i++){
		printf("l1(%d)=%d l2(%d)=%d l3(%d)=%d\n",i,getdata->l1[i],i,getdata->l2[i],i,getdata->l3[i]);
	}
	for(int i=0;i<7;i++){
		printf("xjjz(%d)=%f xjjxy(%d)=%f xtt(%d)=%f\n",i,getdata->xjjz[i],i,getdata->xjjxy[i],i,getdata->xtt[i]);
	}

    return 0;
}
*/
