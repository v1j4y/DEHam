#include "functions.h"

// Function to calculate factorial of a given number
int factorial(int n)
{
    int fact = 1;
    for (int i = 1; i <= n; i++)
    {
        fact *= i;
    }
    return fact;
}


char* intToCharPtr(int num)
{
    // Allocate memory for the character pointer
    char* ptr = (char*)malloc(sizeof(char) * 10);

    // Convert the integer to a string
    sprintf(ptr, "%d", num);

    // Return the character pointer
    return ptr;
}

int charPtrToInt(char *str)
{
    // Initialize result
    int res = 0;

    // Iterate through all characters of input string and
    // update result
    for (int i = 0; str[i] != '\0'; ++i){
        res = res*10 + str[i] - '0';
    }

    // return result.
    return res;
}

// Create a new dictionary
struct Dictionary *dictionary_new() {
    struct Dictionary *d = malloc(sizeof(struct Dictionary));
    d->size = 0;
    d->entries = NULL;
    return d;
}

// Add a key-value pair to the dictionary
void dictionary_add(struct Dictionary *d, char *key, char *value) {
    // Allocate memory for the new entry
    struct Entry *e = malloc(sizeof(struct Entry));
    e->key = key;
    e->value = value;

    // Allocate memory for the new array of entries
    struct Entry **new_entries = malloc(sizeof(struct Entry *) * (d->size + 1));

    // Copy the old entries into the new array
    for (int i = 0; i < d->size; i++) {
        new_entries[i] = d->entries[i];
    }

    // Add the new entry to the end of the array
    new_entries[d->size] = e;

    // Free the old array and replace it with the new one
    free(d->entries);
    d->entries = new_entries;

    // Increment the size of the dictionary
    d->size++;
}

// Look up a value in the dictionary
char *dictionary_get(struct Dictionary *d, char *key) {
    for (int i = 0; i < d->size; i++) {
        if (strcmp(d->entries[i]->key, key) == 0) {
            return d->entries[i]->value;
        }
    }
    return NULL;
}

// Free the memory used by the dictionary
void dictionary_free(struct Dictionary *d) {
    for (int i = 0; i < d->size; i++) {
        free(d->entries[i]);
    }
    free(d->entries);
    free(d);
}

// Function to calculate outer product of two vectors
// Input:
//  a - Vecetor
//  b - Vector
//  m - Dim(a)
//  n - Dim(b)
// Output:
//  matrix - Outer product
void outerProduct(double a[], double b[], int m, int n, double matrix[][n])
{
    int i, j;

    // Outer product of two vectors
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            matrix[i][j] = a[i] * b[j];
    }
}

// Function to flatten a matrix to an array
// Input:
//  rows   - rows of A
//  rows   - rows of A
//  matrix - Matrix A
//  array  - Flattened array
void flattenMatrix(int rows, int columns, double matrix[rows][columns], double array[rows * columns])
{
    int i, j, k = 0;

    // Iterate over the matrix
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < columns; j++)
        {
            // Store the elements of matrix in array
            array[k++] = matrix[i][j];
        }
    }
}

// Function to print array
void printArray(double arr[], int size)
{
    int i;
    for (i = 0; i < size; i++)
        printf("%15f ", arr[i]);
    printf("\n");
}

// Function to print matrix
void printMatrix(int row, int col, double matrix[][col])
{
    int i, j;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
            printf("%15f ", matrix[i][j]);
        printf("\n");
    }
}


// Function to convert base 6 to base 10
int base6ToBase10(int n)
{
    int base = 1;
    int result = 0;

    while (n > 0)
    {
        result += n % 10 * base;
        n /= 10;
        base *= 6;
    }

    return result;
}

//Function to convert base 10 to base 6
int base10ToBase6(int n)
{
    int rem, i=1, octal=0;
    while (n!=0)
    {
        rem=n%6;
        n/=6;
        octal+=rem*i;
        i*=10;
    }
    return octal;
}

/* This function takes an integer as an argument and prints the digits of the number */
void printDigits(int num)
{
    // Declare a variable to store the remainder
    int remainder;

    // If the number is 0, print 0 and return
    if (num == 0)
    {
        printf("0");
        return;
    }

    // Loop until the number is 0
    while (num > 0)
    {
        // Get the remainder of the number
        remainder = num % 10;

        // Print the remainder
        printf("%d", remainder);

        // Divide the number by 10
        num /= 10;
    }
}

/*
This function takes an integer as an argument and stores its digits in an array.
The array is passed as a parameter to the function.

@param int num: The number whose digits are to be stored
@param int arr[]: The array in which the digits are to be stored
@return void
*/
void storeDigits(int num, int arr[], int *sze)
{
    int i = 0;
    while (num > 0) {
        arr[i] = num % 10;
        num /= 10;
        i++;
    }
    (*sze) = i;
}

// Find dimension of model space
int get_nstates(int sze_in, int ms)
{
    int i, j, idx, sze;
    int arr[SZE_MAX];
    int arr_all[DIM_MAX * SZE_MAX];
    int sum, nstate;

    idx = 0;
    arr_all[0]=0;
    for(j=0;j<DIM_MAX;++j){
        sum = 0;
        nstate = 1;
        storeDigits(base10ToBase6(j),arr,&sze);
        for(i=0;i<sze;++i){
            sum += arr[i];
            nstate *= factorial(5)/(factorial(arr[i])*factorial(5-arr[i]));
        }
        if(sum == ms && sze <= sze_in){
            for(i=0;i<sze;++i)
                arr_all[idx*SZE_MAX + i]=0;
            for(i=0;i<sze;++i)
                arr_all[idx*SZE_MAX + i]=arr[i];
            idx++;
        }
    }

    return(idx);
}

// Prepare the factors
// for the spatial part i.e. the Hueckel factors
int prepareHueckelFactors(int nsite, double *factor, int size){
    double *a;
    double *b;

    int m = 3;
    int n = 3;
    // Initialize
    a = (double *)malloc(m * sizeof(double));
    b = (double *)malloc(n * sizeof(double));
    a[0] = 1./2.;
    a[1] = sqrt(2.)/2.;
    a[2] = 1./2.;
    b[0] = 1./2.;
    b[1] = sqrt(2.)/2.;
    b[2] = 1./2.;

    double c[m][n];
    double e[m*n];

    int i;
    int isite = 1;
    while(isite < nsite){

        double c[m][n];
        double e[m*n];

        // First outer product
        outerProduct(a, b, m, n, c);
        flattenMatrix(m, n, c, e);

        free(a);
        // Copy vector
        m = m*n;
        a = (double *)malloc(m * sizeof(double));
        memcpy (a, e, m*sizeof(double));

        isite++;
    }


    // Copyt to factor
    assert ( m == size );
    memcpy (factor, a, size*sizeof(double));

    free(a);
    free(b);
    return 0;
}

// Main function to prepare dictionary
int prepare_dictionary(struct Dictionary *d, int sze_in, int ms)
{
    int i, j, idx, sze;
    int arr[SZE_MAX];
    int arr_all[DIM_MAX * SZE_MAX];
    int sum, nstate;

    idx = 0;
    arr_all[0]=0;
    for(j=0;j<DIM_MAX;++j){
        sum = 0;
        nstate = 1;
        storeDigits(base10ToBase6(j),arr,&sze);
        for(i=0;i<sze;++i){
            sum += arr[i];
            //printf(" %d ",arr[i]);
            nstate *= factorial(5)/(factorial(arr[i])*factorial(5-arr[i]));
        }
        //printf("\n");
        if(sum == ms && sze <= sze_in){
          //printf("\t %s \n",intToCharPtr(base10ToBase6(j)));
            for(i=0;i<sze;++i)
                arr_all[idx*SZE_MAX + i]=0;
            for(i=0;i<sze;++i)
                arr_all[idx*SZE_MAX + i]=arr[i];
            // Add some key-value pairs
            //dictionary_add(d, intToCharPtr(idx), intToCharPtr(base10ToBase6(j)));
            dictionary_add(d, intToCharPtr(base10ToBase6(j)), intToCharPtr(nstate));
            idx++;
        }
    }

    return(idx);
}

// Main function to prepare adressing dictionary
int prepare_dictionary_for_adressing(struct Dictionary *d, int sze_in, int ms)
{
    int i, j, idx, sze;
    int arr[SZE_MAX];
    int arr_all[DIM_MAX * SZE_MAX];
    int sum, nstate;

    idx = 0;
    arr_all[0]=0;
    for(j=0;j<DIM_MAX;++j){
        sum = 0;
        nstate = 1;
        storeDigits(base10ToBase6(j),arr,&sze);
        for(i=0;i<sze;++i){
            sum += arr[i];
            nstate *= factorial(5)/(factorial(arr[i])*factorial(5-arr[i]));
        }
        if(sum == ms && sze <= sze_in){
            for(i=0;i<sze;++i)
                arr_all[idx*SZE_MAX + i]=0;
            for(i=0;i<sze;++i)
                arr_all[idx*SZE_MAX + i]=arr[i];
            // Add some key-value pairs
            dictionary_add(d, intToCharPtr(base10ToBase6(j)), intToCharPtr(idx));
            idx++;
        }
    }

    return(idx);
}
