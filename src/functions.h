#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define SZE_MAX 8
//#define MS 10
#define SZE_HERE 3
#define DIM_MAX 216

// This function takes in a number and returns its factorial
int factorial(int n);

// Structure to store a key-value pair
struct Entry {
    char *key;
    char *value;
};

// Structure to store a dictionary
struct Dictionary {
    int size;
    struct Entry **entries;
};

char* intToCharPtr(int num);
int charPtrToInt(char *str);

// Create a new dictionary
struct Dictionary *dictionary_new();

// Add a key-value pair to the dictionary
void dictionary_add(struct Dictionary *d, char *key, char *value);

// Look up a value in the dictionary
char *dictionary_get(struct Dictionary *d, char *key);

// Free the memory used by the dictionary
void dictionary_free(struct Dictionary *d);

// Function to calculate outer product of two vectors
// Input:
//  a - Vecetor
//  b - Vector
//  m - Dim(a)
//  n - Dim(b)
void outerProduct(double a[], double b[], int m, int n, double matrix[][n]);

// Function to flatten a matrix to an array
// Input:
//  rows   - rows of A
//  rows   - rows of A
//  matrix - Matrix A
//  array  - Flattened array
void flattenMatrix(int rows, int columns, double matrix[rows][columns], double array[rows * columns]);

// Function to print array
void printArray(double arr[], int size);

// Function to print matrix
void printMatrix(int row, int col, double matrix[][col]);

// Function to convert base 6 to base 10
int base6ToBase10(int n);

//Function to convert base 10 to base 6
int base10ToBase6(int n);

/* This function takes an integer as an argument and prints the digits of the number */
void printDigits(int num);

/*
This function takes an integer as an argument and stores its digits in an array.
The array is passed as a parameter to the function.

@param int num: The number whose digits are to be stored
@param int arr[]: The array in which the digits are to be stored
@return void
*/
void storeDigits(int num, int arr[], int *sze);

// Find dimension of model space
int get_nstates(int sze, int ms);

// Prepare the factors
// for the spatial part i.e. the Hueckel factors
int prepareHueckelFactors(int nsite, double *factor, int size);

// Main function to prepare dictionary
int prepare_dictionary(struct Dictionary *d, int sze, int ms);

// Main function to prepare adressing dictionary
int prepare_dictionary_for_adressing(struct Dictionary *d, int sze, int ms);
