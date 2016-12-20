# DEHam

Double Exchange Hamiltonian: Complete Version
=============================================

O. Dependencies
---------------

1. [PETSc](https://www.mcs.anl.gov/petsc/documentation/installation.html) and [SLEPc](http://slepc.upv.es/documentation/current/docs/instal.htm)

2. [IRPF90](https://github.com/scemama/irpf90)

I. Compiling
------------

1. Export environment variables for PETSc and SLEPc
```shell
export PETSC_DIR=${PATH_TO_PETSC_INSTALLATION}
export SLEPC_DIR=${PATH_TO_SLEPC_INSTALLATION}
```

2. Make the executable
```shell
make ex1
```

II. Using DEHam
---------------

1. The DEHam program requires an input file which 
   has the topology of the Hamiltonian and the various parameters
   as explained below in a sample inputfile:

```python
140					# The total number of determinants
7					# The largest number of non-zero elements per row
2					# The number of processors used in parallel
1					# The number of holes
0					# The isz (ms-1/2) value
1,2,3,1,2,3,4,5,6,7	# The topology of the system is specified here
2,3,4,8,7,6,5,6,7,8	# first and second line contain the two sites linked
1,1,1,2,2,2,2,3,3,3	# third line contains the type of link (1 for t, J 2 for K and 3 for none)
.1430,-0.20,0.0000	# The three types of links this line gives J, K
.1430,-0.20,0.0000	# 
-1.00,0.0,0.00		# This line gives t
```
