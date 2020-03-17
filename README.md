[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.20450.svg)](http://dx.doi.org/10.5281/zenodo.20450)

# DEHam

Double Exchange Hamiltonian: Complete Version
=============================================

(under GNU GENERAL PUBLIC LICENSE v2)

_Dependencies_
---------------

  1. [PETSc](https://www.mcs.anl.gov/petsc/documentation/installation.html) and [SLEPc](http://slepc.upv.es/documentation/instal.htm)

  2. [IRPF90](https://github.com/scemama/irpf90)

_Compiling_
------------

  1. Export environment variables for PETSc and SLEPc

```shell
export PETSC_DIR=${PATH_TO_PETSC_INSTALLATION}
export SLEPC_DIR=${PATH_TO_SLEPC_INSTALLATION}
export C_INCLUDE_PATH+=:$PETSC_DIR/include/:$SLEPC_DIR/include:$PETSC_DIR/arch-linux2-c-debug/include/:$SLEPC_DIR/arch-linux2-c-debug/include
# The "arch-linux2-c-debug" directory can have different names depending on PETSC and SLEPC installation procedure.
```


  2. Make the executable

```shell
make ex1
```

[![asciicast](https://asciinema.org/a/Ng3tSNDoWBkV5C9ZYvbxCW43B.png)](https://asciinema.org/a/Ng3tSNDoWBkV5C9ZYvbxCW43B)


_Using DEHam_
---------------

  1. The DEHam program requires an input file which 
   has the topology of the Hamiltonian and the various parameters
   as explained below in a sample inputfile:

```python
8     # The number of orbitals (total)
140   # The largest number of non-zero elements per row (Multiple of Ndet)
1     # The total number of processors used in parallel (Multiple of Ndet)
1     # The number of holes
0     # The isz (ms-1/2) value
true     # Restrict the hole to the 1'st (i.e. half of natom) Family of states. *false* for no restrictions
1,2,3,1,2,3,4,5,6,7	# The topology of the system is specified here
2,3,4,8,7,6,5,6,7,8	# first and second line contain the two sites linked
1,1,1,2,2,2,2,3,3,3	# third line contains the type of link (1 for t or J, 2 for K and 3 for none)
.1430,-0.20,0.0000	# The three types of links this line gives J, K
.1430,-0.20,0.0000	# 
-1.00,0.0,0.00		# This line gives t
0.,0.,0.,0.,0.,0.,0.,0.,0. # Energy of each orbital + one extra term
2                   # The total number of roots
1                   # I   The position of the first
1                   # I   SBox
1                   # I
1                   # I
1                   # II  The positions of the second
1                   # II  SBox
1                   # II
1                   # II
1                   # III
1                   # III The positions of the third
1                   # III SBox
1                   # III
1                   # positio of the hole
0                   # fix the position of the first hole during the CI
0                   # fix the position of the second hole during the CI
0                   # Print the wavefunction. It is stored in the FIL666 file after the run
```

  2. running DEHam

```shell
mpiexec -n [nprocs] ./ex1 inpfile 
```

_Sample Application_
--------------------

A 2D t-J model Hamiltonian description and setup for using DEHam to solve for few low lying states
is provided in the notbooks folder. Please have a look about the details of using DEHam to study
t-J Hamiltonians.

![4x4 Lattice](https://raw.githubusercontent.com/v1j4y/DEHam/master/notebook/graph.png)


_Publications using this code_
-------------------------------

  1. High-Spin Chains and Crowns from Double-Exchange Mechanism [doi:10.3390/cryst6040039](http://www.dx.doi.org/10.3390/cryst6040039)
