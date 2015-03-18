# Installation #

## Platform ##

The program has been tested on Linux, HP Unix and Mac OS X, but should also be possible to migrate to other platforms. The instructions here assume you are using Ubuntu Linux.

## Required Software ##
  1. BLAS and LAPACK
  1. Lapack++
  1. ARPACK

## Step by Step ##
### 1.Install BLAS and LAPACK ###
BLAS and LAPACCK are highly optimized FORTRAN routines for linear algebra operations. If you are using a HPC cluster, normally they are already installed within packages like: Intel MKL and ATLAS. We can go to next step directly.

If you have Ubuntu or Linux Mint on your computer you can install ATLAS by running this command:
```
sudo apt-get install libatlas3gf-base libatlas-base-dev
```
You can also follow [this guide](http://math-atlas.sourceforge.net/atlas_install/) to compile by yourself.

### 2. Install gfortran ###
```
sudo apt-get install gfortran
```

### 3. Install Lapack++ ###
Lapack++ provides some c++ classes based on BLAS and LAPACK.

  1. Download the source package from [here](http://sourceforge.net/projects/lapackpp/files/).
  1. Extract the files and enter the directory.
  1. Suppose you want to install lapack++ in your home directory under usr, run the following commands:
```
./configure --enable-static=yes --enable-shared=no --prefix=$HOME/usr/
make
make install
```
> if you want to use Intel MKL library run the flowing commands before the running the above commands:
```
export CC=icc
export CXX=icpc
export LIBS=$MKL_LIB
```
> $MKL\_LIB is the linking of tags for the MKL libraries, and you should consult your system administrators or your system wiki's for this information.

### 4. Install ARPACK ###
ARPACK is a FORTRAN package for spare matrix diagnalization.
On Ubuntu and Linux Mint you can install it by running
```
sudo apt-get install libarpack2 libarpack2-dev
```