# algencan_adolc
Algencan 3.1.1 using Automatic Differenciation library ADOL-C 2.7.3

In this file we present the complete instalation of Algencan 3.1.1 and ADOL-C 2.7.3 and how to use those libraries together in Linux Ubuntu.

Algencan 3.1.1 is a nonlinear programming solver and ADOL-C 2.7.3 is used to evaluate sparses Hessians in Algencan methods.

# INSTALLATION

Before installing ADOL-C, we first need to install two libraries: 1. Boost 1.75.0; 2. Colpack 1.0.10.

## ADOL-C 2.7.3 INSTALLATION

### Boost 1.75.0 INSTALLATION

1 Download boost_1_75_0.tar.bz2 in https://www.boost.org/users/download/
2 Extract the file in the desired directory
> tar --bzip2 -xf /path/to/boost_1_75_0.tar.bz2
3 Open boost directory > cd path/to/boost_1_75_0
4 Run bootstrap.sh > ./bootstrap.sh --prefix=path/to/installation/prefix
5 Run b2 to generate external libraries in \lib\ >./b2 install

### Colpack 1.0.10 INSTALLATION

Colpack is used by ADOL-C to treat sparse matrices.

1 Download Colpack in https://github.com/CSCsw/ColPack/releases
2 Extract the file
3 Run the following commands in Colpack directory:
    > autoreconf -vif
    > ./configure --prefix=/path/to/install/
    > make
    > make install

### ADOL-C INSTALLATION

If Boost and Colpack are installed successfully, we may now install ADOL 2.7.3

1 Download ADOL-C in https://github.com/coin-or/ADOL-C
2 Extract the file
3 Run the following commands in ADOL-C directory:
    > autoreconf -fi
    > ./configure --enable-sparse --with-boost=PATHTOBOOST --with-colpack=PATHCOLPACK
    > make
    > make install

where PATHTOBOOST and PATHCOLPACK is the path to Boost and Colpack directories, respectively.


If the installation is successfull, ADOL-C external libraries are generated in /home/USER/adolc_base, where /home/USER is the user directory.

## ALGENCAN 3.1.1 INSTALLATION

1 Download ALGENCAN 3.1.1 in https://www.ime.usp.br/~egbirgin/tango/codes.php
2 Extract the file 
3 Run in ALGENCAN directory > make

If the installation is succeed, external library libalgencan.a is generated in PATHTOALGENCAN/lib/

# Linking and compiling

## Compiling codade via terminal

To compile an ADOL-C code main.cpp use the following command:

> g++ -w -I/home/USER/adolc_base/include -o main main.cpp -Wl,--rpath -Wl,/home/USER/adolc_base/lib64 -L/home/USER/adolc_base/lib64 -ladolc}

To compile an ALGENCAN-C code main.c use the following command:

> gcc -O3 main.c -L\$ALGENCAN/lib -lalgencan -lgfortran -lm -o algencan}

To compile an ALGENCAN-C code algencan.c together with ADOL-C code adolc.cpp use the following commands:

> g++ -w -c -I/home/USER/adolc_base/include adolc.cpp \ -Wl,--rpath -Wl,/home/USER/adolc_base/lib64 \ -L/home/USER/adolc_base/lib64/ -ladolc

> gcc -w -c algencan.c \ -L/$ALGENCAN/lib/ -lalgencan -lm \ -L/usr/lib/gcc/x86_64-linux-gnu/5/ -lstdc++ -lgfortran

> g++ -O3 -Wall -I/home/USER/adolc_base/include -o main algencan.o adolc.o \ -Wl,--rpath -Wl,/home/USER/adolc_base/lib64 \ -L/home/USER/adolc_base/lib64/ -ladolc \ -L/$ALGENCAN/lib/ -lalgencan -lm \ -L/usr/lib/gcc/x86_64-linux-gnu/5/ -lstdc++ -lgfortran

Run the main program:

> ./main

## Linking and compiling using Codeblocks 20.03

1 Create C/C++ Console Application project.
2 In Build Options >> Linker Settings, add the following libraries:

> adolc.so
> libColPack.so
> libalgencan.a
> gfortran
3 In Build Options$>>$Search Directories, add paths:
> /home/USER/adolc\_base/include
> /home/USER/adolc\_base/lib64
in Compiler, Linker, and Resource compiler
