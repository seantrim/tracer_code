# tracer_code
Mantle convection code featuring the tracer ratio method.
Examples of how to compile the code using gfortran and ifort are included in the make file (Makefile). In the case of gfortran, intel's mkl libraries must first be installed (free through intel's website) and linked to the compiler. Input parameters must be specified in the file main.f90 before compilation.
Output files are labelled T10000, T10001 and so on and are written at fixed model time intervals. These files contain temperatures, velocities, and compositions. The file stats.dat, contains statistics output at every model time step.
