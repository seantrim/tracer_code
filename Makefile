#!/bin/bash

#CEDAR
#ifort -O3 -mkl=parallel -qopenmp *.f90
#ifort -O3 -mkl=parallel -qopenmp *.f90

#SCINET
#ifort -O3 -mkl=parallel -openmp *.f90
#ifort -O3 -mkl=parallel -openmp *.f90

#LINUX MINT
export OMP_STACKSIZE=10M
gfortran -fcheck=bounds -O3 -fopenmp -m64 -I${MKLROOT}/include *.f90  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
gfortran -fcheck=bounds -O3 -fopenmp -m64 -I${MKLROOT}/include *.f90  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

#gfortran -O3 -fopenmp -m64 -I${MKLROOT}/include *.f90  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
#gfortran -O3 -fopenmp -m64 -I${MKLROOT}/include *.f90  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

echo "Make Complete"
