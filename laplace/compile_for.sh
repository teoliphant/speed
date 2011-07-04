FFLAGS="-O3 -ffast-math -funroll-loops"
gfortran $FFLAGS -o laplace_for.o -c laplace_for.f90
gfortran -o laplace_for laplace_for.o
