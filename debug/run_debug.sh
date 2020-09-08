ifort -g -c debug.f90 -traceback -mkl
ifort -g -o ./test.out debug.o -traceback -mkl
./test.out
rm *.mod
