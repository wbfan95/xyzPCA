FC  = ifort
OPT = -mkl -qopenmp
EXE = xyzPCA

make:
	$(FC) -c module.f90 $(OPT)
	$(FC) -c main.f90 $(OPT)
	$(FC) -o $(EXE) main.o module.o $(OPT)

clean:
	rm -f *.o *.mod *.out $(EXE)
