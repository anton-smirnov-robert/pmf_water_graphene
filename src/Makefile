FILENAME= main
FC= gfortran
OBJ= constants.o input.o declaration.o matrix_operation.o spce.o graphite.o rpa.o $(FILENAME).o

#FLAGS= -g -fcheck=all -Wall ${MKL_LDFLAGS}

LIBLAPACK=/Users/robert/lapack-3.9.0
LIBBLAS=/Users/robert/lapack-3.9.0
FLAGS= -g -fcheck=all,no-recursion -Wall -L$(LIBLAPACK) -llapack -L$(LIBBLAS) -lrefblas


$(FILENAME).out: $(OBJ)
	$(FC)  $(OBJ) $(FLAGS) -o main 

constants.o: constants.f90
	$(FC)  constants.f90 $(FLAGS) -c

input.o: input.f90
	$(FC)  input.f90 $(FLAGS) -c 

declaration.o:declaration.f90
	$(FC)  declaration.f90 $(FLAGS) -c

matrix_operation.o: matrix_operation.f90
	$(FC)  matrix_operation.f90 $(FLAGS) -c

spce.o: spce.f90
	$(FC)  spce.f90 $(FLAGS) -c

graphite.o: graphite.f90
	$(FC)  graphite.f90 $(FLAGS) -c


rpa.o: rpa.f90
	$(FC)  rpa.f90 $(FLAGS) -c

$(FILENAME).o: $(FILENAME).f90
	$(FC)  $(FILENAME).f90 $(FLAGS) -c


clean:
	@echo "CLEANING..."
	rm -f *.o *.mod main 
