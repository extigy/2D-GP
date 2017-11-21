#If you do not have a system installed version of libnetcdf-dev or libnetcdff-dev
#replace these with your local install directories
NETCDF = /usr
NETCDF-FORTRAN = /usr

FC = gfortran
FCFLAGS = -O3 -Wall -march=native -fopenmp -I$(NETCDF)/include -I$(NETCDF-FORTRAN)/include
LDFLAGS = -lm -lnetcdff -lnetcdf -L$(NETCDF)/lib -L$(NETCDF-FORTRAN)/lib

PROGRAMS = gp

all: $(PROGRAMS) params.in ic.in

gp: params.o output.o bitmap.o output.o utils.o rhs.o potential.o 

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) 

%.o: %.f95
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS)

%.o: %.F95
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS)

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD *.in

veryclean: clean
	rm -f *~ $(PROGRAMS)
