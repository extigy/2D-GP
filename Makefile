#If you do not have a system installed version of libnetcdf-dev or libnetcdff-dev
#replace these with your local install directories
NETCDF = /data/.fs/netcdf
NETCDF-FORTRAN = /data/.fs/netcdf-fortran

FC = gfortran
FCFLAGS = -O3 -march=native -fopenmp -I$(NETCDF)/include -I$(NETCDF)-fortran/include
LDFLAGS = -lm -lnetcdff -lnetcdf -L$(NETCDF)/lib -L$(NETCDF)-fortran/lib

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
