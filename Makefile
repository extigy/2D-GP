
FC = gfortran

FCFLAGS = -O3 -fopenmp

LDFLAGS = -I/usr/include/ -L/usr/lib/

PROGRAMS = gp mu


all: $(PROGRAMS)

gp: params.o bitmap.o utils.o rhs.o potential.o

mu: mu.o params.o bitmap.o utils.o rhs.o potential.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ -lnetcdff -lnetcdf

%.o: %.f95
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS)

%.o: %.F95
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS)

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD *.in

veryclean: clean
	rm -f *~ $(PROGRAMS)
