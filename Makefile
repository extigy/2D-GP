FC = $(shell nf-config --fc)
INCLUDES = $(shell nf-config --fflags) $(foreach d,$(subst :, ,$(CPATH)),-I$d) -I/usr/include
FCFLAGS = -O3 -Wall -march=native -Wunused -fopenmp $(INCLUDES)
LDFLAGS = -lm $(shell nf-config --flibs)

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
