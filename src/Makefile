FCOMP = nvfortran
FCFLAGS = -O2 -acc=gpu -Minfo=all -gpu=cc80
PROGRAM =  mcgrid
PROGRAM_PATH = ../

SRCS =      ran2.f90 \
			iarray.f90 \
			density.f90 \
			gridset.f90 \
			photon.f90 \
			sourceph.f90 \
			stokes.f90 \
			tauint2.f90 \
			mcpolar.f90 
			
OBJECTS = $(SRCS:.f90=.o)

all:    $(PROGRAM) install

$(PROGRAM): $(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^ 
	
%.o:  %.f90
	$(FCOMP)  $(FCFLAGS) -c $<

%.o:  %.f
	$(FCOMP)  $(FCFLAGS) -c $<

install: $(PROGRAM)
	mv mcgrid $(PROGRAM_PATH)

.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD mcgrid
