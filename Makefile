WIGDIR =/Users/wcoulton/Downloads/Software/wigxjpf-1.9/
#WIGDIR =/mnt/raid-cita/meerburg/WignerJ/wigxjpf-1.9/
#all: $(WIGDIR)/lib/libwigxjpf.a $(WIGDIR)/bin/wigxjpf

IFLAG = -I
MPILIB = #/opt/openmpi/lib


# The compiler

#When using Intel:
#FC     = ifort
FC      = gfortran
FFLAGS =  -fopenmp -fPIC -O3 -ffree-line-length-none -Wno-tabs 
FCFLAGS += -I $(WIGDIR)/mod/
FCFLAGS += -I $(WIGDIR)/inc/
FCFLAGS += -I /Users/wcoulton/Downloads/Research-Projects/Non-Gaussianity/BinnedPrimordialTemplates/BinnedPrimordialTemplates/
SepBispectrumDir =/Users/wcoulton/Downloads/Research-Projects/Non-Gaussianity/BinnedPrimordialTemplates/BinnedPrimordialTemplates/SepBispectrum.o
FCFLAGS += -I $(WIGDIR)/lib/
FCFLAGS += -L $(WIGDIR)/lib/ -lwigxjpf 

# "make" builds all
default: Bispectrum_NGcovariance_v1

Bispectrum_NGcovariance_v1.o:

OBJFILES = Bispectrum_NGcovariance_v1.o 

%.o: %.f90
	$(FC) $(FFLAGS) $(FCFLAGS) -c $*.f90 -o $*.o

Bispectrum_NGcovariance_v1:  $(OBJFILES)
	$(FC) -o BCV $(OBJFILES) $(SepBispectrumDir) $(FFLAGS) $(FCFLAGS) 

clean: cleantheory

cleantheory:
	rm -f *.o *.mod BCV
