#WIGDIR =/Users/wcoulton/Downloads/Software/wigxjpf-1.9/
WIGDIR =/mnt/raid-cita/meerburg/WignerJ/wigxjpf-1.9/
#all: $(WIGDIR)/lib/libwigxjpf.a $(WIGDIR)/bin/wigxjpf

IFLAG = -I
MPILIB = /opt/openmpi/lib
BINNEDBDIR = /mnt/raid-cita/meerburg/BinnedPrimordialTemplates/BinnedPrimordialTemplates

# The compiler

#When using Intel:
FC     = ifort
#FC      = gfortran
#FFLAGS =  -O3 -ffree-line-length-none -Wno-tabs   -qopenmp
FFLAGS =  -O3 -fast -qopenmp
FCFLAGS += -I $(WIGDIR)/mod/
FCFLAGS += -I $(WIGDIR)/inc/
FCFLAGS += -I $(WIGDIR)/lib/
FCFLAGS += -L $(WIGDIR)/lib/ -lwigxjpf 
FCFLAGS += -L -lBisRoutines
FCFLAGS += -I $(BINNEDBDIR)/
# "make" builds all
default: BCV

SepBispectrum.o:  
theoryBispec.o: SepBispectrum.o
Bispectrum_NGcovariance_v1.o: SepBispectrum.o

#OBJFILES = BisVar_FlatSky.o
#OBJFILES = bisvar_noLmaxLoop_v2.o
OBJFILES = theoryBispec.o SepBispectrum.o Bispectrum_NGcovariance_v1.o
%.o: %.f90
	$(FC) $(FFLAGS) $(FCFLAGS) -c $*.f90 -o $*.o

BCV:  $(OBJFILES)
	$(FC) -o BCV $(OBJFILES) $(FFLAGS) $(FCFLAGS) 

clean: cleantheory

cleantheory:
	rm -f *.o *.mod BCV
