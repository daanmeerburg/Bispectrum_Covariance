WIGDIR =/Users/wcoulton/Downloads/Software/wigxjpf-1.9/
#WIGDIR =/mnt/raid-cita/meerburg/WignerJ/wigxjpf-1.9/
#all: $(WIGDIR)/lib/libwigxjpf.a $(WIGDIR)/bin/wigxjpf

IFLAG = -I
MPILIB = /opt/openmpi/lib


# The compiler

#When using Intel:
#FC     = ifort
FC      = gfortran
FFLAGS =  -O3 -ffree-line-length-none -Wno-tabs   #-qopenmp
FCFLAGS += -I $(WIGDIR)/mod/
FCFLAGS += -I $(WIGDIR)/inc/
FCFLAGS += -I $(WIGDIR)/lib/
FCFLAGS += -L $(WIGDIR)/lib/ -lwigxjpf 

# "make" builds all
default: BCV

bisvar.o:

OBJFILES = bisvar.o

%.o: %.f90
	$(FC) $(FFLAGS) $(FCFLAGS) -c $*.f90 -o $*.o

BCV:  $(OBJFILES)
	$(FC) -o BCV $(OBJFILES) $(FFLAGS) $(FCFLAGS) 

clean: cleantheory

cleantheory:
	rm -f *.o *.mod BCV
