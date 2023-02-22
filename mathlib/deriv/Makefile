comp=gnu
ifeq ($(comp),intel)
FF = ifort 
FFLAGS =-shared-intel -mcmodel=medium -qopenmp -Ofast
else
ifeq ($(comp),gnu)
FF = gfortran
FFLAGS = -mcmodel=medium -fopenmp -shared-libgcc -Ofast
endif
endif

auto_deriv.o:
	$(FF) $(FFLAGS) -c auto_deriv.f90

clean:
	rm -f *.o *.mod
