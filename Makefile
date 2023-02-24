comp=gnu
ifeq ($(comp),intel)
FF = ifort 
MPIFF = ifort 
FFLAGS =-shared-intel -mcmodel=medium -qopenmp -Ofast
else

ifeq ($(comp),gnu)
export OMPI_FC=$(FF)
FF = gfortran
MPIFF = mpifort
FFLAGS = -shared-libgcc  -Ofast  -ffree-line-length-512 -fopenmp
endif
endif
LIB = 
minpackDIR := mathlib/minpack
minpack := $(addprefix $(minpackDIR)/,dogleg.o enorm.o \
	hybrd1.o hybrj1.o qform.o \
	r1mpyq.o dpmpar.o fdjac1.o hybrd.o hybrj.o \
	qrfac.o r1updt.o) 

quadpackDIR := mathlib/quadpack
quadpack := $(addprefix $(quadpackDIR)/,dqage.o dqawoe.o \
	dqcheb.o dqk15w.o dqk41.o dqag.o \
	dqawo.o dqelg.o dqk21.o dqk51.o dqpsrt.o dqc25f.o \
	dqk15.o dqk31.o dqk61.o dqwgtf.o d1mach.o  dgtsl.o \
	dqagp.o dqagpe.o dqaws.o dqawse.o dqc25s.o \
	dqmomo.o dqwgts.o)

derivDIR := mathlib/deriv
deriv := $(addprefix $(derivDIR)/, auto_deriv.o)

chorus : $(minpack) $(quadpack) $(deriv) chorus_IO.o chorus_TYPE.o \
	chorus_LIB.o chorus_INIT.o chorus_WAVE.o chorus_EON.o chorus_EON_MOV.o chorus_MAIN.o
	$(MPIFF) $(FFLAGS) -o chorus $(minpack) $(quadpack) $(deriv) \
	chorus_IO.o chorus_TYPE.o chorus_LIB.o chorus_INIT.o chorus_WAVE.o \
	chorus_EON.o chorus_MAIN.o $(LIB)

chorus_IO.o : chorus_IO.f90
	$(FF) $(FFLAGS) -c chorus_IO.f90

chorus_TYPE.o : chorus_TYPE.f90
	$(FF) $(FFLAGS) -c chorus_TYPE.f90

chorus_LIB.o : chorus_LIB.f90
	$(FF) $(FFLAGS) -c chorus_LIB.f90

chorus_INIT.o : chorus_INIT.f90 auto_deriv.f90
	$(FF) $(FFLAGS) -c auto_deriv.f90
	$(FF) $(FFLAGS) -c chorus_INIT.f90

chorus_WAVE.o : chorus_WAVE.f90
	$(FF) $(FFLAGS) -c chorus_WAVE.f90

chorus_EON.o : chorus_EON.f90
	$(FF) $(FFLAGS) -c chorus_EON.f90

chorus_EON_MOV.o : chorus_EON_MOV.f90
	$(FF) $(FFLAGS) -c chorus_EON_MOV.f90

chorus_MAIN.o : chorus_MAIN.mpi.f90
	$(MPIFF) $(FFLAGS) -c chorus_MAIN.mpi.f90 -o chorus_MAIN.o

clean :
	$(RM) *.o *.mod
	$(RM) dist.out
