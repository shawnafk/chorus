comp=gnu
ifeq ($(comp),intel)
	FF = ifort -Ofast
	FPP = -DINTEL -fpp
else
ifeq ($(comp),gnu)
	FF = gfortran  -Ofast
	FPP = -DGNU -cpp
endif
endif
FFLAGS:= 
objects = dqage.o dqawoe.o dqcheb.o dqk15w.o dqk41.o dqag.o \
	dqawo.o dqelg.o dqk21.o dqk51.o dqpsrt.o dqc25f.o \
	dqk15.o dqk31.o dqk61.o dqwgtf.o d1mach.o  dgtsl.o \
	dqagp.o dqagpe.o dqaws.o dqawse.o dqc25s.o \
	dqmomo.o dqwgts.o test.o

test : $(objects)
	$(FF) $(FFLAGS) -o test $(objects)

test.o : test.f90
	$(FF) $(FPP) -c test.f90 
dqage.o : dqage.f  
	$(warning  $(FF))
	$(FF) -c dqage.f
dqawoe.o : dqawoe.f 
	$(FF) -c dqawoe.f
dqcheb.o : dqcheb.f 
	$(FF) -c dqcheb.f
dqk15w.o : dqk15w.f 
	$(FF) -c dqk15w.f
dqk41.o : dqk41.f  
	$(FF) -c dqk41.f
dqag.o : dqag.f   
	$(FF) -c dqag.f
dqawo.o : dqawo.f  
	$(FF) -c dqawo.f
dqelg.o : dqelg.f  
	$(FF) -c dqelg.f
dqk21.o : dqk21.f 
	$(FF) -c dqk21.f
dqk51.o :  dqk51.f 
	$(FF) -c dqk51.f
dqpsrt.o : dqpsrt.f
	$(FF) -c dqpsrt.f
dqc25f.o : dqc25f.f
	$(FF) -c dqc25f.f
dqk15.o : dqk15.f
	$(FF) -c dqk15.f
dqk31.o : dqk31.f 
	$(FF) -c dqk31.f
dqk61.o : dqk61.f  
	$(FF) -c dqk61.f
dqwgtf.o : dqwgtf.f
	$(FF) -c dqwgtf.f
d1mach.o : d1mach.f
	$(FF) -c d1mach.f
dgtsl.o : dgtsl.f
	$(FF) -c dgtsl.f
dqagp.o : dqagp.f
	$(FF) -c dqagp.f
dqagpe.o : dqagpe.f
	$(FF) -c dqagpe.f
dqaws.o : dqaws.f
	$(FF) -c dqaws.f
dqawse.o : dqawse.f
	$(FF) -c dqawse.f
dqc25s.o : dqc25s.f
	$(FF) -c dqc25s.f
dqmomo.o : dqmomo.f
	$(FF) -c dqmomo.f
dqwgts.o : dqwgts.f
	$(FF) -c dqwgts.f

.PHONY : clean
clean :
	rm $(objects)
 
