comp=gnu
ifeq ($(comp),intel)
	FF = ifort -Ofast
	FPP = -DINTEL -fpp
else
ifeq ($(comp),gnu)
	FF = gfortran -Ofast
	FPP = -DGNU -cpp
endif
endif

FFLAGS:=


objects = dogleg.o enorm.o hybrd1.o hybrj1.o qform.o \
	  r1mpyq.o dpmpar.o fdjac1.o hybrd.o hybrj.o \
	  qrfac.o r1updt.o test.o

test : $(objects)
	$(FF) $(FFLAGS) -o test $(objects)

test.o : test.f90
	$(FF) $(FPP) -c test.f90
dogleg.o : dogleg.f
	$(FF) -c dogleg.f
enorm.o : enorm.f
	$(FF) -c enorm.f
hybrd1.o : hybrd1.f
	$(FF) -c hybrd1.f
hybrj1.o : hybrj1.f
	$(FF) -c hybrj1.f
qform.o : qform.f
	$(FF) -c qform.f
r1mpyq.o : r1mpyq.f
	$(FF) -c r1mpyq.f
dpmpar.o : dpmpar.f
	$(FF) -c dpmpar.f
fdjac1.o : fdjac1.f
	$(FF) -c fdjac1.f
hybrd.o : hybrd.f
	$(FF) -c hybrd.f
hybrj.o : hybrj.f
	$(FF) -c hybrj.f
qrfac.o : qrfac.f
	$(FF) -c qrfac.f
r1updt.o : r1updt.f
	$(FF) -c r1updt.f

.PHONY : clean
clean : 
	rm test $(objects)
