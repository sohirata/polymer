OBJS = chpco.o cqrsl.o dgbdi.o dppfa.o schdc.o spodi.o ssvdc.o zhidi.o zppsl.o cchdc.o chpdi.o csico.o dgbfa.o dppsl.o schdd.o spofa.o strco.o zhifa.o zptsl.o cchdd.o chpfa.o csidi.o dgbsl.o dptsl.o schex.o sposl.o strdi.o zhisl.o zqrdc.o cchex.o chpsl.o csifa.o dgeco.o dqrdc.o schud.o sppco.o strsl.o zhpco.o zqrsl.o cchud.o cpbco.o csisl.o dgedi.o dqrsl.o sgbco.o sppdi.o zchdc.o zhpdi.o zsico.o cgbco.o cpbdi.o cspco.o dgefa.o dsico.o sgbdi.o sppfa.o zchdd.o zhpfa.o zsidi.o cgbdi.o cpbfa.o cspdi.o dgesl.o dsidi.o sgbfa.o sppsl.o zchex.o zhpsl.o zsifa.o cgbfa.o cpbsl.o cspfa.o dgtsl.o dsifa.o sgbsl.o sptsl.o zchud.o zpbco.o zsisl.o cgbsl.o cpoco.o cspsl.o dpbco.o dsisl.o sgeco.o sqrdc.o zgbco.o zpbdi.o zspco.o cgeco.o cpodi.o csvdc.o dpbdi.o dspco.o sgedi.o sqrsl.o zgbdi.o zpbfa.o zspdi.o cgedi.o cpofa.o ctrco.o dpbfa.o dspdi.o sgefa.o ssico.o zgbfa.o zpbsl.o zspfa.o cgefa.o cposl.o ctrdi.o dpbsl.o dspfa.o sgesl.o ssidi.o zgbsl.o zpoco.o zspsl.o cgesl.o cppco.o ctrsl.o dpoco.o dspsl.o sgtsl.o ssifa.o zgeco.o zpodi.o zsvdc.o cgtsl.o cppdi.o dchdc.o dpodi.o dsvdc.o spbco.o ssisl.o zgedi.o zpofa.o ztrco.o chico.o cppfa.o dchdd.o dpofa.o dtrco.o spbdi.o sspco.o zgefa.o zposl.o ztrdi.o chidi.o cppsl.o dchex.o dposl.o dtrdi.o spbfa.o sspdi.o zgesl.o zppco.o ztrsl.o chifa.o cptsl.o dchud.o dppco.o dtrsl.o spbsl.o sspfa.o zgtsl.o zppdi.o chisl.o cqrdc.o dgbco.o dppdi.o spoco.o sspsl.o zhico.o zppfa.o

#FC = gfortran
FC = mpif90
#FLAGS = -c -i4
#FLAGS = -g -c -pg
FLAGS = -c -O2

.SUFFIXES: .f .o

.f.o:
	$(FC) $(FLAGS) $*.f

linpack.a: $(OBJS)
	ar lruv linpack.a $(OBJS)
