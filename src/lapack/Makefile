OBJS = dgetrf.o dlaln2.o dlarfg.o dorg2r.o sgehd2.o slacpy.o slanv2.o slascl.o xerbla.o dgebak.o dgetrs.o dlamch.o dlarft.o dorghr.o sgehrd.o sladiv.o slapy2.o slaset.o dgebal.o dhseqr.o dlange.o dlarfx.o dorgqr.o sgesv.o slahqr.o slarf.o  slassq.o dgeev.o dlabad.o dlanhs.o dlartg.o dtrevc.o sgetf2.o slahrd.o slarfb.o slaswp.o dgehd2.o dlacpy.o dlanv2.o dlascl.o ilaenv.o sgetrf.o slaln2.o slarfg.o sorg2r.o dgehrd.o dladiv.o dlapy2.o dlaset.o sgebak.o sgetrs.o slamch.o slarft.o sorghr.o dgesv.o dlahqr.o dlarf.o dlassq.o sgebal.o shseqr.o slange.o slarfx.o sorgqr.o dgetf2.o dlahrd.o dlarfb.o dlaswp.o sgeev.o slabad.o slanhs.o slartg.o strevc.o

#FC = gfortran
FC = mpif90
#FLAGS = -c -i4
#FLAGS = -g -c -pg
FLAGS = -c -O2

.SUFFIXES: .f .o

.f.o:
	$(FC) $(FLAGS) $*.f

lapack.a: $(OBJS)
	ar lruv lapack.a $(OBJS)
