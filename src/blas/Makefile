OBJS = dasum.o dgemv.o drotmg.o dswap.o dsyrk.o dtrmv.o daxpy.o dger.o dsbmv.o dsymm.o dtbmv.o dtrsm.o dcopy.o dnrm2.o dscal.o dsymv.o dtbsv.o dtrsv.o ddot.o drot.o dspmv.o dsyr.o dtpmv.o idamax.o dgbmv.o drotg.o dspr.o dsyr2.o dtpsv.o lsame.o dgemm.o drotm.o dspr2.o dsyr2k.o dtrmm.o xerbla.o

#FC = gfortran
FC = mpif90
#FLAGS = -c -i4
#FLAGS = -g -c -pg
FLAGS = -c -O2

.SUFFIXES: .o .f

.f.o:
	$(FC) $(FLAGS) $*.f

blas.a: $(OBJS)
	ar lruv blas.a $(OBJS)
