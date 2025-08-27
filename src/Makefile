OBJS = pmodules.o pmain.o plegend.o pread.o pbasis.o psinglepoint.o pdump.o pmachinedependent.o ponee.o pshell.o pfock.o porbitals.o pdensitymatrix.o pfourindexdirect.o pfourindexgradient.o pselfconsistentfield.o pfourindexstore.o pexchangecorrelation.o pgrid.o phellmannfeynman.o pmp2energy.o pdyson.o poptimization.o plebedev.o psinanoglu.o pfci.o pionizationpotential.o pelectronaffinity.o pccpt.o pthermalfci.o pmbgf.o pfullmbgf.o pdiagon.o pcc.o pdegenerate.o ptce.o

#pthreeindexdirect.o pthreeindexstore.o ptwoindex.o pfourindexstoreprimitive.o pselfconsistentdyson.o
#pcis.o pfourindexforcis.o poep.o ppolar.o pcartesiangridhf.o pcartesiangridmp2.o prelativistic.o

# When "unresolved subroutinenames error" occurs in linking stage, 
# just change the orders of the following libraries.

LIBS = lapack/lapack.a blas/blas.a linpack/linpack.a tcedir/tcedir.a

#FC      = gfortran
FC	= mpif90
CFLAGS1 = -c -O2
#CFLAGS1 = -g -c 
LFLAGS1 = -O2
#LFLAGS1 = 

.SUFFIXES: .f90 .f .o

polymer: $(OBJS)
	cd lapack ; make
	cd linpack ; make
	cd blas ; make
	cd tcedir ; make
	ranlib $(LIBS)
	$(FC) $(LFLAGS1) -o polymer $(OBJS) $(LIBS)

.f90.o:
	$(FC) $(CFLAGS1) $*.f90

.f.o:
	$(FC) $(CFLAGS1) $*.f

veryclean:
	rm -f *.o *.mod polymer lapack/*.o linpack/*.o blas/*.o tcedir/*.o lapack/*.a linpack/*.a blas/*.a tcedir/*.a

clean:
	rm -f *.o *.mod polymer
