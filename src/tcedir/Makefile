OBJS = p_ccd_e.o p_ccd_t2.o p_ccsd_e.o p_ccsd_t1.o p_ccsd_t2.o p_ccsdt_e.o p_ccsdt_t1.o p_ccsdt_t2.o p_ccsdt_t3.o \
       np_ccd_e.o np_ccd_t2.o np_ccsd_e.o np_ccsd_t1.o np_ccsd_t2.o np_ccsdt_e.o np_ccsdt_t1.o np_ccsdt_t2.o np_ccsdt_t3.o \
       p_ccd_t2_disconnected.o p_ccsd_t1_disconnected.o p_ccsd_t2_disconnected.o \
       p_ccsdt_t1_disconnected.o p_ccsdt_t2_disconnected.o p_ccsdt_t3_disconnected.o \
       np_ccd_t2_disconnected.o np_ccsd_t1_disconnected.o np_ccsd_t2_disconnected.o \
       np_ccsdt_t1_disconnected.o np_ccsdt_t2_disconnected.o np_ccsdt_t3_disconnected.o \
       p_ccsd_t2_overlap.o np_ccsd_t2_overlap.o p_ccsdt_t3_overlap.o np_ccsdt_t3_overlap.o \
       p_m_ccsd_e.o p_m_ccsd_t1_disconnected.o p_m_ccsd_t2_disconnected.o p_m_ccsd_t2_overlap.o \
       p_m_ccs_e.o p_m_ccs_t1.o p_m_ccs_t1_disconnected.o \
       p_ccs_e.o np_ccs_e.o p_ccs_t1.o np_ccs_t1.o p_ccs_t1_disconnected.o np_ccs_t1_disconnected.o

#FC = gfortran
FC = mpif90
#FLAGS = -c -i4
#FLAGS = -g -c -pg
FLAGS = -c -O2

.SUFFIXES: .f90 .f .o

.f90.o:
	$(FC) $(FLAGS) $*.f90

.f.o:
	$(FC) $(FLAGS) $*.f

tcedir.a: $(OBJS)
	ar lruv tcedir.a $(OBJS)
