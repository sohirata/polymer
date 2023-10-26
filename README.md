# polymer
Gaussian-basis-set <i>ab initio</i> molecular-orbital and density-functional methods for molecules and crystals of one-, two-, or three-dimensional periodicity (two- and three-dimensional capabilities are untested). For molecules and helical polymers, supported methods are

* HF energy and gradients at zero & nonzero temperature
* DFT energy and gradients at zero & nonzero temperature
* MP2
* GF2
* CCSD
* CIS & CIS-GF2
* TDHF
* TDDFT

For molecules only, additionally

* General-order MP at zero & nonzero temperature
* General-order CI
* General-order CC
* General-order EOM-CC
* General-order IP/EA EOM-CC
* General-order hole/particle CI
* General-order GF
* Grid-based HF
* Grid-based MP2
* OEP & TDOEP
* Dynamic polarizability at TDHF, TDDFT or TDOEP

# author
This software has been developed by So Hirata (sohirata@illinois.edu) with financial support from National Science Foundation and Department of Energy, Office of Science.

# how to compile

Modify MPI Fortran90 compiler and compile options in "Makefile," and then

     $ cd src
     $ make

When pmodules.f90 is modified, 

     $ cd src
     $ make clean
     $ make

# how to execute interactively

Modify directory information in "debug," and then

     $ cd data
     $ debug input nprocs > output &

# how to submit a batch job

Modify directory information in "release" and "releasescript," and then

     $ cd data
     $ release input nprocs

# input keywords

Mandatory input in bold letters.

* <b>GEOMETRY</b> (atomic symbol X, Y, Z coordinates; X 0.0 0.0 0.0 as terminating line)
* JOBNAME (checkpoint file with the same name will be read)
* <b>BASISSET</b> (cc-pVDZ, 6-31Gss, etc.)
* <b>PERIOD</b> (X, Y, Z lattice constants)
* HELIX (helical angle in degrees)
* CHARGE (integer electric charge)
* <b>KPOINTS</b> (X, Y, Z numbers of k points per reciprocal unit cell; 0 0 0 for a molecular calculation)
* JOBTYPE (0 = single point; 1 = gradient; 2 = geometry optimization)
* UNITS (BOHR or ANGSTROM)
* PRINT (0 = default; 1 = verbose; 2 = debug)
* <b>CUTOFF1</b> (Namur short-range cutoff)
* <b>CUTOFF2</b> (Namur long-range cutoff)
* INTTOL (two-electron integrals tolerance)
* SCFCYCLES (maximum number of SCF cycles)
* SCFCONV (SCF convergence criterion)
* SCFRELAX (damping parameter for SCF density matrix)
* DIIS (SCF DIIS order)
* AUXILIARY (auxiliary basis set for RI approximation)
* MULTIPOLE (0 = off; 1 = dipole)
* HFEXCHANGE (HF exchange mixing ratio)
* SLATER (Slater exchange mixing ratio)
* VWN (Vosko-Wilk-Nusair5 correlation mixing ratio)
* BECKE88 (Becke88 exchange mixing ratio)
* LYP (Lee-Yang-Parr correlation mixing ratio)
* MP2 (MP2 correlation mixing ratio)
* DIRECT (T = direct algorithm for two-electron integrals; F = disk-based algorithm)
* MAXDISK (disk space)
* FROZENCORE (number of frozen cores; -1 = 1s cores frozen)
* MAXMEMORY (memory space)
* CUTOFF3 (super-long-range cutoff for asymptotic expansion)
* VAPPROX (RI V algorithm)
* SAPPROX (RI S algorithm)
* DYSON (T = frequency-dependent self-energy; F = frequency-independent self-energy)
* RI_SCF (T = RI SCF)
* RI_MP2 (T = RI MP2)
* THEORY (HF; SNULL; BNULL; SVWN; BLYP; SLYP; BVWN; B3LYP; MP2FC; MP2FULL)
* CIS_ROOTS (number of CIS roots sought)
* RPA_ROOTS (number of TDHF roots sought)
* DAVIDSON (trial vector CI convergence criterion)
* SINGLET (T = singlet excited states sought)
* TRIPLET (T = triplet excited states sought)
* OPTCONV (geometry optimization convergence criterion)
* GDIIS (Geometry DIIS order)
* CISQP (T = CIS-GF2 approximation)
* QPBANDS (number of GF2 bands sought)
* DYSONDAMP (damping parameter for iterative solution of Dyson equation)
* OEP (T = OEP calculation)
* MP2KJOBS (start k point in partial MP2 calculation)
* MP2KJOBE (end k point in partial MP2 calculation)
* HIGHMP (end perturbation order of general-order MP)
* HIGHCIS (start CI order of general-order CI)
* HIGHCIE (end CI order of general-order CI)
* HIGHCCS (start CC order of general-order CC)
* HIGHCCE (end CC order of general-order CC)
* CCDIIS (CC DIIS order)
* FROZENVIRT (number of frozen virtuals)
* HIGHCIROOT (number of excited states sought in general-order CI; >=999 full diagonalization)
* EOMORDERS (start EOM-CC order of general-order EOMCC)
* EOMORDERE (end EOM-CC order of general-order EOMCC)
* EOMROOT (number of excited states sought in general-order EOM-CC; >=999 full diagonalization)
* IP (T = hole-particle CI or IP-EOM-CC)
* EA (T = particle-hole CI or EA-EOM-CC)
* CCCONV (CC convergence criterion)
* CCRELAX (damping parameter for iterative solution of CC)
* CCALG (1 = H exp(T) algorithm; 2 = exp(-T) H exp(T) algorithm)
* XCC (T = expectation-value CC)
* CCPTROOT (number of EOM-CCPT excited states sought)
* CICONV (CI convergence criterion)
* CCPTCONV (EOM-CCPT convergence criterion)
* SPHERICAL (T = spherical d,f; F = Cartesian d,f)
* CCTHEORY (CCD; LCCD; CCSD; LCCSD; QCISD; D123; D45) 
* OEPALG (1 = V; 2 = S; 3 = Slater; 4 = KLI)
* SLATER51 (T = Slater potential calculation)
* KLI (T = Krieger-Li-Iafrate potential calculation)
* POTDUMP (-1 = dump xc potential)
* DENDUMP (-1 = dump density)
* KERNEL_HF (HF exchange kernel mixing ratio)
* KERNEL_S (Slater exchange kernel mixing ratio)
* KERNEL_VWN (Vosko-Wilk-Nusair5 correlation kernel mixing ratio)
* KERNEL_B88 (Becke88 exchange kernel mixing ratio)
* KERNEL_LYP (Lee-Yang-Parr correlation kernel mixing ratio)
* RADIAL (radial grid points per atom)
* ANGULAR (angular Lebedev grid points per atom)
* AC (T = asymptotic correction to the xc potential)
* KERNEL_OEP (OEP exchange kernel mixing ratio)
* LINDEP (basis set linear dependency threshold)
* POLAR (T = dynamic polarizability calculation)
* OMEGA (frequency for dynamic polarizability and GF)
* POLARALG (1 = Pople-type trial vectors in dynamic polarizabiity; 2 = Pople-Davidson-type; 3 = Gauss-Chevyshev Pople-type; 4 = Gauss-Chebyshev Pople-Davidson-type ) 
* POLARX (T = xx polarizability)
* POLARY (T = yy polarizability)
* POLARZ (T = zz polarizability)
* SOS (T = Sum-over-states approximation in dynamic polarizabiity)
* RELATIVITY (deprecated) (T = relativistic full CI)
* SINANOGLU (T = Becke-grid-based MP2 calculation)
* SALG (1 = 7-point linear radial grid in Becke-grid-based MP2; 2 = 2-point linear; 3 = 7-point log; 4 = 2-point log; 5 = 7-point Becke; 6 = 2-point Becke)
* SINGULAR (SVD threshold for OEP)
* HIGHGF (end order of general-order GF)
* TEMP (temperature in K)
* MODULO (integer n in modulo-n approximation)
* NGRID (deprecated) (number of Cartesian grid points in Cartesian grid-based HF & MP2)
* DELTAH (deprecated) (grid spacing in Cartesian grid-based HF & MP2)
* GEMINAL (F12 = Slater-type correlation factor in Becke-grid-based MP2; R12 = linear correlation factor)
* LMAX (angular momentum quantum number of spherical harmonic basis in Becke-grid-based MP2)
* MBGFALG (1 = recursion 1 & $\lambda$ variation; 2 = recursions 1 & 2; 3 = recursions 1 & 2 & 3; -1 = recursion 1 $\omega$ scan; 4 = recursions 1 & 2 & 3 & $\lambda$ variation; 99 = exact; -99 = exact scan; 999 = exact all poles & residues)
* MP2ALG (1 = fast big memory; 2 = medium speed medium memory; 3 = slow small memory)
* DIAGONAL (T = diagonal approximation to self-energy)
* DYSONCONV (GF convergence criterion)                                       
     
