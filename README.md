# polymer
Gaussian-basis-set <i>ab initio</i> molecular-orbital and density-functional methods for molecules and crystals of one-, two-, and three-dimensional periodicity.

# author
So Hirata (sohirata@illinois.edu)

# how to compile

     $ cd src
     $ make

# how to execute interactively

     $ cd data
     $ debug input nprocs > output &

Modify directory information in "debug"

# how to submit a batch job

     $ cd data
     $ submit input nprocs

Modify directory information in "submitscript"

# input syntax

* JOBNAME (checkpoint file with the same name will be read)
* BASISSET (cc-pVDZ, 6-31Gss, etc.)
* PERIOD (X, Y, Z lattice constants)
* HELIX (helical angle in degrees)
* CHARGE (integer electric charge)
* KPOINTS (X, Y, Z numbers of k points per reciprocal unit cell; 0 0 0 for a molecular calculation)
* JOBTYPE (0 = single point; 1 = gradient; 2 = geometry optimization)
* UNITS (BOHR or ANGSTROM)
* PRINT (0 = default; 1 = verbose; 2 = debug)
* CUTOFF1 (Namur short-range cutoff)
* CUTOFF2 (Namur long-range cutoff)
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
        'OEP       ','MP2KJOBS  ','MP2KJOBE  ','HIGHMP    ','HIGHCIS   ',& !  46 -  50
        'HIGHCIE   ','HIGHCCS   ','HIGHCCE   ','CCDIIS    ','FROZENVIRT',& !  51 -  55
        'HIGHCIROOT','EOMORDERS ','EOMORDERE ','EOMROOT   ','IP        ',& !  56 -  60
        'EA        ','CCCONV    ','CCRELAX   ','CCALG     ','XCC       ',& !  61 -  65
        'CCPTROOT  ','CICONV    ','CCPTCONV  ','SPHERICAL ','CCTHEORY  ',& !  66 -  70
        'OEPALG    ','SLATER51  ','KLI       ','POTDUMP   ','DENDUMP   ',& !  71 -  75
        'KERNEL_HF ','KERNEL_S  ','KERNEL_VWN','KERNEL_B88','KERNEL_LYP',& !  76 -  80
        'RADIAL    ','ANGULAR   ','AC        ','KERNEL_OEP','LINDEP    ',& !  81 -  85
        'POLAR     ','OMEGA     ','POLARALG  ','POLARX    ','POLARY    ',& !  86 -  90
        'POLARZ    ','SOS       ','RELATIVITY','SINANOGLU ','SALG      ',& !  91 -  95
        'SINGULAR  ','HIGHGF    ','TEMP      ','MODULO    ','NGRID     ',& !  96 - 100
        'DELTAH    ','GEMINAL   ','LMAX      ','MBGFALG   ','MP2ALG    ',& ! 101 - 105
        'DIAGONAL  ','DYSONCONV '                                       
     
