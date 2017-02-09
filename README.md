### EVSL:  EigenValues Slicing Library (Version 1.0)
              ___      __   __      ___       _    
             | __|     \ \ / /     / __|     | |   
             | _|   _   \ V /   _  \__ \  _  | |__ 
             |___| (_)   \_/   (_) |___/ (_) |____|       
   ChebLanTR, ChebLanNR,  ChebSI, RatLanTr, and RatLanNr
   Polynomial  and   Rational  Filtered  Lanczos   and  subspace
   iteration algorithms For Symmetric Eigenvalue problems

Welcome to EVSL. EVSL is a C library for computing the eigenvalues of
a symmetric matrix  that are located in a given  interval.  This first
release includes the routines listed above and does not yet offer full
parallel implementations  (trivial openMP test programs  are available
in among  the test  drivers).  EVSL also  provides tools  for spectrum
slicing, i.e.,  the technique of  subdividing a given interval  into p
smaller subintervals and computing the eigenvalues in each subinterval
independently.  EVSL  implements a polynomial filtered  Lanczos (thick
restart, no  restart) a rational  filtered Lanczos (thick  restart, no
restart), and a polynomial filtered subspace iteration.

For questions/feedback send e-mail to Yousef Saad [saad@umn.edu]

===
    
### Description of contents:

FILES:
  INC/

    evsl.h           : user-level function prototypes and constant definitions
    blaslapack.h     : C API for BLAS/LAPACK functions used in evsl
    def.h            : miscellaneous macros 
    struct.h         : miscellaneous structs used in evsl
    internal_proto.h : internal function prototypes for SRC/
    
  SRC/
    chelanNr.c    :  Polynomial Filtered no-restart Lanczos
    chelanTr.c    :  Polynomial Filtered thick restart Lanczos
    chebpoly.c    :  Computing and applying polynomial filters
    chebsi.c      :  Polynomial Filtered Subspace iteration
    dumps.c       :  Miscellaneous functions for I/O and for debugging 
    evsl.c        :  Set global variable evslData
    lanbounds.c   :  Lanczos alg. to give a bound of spectrum
    mactime.c     :  Timer for mac os    
    misc_la.c     :  Miscellaneous linear algebra functions
    ratfilter.c   :  Computing and applying rational filters
    ratlanNr.c    :  Rational Filtered no-restart Lanczos
    ratlanTr.c    :  Rational Filtered thick restart Lanczos
    spmat.c       :  Sparse matrix routines
    spslicer.c    :  Spectrum slicing
    suitesparse.c :  Interface to SuiteSparse
    timing.c      :  Timer
    vect.c        :  Vector operations

  liblancheb.a  : library

  TESTS_Gen/         test drivers for Generalized Eigenvalue problem
    LapPLanR_Gen.c : Laplacian matrices. Polynomial filtering T-R Lanczos
    LapRLanN_Gen.c : Laplacian matrices. Rational filtering Lanczos
    lapl.c         : functions to set up a laplacean matrix and also to 
                     compute the exact eigenvalues of Laplacians
    io.c           : parse command-line input parameters

  TEST_Lap/              test drivers for Laplacean matrices
    LapPLanN.c         : Polynomial filtering non-restart Lanczos
    LapPLanN_MatFree.c : "matrix-free" version: not forming matrix but
                         passing mat-vec function
    LapPLanR.c         : Polynomial filtering T-R Lanczos
    LapPSI.c           : Polynomial filtering subspace iterations
    LapRLanN.c         : Rational filtering non-restart Lanczos
    LapRLanR.c         : Rational filtering T-R Lanczos
    lapl.c             : Same as above
    io.c               : Same as above

  TESTS_Mat/     general matrices in sparse format read from a file
    GenPLanN.c : Polynomial filtering non-restart Lanczos
    GenPLanR.c : Polynomial filtering T-R Lanczos
    GenPSI.c   : Polynomial filtering subspace iterations
    GenRLanN.c : Rational filtering non-restart Lanczos
    GenRLanR.c : Rational filtering T-R Lanczos
 
  See below and 00README files therein for further information

-----------------------------------------------------------------------
  INSTALLATION
-----------------------------------------------------------------------

  Library:
     The users only need to modify the file makefile.in 
     [see makefile.in.example for samples of files makefile.in 
     that are given for mac-os and for Linux].
      
     1. cp makefile.in_Linux/MacOS.example makefile.in. 
     2. modify makefile.in 
        [provide BLAS/LAPACK path, and optionally SuiteSparse path]
     3. make clean; make

     make in the directory SRC will create the library liblancheb.a
      
  Test programs:
      In directories TESTS_Lap, TESTS_Gen and TESTS_Mat you will find makefiles to 
      run sample drivers that test a few different situations.

  SuiteSparse:
      SuiteSparse is the default direct linear solver of EVSL, for the
      rational filtering and generalized eigenvalue problem.
      If EVSL is compiled with SuiteSparse and the users do not provide
      function for solving systems with (A-SIGMA*I) and (A-SIGMA*B), 
      EVSL will use UMFPACK to factorize (A-SIGMA*I) or (A-SIGMA*B).
      If EVSL is not compiled with SuitSparse, these functions must be
      provided by the users.
      For generalized eigenvalue problem, the current version of EVSL
      will factor B with CHOLMOD. So, CHOLMOD is required.
 
  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
  NOTE:  SuiteSparse is NOT distributed with EVSL, and is Copyrighted
         by Timothy Davis.  Please refer to that package for its License.
         [http://faculty.cse.tamu.edu/davis/suitesparse.html]
  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

-----------------------------------------------------------------------
  LINKING  WITH  UMFPACK (SuiteSparse 4.5.3)
-----------------------------------------------------------------------
  UMFPACK and CHOLMOD requires AMD, ,  COLAMD, CCOLAMD
  and  CAMD,  and  optionally  METIS 5.1.0.   Compile  each  of  these
  packages  to  have  the  library  file in  the  Lib  directory.   If
  SuiteSparse  is configured  with METIS,  give the  path to  METIS (v
  5.1.0)  as  well   (to  make  libmetis.a,  type   `make  config'  in
  metis-5.1.0/ and then `make') 
  Please  refer to SuiteSparse and METIS for installation details.

-----------------------------------------------------------------------
  RATIONAL FILTERING
-----------------------------------------------------------------------
  Rational  filtering  requires  solving  linear  systems  (where  the
  coefficient matrix  is the  original matrix  shifted by  a *complex*
  shift).  A linear solver routine  must be provided.  

  After  having  computed  the  rational  filter  by  "find_ratf(intv,
  &rat)", users  can call "set_ratf_solfunc(&rat, &Acsr,  NULL, NULL)"
  to tell  EVSL to use UMFPACK  to factorize the shifted  matrices and
  solve the linear systems. The factors will be saved in struct rat.

  EVSL  can  also  take  user-specified  solvers.  To  do  this,  call
  "set_ratf_solfunc(&rat,  &Acsr,  func,  data)" to  pass  the  solver
  functions and the associated data for each pole.  "func" is an array
  of function  pointers of  length num of  poles, i.e.,  rat->num. So,
  func[i]  is the  function  to solve  the systems  with  pole i,  the
  coefficient matrix of  which is $A - s_i I(or, B)$,  where s_i = rat->zk[i]
  is the  complex shift.  "data"  is an array  of (void*) of  the same
  length,  where data[i]  is the  data needed  by func[i].   

  All "func" must be of the following prototype

      void linsolFunc(int n, double *br, double *bz, double *xr, 
                      double *xz, void *data);

  where n  is the size  of the system,  br, bz are  the right-hand
  side (real and  imaginary parts of complex vector),  xr, xz will
  be the  solution (complex vector),  and "data" contains  all the
  data  needed  for  the  solver.    

  Once "set_ratf_solfunc" is done, rational filtering Lanczos methods 
  should be ready to use.
      
-----------------------------------------------------------------------
  MATRIX-FREE SOLVERS
-----------------------------------------------------------------------
  All  the  iterative solvers  in  EVSL  can  be used  in  matrix-free
  ways. Users need only to  provide the matrix-vector product function
  of the following prototype:
  
      void Matvec(double *x, double *y, void *data);

  where y  = A *  x and data  is the pointer  to the associated  data to
  perform the  matvec. The  (void *)  argument is  to provide  a uniform
  interface  to all  user-specific  functions. For  a particular  Matvec
  function, one can pack all data  needed by this function into a struct
  and pass the  pointer of this struct  to EVSL (after cast  it to (void
  *)). This  function needs to  be passed to EVSL  as well, so  EVSL can
  call  this  function  to  perform  all matvecs.  Note  that  when  the
  user-input  matvec  function  is  set,  the CSR  matrix  will  not  be
  referenced even if it  is provided (so it is safe to give NULL for the
  matrix input in this case).

  In  TEST_LAP,  an example  of  matvec  functions for  2D/3D  Laplacian
  matrices is  provided, where the  matrix is not explicitly  formed but
  5pt/7pt stencil is used instead. In  this example, a struct for matvec
  is first defined:

  typedef struct _lapmv_t {  
    int  nx,  ny,  nz; 
    double  *stencil;  
  } lapmv_t;

  and the matvec function is implemented

  void Lap2D3DMatvec(double  *x, double  *y,  void  *data) {  
    lapmv_t *lapmv = (lapmv_t *) data; 
    int nx = lapmv->nx; 
    int ny = lapmv->ny;
    int nz = lapmv->nz; 
    double *stencil = lapmv->stencil; 
    ...  
  }

  in  which   the  pointer  is  first   casted  and  all  the   data  is
  unpacked. Once these are ready, they can be passed to EVSL by calling
    
  SetMatvecFunc(n, &Lap2D3DMatvec, (void*) &lapmv);

  where the first input is the size of the "matrix", the second input is
  the  function pointer  and the  third one  is the  data pointer.  Once
  "SetMatvecFunc"  is  called,  EVSL  will  use  the  registered  matvec
  function to perform all matvecs with A.

  Users should first create a function  wrapper of the above type for an
  external matvec routine. Then, following  the steps in the example, it
  will be straightforward to use it in EVSL.

-----------------------------------------------------------------------
  GENERALIZED EIGENVALUE PROBLEM
-----------------------------------------------------------------------
  For solving A * x = \lambda * B * x, the users must first call
  
    SetRhsMatrix(&Bcsr)

  to pass B as a CSR matrix to EVSL. Internally, EVSL will use CHOLMOD
  to factorize B and save the Cholesky factors.
