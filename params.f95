module params
  !DEFAULT VALUES---------------------------------------------------------------
  !Iterations to run
  integer :: NSTEPS=10000
  integer :: ISTEPS=2000

  !Resolution
  integer :: NX = 512
  integer :: NY = 512
  double precision :: DSPACE = 0.05d0, DTSIZE = 0.001d0

  !Dump frequency - Density - Wavefunction - Misc Utils
  integer :: dumpd = 100, dumpwf = 100, dumputil = 100

  !GPE Type - 0 Natural Units - 1 Hamonic Oscillator Units
  integer :: RHSType = 1
  double precision :: harm_osc_C = 100.0d0
  double precision :: harm_osc_mu = 5.64d0
  complex*16 :: GAMMAC = 0.0d0
  logical :: rtNorm = .false.

  !Boundary Conditions - 0 reflective - 1 periodic
  integer :: BCX = 0
  integer :: BCY = 0

  !Noise Amplitude - 0.001d0 works well
  integer :: noiseamp = 0.000d0

  !Flow Speed in X Dir - Start End Stride
  integer :: VOBS = 0, VOBE = 0, VOBST = 1

  !Potential types - 0 object - 1 free rotational obj - 2 oscillating obj - 3 afm-img
  logical :: enablePot = .true.
  logical :: enableTrap = .true.
  integer :: potType = -1
  !Enable if you need to constantly recalculate the potential
  integer :: potRep = 0
  logical :: doShin = .false.

  !Object properties
  double precision :: RRX=2.0d0
  double precision :: RRY=2.0d0
  double precision :: OBJXDASH=0.0d0
  double precision :: OBJYDASH=0.0d0
  double precision :: OBJHEIGHT=0.0d0
  !Rotating Obj
  double precision :: OBJANGLE=0.1d0
  double precision :: OBJANGLEV=0.1d0
  double precision :: MOMINERTIA=0.000001d0
  !Oscillating Obj
  double precision :: FVAL = 0.0d0, WF = 0.0d0,W0 = 0.05d0 !force amp,force feq, nat freq
  !Trap
  double precision :: TXDASH=0.0d0
  double precision :: TYDASH=0.0d0
  double precision :: TVXDASH=0.0d0
  double precision :: TVYDASH=0.0d0
  double precision :: TTM=0.0d0

  !AFM-IMAGE
  character(2048) :: afm_filename
  integer :: afmRES = 256
  integer :: afmSlice = 120
  double precision :: xi1 = 0.066d0
  double precision :: afmXScale=0.035d0
  double precision :: afmYscale=1.0d0
  double precision :: TRUNCPARAM = 1.0d0
  !GLOBALS----------------------------------------------------------------------
  double precision :: VOB
  integer :: LOOPNO=5
  double precision,parameter :: PI = 4.0d0*ATAN(1.0d0)
  complex*16 :: DT,EYE = (0.0d0,1.0d0)
  double precision :: NORM,OLDNORM = 1.0d0
  complex*16, dimension(:,:), ALLOCATABLE :: GRID,OBJPOT
  double precision :: TIME
  double precision, dimension(2) :: FVECOLD = 0.0d0
  double precision :: OBJYVEL = 0.0d0,OBJXVEL = 0.0d0

contains
  SUBROUTINE init_params
    IMPLICIT NONE
    afm_filename = repeat(" ", 2048) !Clear memory so entire string is blank
    include 'params.in'
    ALLOCATE(GRID(-NX/2:NX/2,-NY/2:NY/2))
    ALLOCATE(OBJPOT(-NX/2:NX/2,-NY/2:NY/2))
   END SUBROUTINE
end module
