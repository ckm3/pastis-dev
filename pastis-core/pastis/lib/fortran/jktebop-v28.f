!=======================================================================
!     PROGRAM JKTEBOP           John Southworth  (jkt~astro.keele.ac.uk)
!                               Astrophysics Group, Keele University, UK
!-----------------------------------------------------------------------
! V(1) = surface brightness ratio      V(15) = third light
! V(2) = sum of fractional radii       V(16) = phase correction
! V(3) = ratio of stellar radii        V(17) = light scaling factor
! V(4) = linear LD for star A          V(18) = integration ring size (o)
! V(5) = linear LD for star B          V(19) = orbital period (days)
! V(6) = orbital inclination           V(20) = ephemeris timebase (days)
! V(7) = e cos(omega) OR ecentricity   V(21) = nonlinear LD for star A
! V(8) = e sin(omega) OR omega         V(22) = nonlinear LD for star B
! V(9) = gravity darkening 1           VEXTRA(1) = primary star radius
! V(10) = gravity darkening 2          VEXTRA(2) = secondary star radius
! V(11) = primary reflected light      VEXTRA(3) = stellar light ratio
! V(12) = secondary reflected light    VEXTRA(4) = eccentricity
! V(13) = stellar mass ratio           VEXTRA(5) = periastron longitude
! V(14) = tidal lead/lag angle (deg)   VEXTRA(6) = reduced chi-squared
! V(23-37) five lots of sine curve [T0,P,amplitude]
! V(38-67) five lots of polynomial [pivot,x,x2,x3,x4,x5]
!-----------------------------------------------------------------------
! Version 1:  Simplex minimisation algorithm and new input / output used
! Version 2:  Monte Carlo simulation and parametr perturbation algorithm
! Version 3:  Adjustment to Monte Carlo LD coeffs and input/output files
! Version 4:  Solves for sum of radii and convergence criterion modified
! Version 5:  Added TASK0 to find LD and GD coeffs.  Minor modifications
! Version 6:  Reflection and  scale factor  can all be fixed or adjusted
! Version 7:  Can use (e,w) or (ecosw,esinw).    SFACT modified to be in
!             magnitudes; observ'nal errors found; spherical star option
! Version 8:  Bootstrapping error analysis added and the output modified
! Version 9:  Command-line arguments allowed, Monte Carlo without param-
!             eter kicking option and fitting for period and Tzero added
! Version 10: Can now use 99999 datapoints. Whole code now in magnitudes
! Version 11: Bug fixes, tasks renumbered,  added sigma clipping and the
!             global fit procedures, but not thoroughly tested these yet
! Version 12: Nonlinear limb darkening law, fitting for times of minimum
!             light, and FMAX corrections included (from Alvaro Gimenez)
! Version 13: Removed  BILINEAR  and modified  TASK1  to just call JKTLD
!             Modified input file and arrays  for the diff types of data
! Version 14: Fixed the requirement for inputtd INTRING to be an integer
!             Fixed formal errors (when observational ones not supplied)
!             Sorted out numerical derivatives problem when  i ~ 90  deg
! Version 15: Added TASK 9 and cubic LD law, modified simulation output,
!             made MRQMIN a factor of 3 faster, and working on red noise
! Version 16: Added ability to include sine perturbations on parameters.
! Version 17: Added ability to specify a third light and its uncertainty
! Version 18: Made possible to fit for r1 and r2 instead of r1+r2 and k.
! Version 19: Added VARY=3 flag + finer fit phasing if r1 or r2 is small
! Version 20: Polynomial, optimisation check, TASK4 debug, DV centralise
! Version 21: Added input of ecosw and esinw observational  constraints.
! Version 22: Added input of  e and omega  as observational constraints.
! Version 23: (did not) fix a minor bug with calculation of third light.
! Version 24: Moved to gfortran compiler.  Adjusted date_and_time usage.
!             Fixed bug with third light. Added ECQUADPHASES subroutine.
! Version 25: Added numerical integration over long exposure times.
! Version 26: Fixed bug with TASK 5 and modified TASK 5 output slightly.
! Version 27: Fixed TASK 8 bug, converted all to real*8, 999999 datapnt.
! Version 28: Corrected ECQUADPHASES. All write(*) or print* => write(6)
! Last modified: 1st March 2012
!-----------------------------------------------------------------------
! Possible modifications in future:
! 1) Extend to WD2003 and WINK
! 2) Port to F90 or F95 to have long lines, modules, and improved output
! 3) Incorporate change of omega (apsidal motion)
! 4) Include radial velocities
! 5) Include a light-time effect
! 6) Allow for multiple light and radial velocity curves
! 7) Try LD power law proposed by Hestroffer (1997A+A...327..199H)
! 8) Add four-parameter LD law (Claret 2000A+A...363.1081C)
!-----------------------------------------------------------------------
! Miscellaneous notes:
! 1) Phase shift has been redefined compared to original EBOP so that it
!    corresponds directly to the phase of primary minimum.
! 2) MRQMIN adjusts coeffs only if VARY (called 'ia' in MRQMIN) is 1
! 3) If VARY=2 then the parameter is fixed during the initial fit but is
!    perturbed by a set amount (flat distribution) for later analyses.
! 4) If VARY(11) and/or  VARY(12) are "-1" then V(11) and/or V(12) are
!    calculated from the system geometry; if 0 they are fixed at the
!    input value and if 1 are freely adjusted to best fit.
! 5) If the mass ratio is <= 0 then both stars are assumed to be spheres
! 6) If ecosw > 5.0 then (ecosw,esinw) will be taken to be  (10+e,omega)
!    and fitting will occur using e and omega as parameters. e and omega
!    can be strongly correlated, but this option is useful if e is known
!    but omega isn't; this can happen for EBs exhibiting apsidal motion.
! 7) Observational errors are looked for in the  input light curve file.
!    If they are not found then  equal weight  is given  to each  point.
! 8) Nonlinear LD is now supported for the two-coefficient  logarithmic,
!    quadratic and square-root laws.  The type of law must be specified.
!    on input.   Star B can also be forced to the same coeffs as star A.
!    BUT: normalisation for logarithmic not possible (not got equations)
! 9) Fitting for times of minimum light is directly possible.  The cycle
!    numbers and times are inputted  on lines  immediately below all the
!     parameter lines in the input file.
! 10) EBOP results are symmetric about 90 degrees for inclination (which
!     means i=89 gives same answer as i=91),  which causes problems with
!     numerical derivativs when i>89.9. In this case some extra is added
!     to the numerical derivative to keep the solution slightly below 90
! 11) If input k is negative then  (r1+r2,k)  is interpreted as  (r1,r2)
! 12) If  VARY=3  then the parameter is optimised during all fits but is
!     not perturbed by a set amount in later analyses (eg. Monte Carlo).
!-----------------------------------------------------------------------
! Task numbers and purposes:
! (1) This outputs LD coefficients for a given Teff, logg, [M/H], Vmicro
! (2) This outputs a model light curve for fixed input parameters.
! (3) This fits a model to an observed light curve and outputs results.
! (4) This fits a model, rejects discrepant observations, and refits.
! (5) This does a pseudo-global minimisation by perturbing input params.
! (6) This investigates how different parameters vary around best fit.
! (7) This conducts bootstrapping simulations to find robust errors.
! (8) This conducts Monte Carlo simulations to find robust errors.
! (9) This conducts residual permutations to deal with correlated noise.
!-----------------------------------------------------------------------
! Language:  JKTEBOP is written in FORTRAN 77, using several extensions
!   to the ANSI standard:   ==   <=   <   >   >=   /=   !   enddo  endif
!   Until version 24 it was only ever compiled in g77.
! g77 compiler: I only occasionally check if this works as the g77 comp-
!   iler is no longer supported.   To compile with g77 you should change
!   the way the DATE_AND_TIME intrinsic function is called.  To do this,
!   simply search for the lines containing "DTTIME*9",  uncomment these,
!   and comment out the lines containing "DTTIME*10".     I successfully
!   compiled JKTEBOP v25 on 29/10/2010 using  gcc version 3.4.6 20060404
!   (Red Hat 3.4.6-4) on a Scientific Linux PC. The compilation command:
!   g77 -O -Wuninitialized -fbounds-check -fno-automatic -o jktebop
! g95 compiler: this compiles successfully but I have not actually tried
!   to run the resulting executable file.    Compiler version used: "G95
!   (GCC 4.1.2 (g95 0.93!) Jun 16 2010)"  running on a kubuntu 10.04 PC.
! gfortran compiler:  this compiles successfully and executes correctly.
!   The compiler version used last time I modified the current text was:
!   "GNU Fortran (Ubuntu 4.4.3-4ubuntu5) 4.4.3" running on kubuntu 10.04
! Intel-Fortran compiler: this is periodically checked and found to work
!   well. JKTEBOP versions v26 and earlier must be compiled with the -r8
!   command-line flag in order to avoid numerical noise arising from the
!   use of single-precision variables. JKTEBOP v27 onwards is all real*8
!=======================================================================
!=======================================================================
      PROGRAM JKTEBOP
      implicit none
      integer VERSION               ! Version of this code
      integer TASK                  ! Task to perform (between 0 and 5)
      character INFILE*30           ! Name of the input parameter file
      real*8  V(67)                 ! Fundamental EBOP model parameters
      integer VARY(67)              ! Params fixed (0) or variable (1)
      integer LDTYPE(2)             ! Type of LD law for each star
      real*8  DATA(3,999999)        ! Data: (time,magnitude,uncertainty)
      integer DTYPE(999999)         ! Type of data (between 1 and 4)
      integer NDATA                 ! Number of data points
      integer NMIN                  ! Number of times of minimum light
      integer NLR                   ! Number of observed light ratios
      integer NL3                   ! Number of obsd third light values
      integer NECW                  ! Number of observed e*cos(omega)'s
      integer NESW                  ! Number of observed e*sin(omega)'s
      integer NSINE                 ! Number of sine curves to include
      integer PSINE(5)              ! Which parameter each sine acts on
      integer NPOLY                 ! Number of polynomials to include
      integer PPOLY(5)              ! Which parameter each poly acts on
      integer NSIM                  ! Number of simulations or refits
      real*8  SIGMA                 ! Sigma rejection value for task 4
      integer NUMINT                ! Number of numerical integrations
      real*8  NINTERVAL             ! Time interval for numerical ints.
      integer i,ERROR               ! Loop counter and error flag
      real*8 getmin
            ! DTYPE=1 for light curve datapoint
            ! DTYPE=2 for light ratio
            ! DTYPE=3 for times of minimum light
            ! DTYPE=4 for third light
            ! DTYPE=5 for e*cos(omega) or e
            ! DTYPE=6 for e*sin(omega) or omega

      VERSION = 28

      ERROR = 0
      do i = 1,999999
        DATA(1,i) = 0.0d0
        DATA(2,i) = 0.0d0
        DATA(3,i) = 0.0d0
        DTYPE(i) = 0
      end do
      do i = 1,67
        V(i) = -100.0d0
        VARY(i) = 0
      end do
      NDATA = 0
      NMIN = 0
      NLR = 0
      NL3 = 0
      NECW = 0
      NESW = 0
      NSIM = 0
      NSINE = 0
      NPOLY = 0
      do i = 1,5
        PSINE(i) = 0
        PPOLY(i) = 0
      end do
      LDTYPE(1) = -100
      LDTYPE(2) = -100
      NUMINT = 1
      NINTERVAL = 0.0d0

      V(14) = 0.0d0                 ! Tidal angle (leave at 0 normally)
      VARY(14) = 0                  ! Don't vary tidal lead/lag angle
      VARY(18) = 0                  ! Don't vary integration ring size


            ! First output the code name and version, then check for
            ! command-line arguments. If one is found, then take it to
            ! be the input file name.

      write(6,*) " "
      write(6,'(A10,I2,A13,A55)') "JKTEBOP  v",VERSION,"      John So",
     &        "uthworth  (Keele University, UK, jkt~astro.keele.ac.uk)"

      if ( iargc() == 1 ) then
        CALL GETARG (1,INFILE)
      else
        write(6,'(A39,A41)') "A package for modelling the light curve",
     &                      "s of well-detached eclipsing binary stars"
        write(6,'(A39,A41)') "Task 1  outputs limb and gravity darken",
     &                      "ing coefficients for given Teff and log g"
        write(6,'(A39,A41)') "Task 2  outputs one model light curve c",
     &                      "alculated using a set of input parameters"
        write(6,'(A39,A41)') "Task 3  finds the best fit of the model",
     &                      " to observations (internal errors quoted)"
        write(6,'(A39,A41)') "Task 4  finds the best fit to the obser",
     &                      "vations, sigma clips, and refits the data"
        write(6,'(A39,A41)') "Task 5  finds global best fit, by pertu",
     &                      "rbing parameters and refitting many times"
        write(6,'(A39,A41)') "Task 6  fits observations and finds goo",
     &                      "dness of fit for several parameter values"
        write(6,'(A39,A41)') "Task 7  finds robust reliable errors by",
     &                      " analysing with a bootstrapping algorithm"
        write(6,'(A39,A41)') "Task 8  finds robust errors by analysin",
     &                      "g with a Monte Carlo simulation algorithm"
        write(6,'(A39,A41)') "Task 9  finds robust errors using Monte",
     &                      " Carlo given significant correlated noise"
        write(6,*) " "
        write(6,'(A39,A41)') "Usage:  'jktebop  [inputfile]' to under",
     &                      "take one of these tasks                  "
        write(6,'(A39,A41)') "Usage:  'jktebop   newfile'    to outpu",
     &                      "t an empty input file for a given task   "
        write(6,'(A39,A41)') "Usage:  'jktebop     1'        to under",
     &                      "take Task 1 (limb darkening coefficients)"
        write(6,*) " "
        stop
      end if

      if ( INFILE == "newfile" .or. INFILE == "NEWFILE" ) then
        CALL NEWFILE ()
        write(6,*) " "
        stop
      else if ( INFILE == "1" ) then
        CALL TASK1 ()
        STOP
      end if

      CALL INPUT (VERSION,TASK,INFILE,V,VARY,LDTYPE,DATA,DTYPE,NDATA,
     &            NLR,NMIN,NSIM,SIGMA,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,
     &            NESW,NUMINT,NINTERVAL,ERROR)
      if ( ERROR /= 0 ) then
        write(6,*) " "
        stop
      end if

      if ( TASK == 2 ) CALL TASK2 (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY)

      if ( TASK == 3 .or. TASK == 4 ) CALL TASK34 (TASK,V,VARY,LDTYPE,
     &                      DATA,DTYPE,NDATA,NLR,NMIN,SIGMA,NSINE,PSINE,
     &                       NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      if ( TASK == 6 ) CALL TASK6 (V,VARY,LDTYPE,DATA,DTYPE,NDATA,NLR,
     &      NMIN,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      if ( TASK == 5 .or. TASK == 7 .or. TASK == 8 .or. TASK == 9)
     &  CALL TASK5789 (TASK,V,VARY,LDTYPE,DATA,DTYPE,NDATA,NLR,NMIN,
     &      NSIM,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      write(6,*) " "

      END PROGRAM JKTEBOP
!=======================================================================
!=======================================================================
      SUBROUTINE NEWFILE ()                ! Makes all empty input files
      implicit none
      integer TASK                         ! Task the file is wanted for
      character*30 OUTFILE                 ! Name of output file to make
      integer i,ERROR                      ! Loop counter and error flag

      ERROR = 0
      write(6,'(A39,$)')   "Enter name of input file to create >>  "
      read (*,*) OUTFILE
      write(6,'(A39,$)')   "Enter task number (between 2 and 9) >> "
      read (*,*,iostat=ERROR) TASK

      if ( ERROR /= 0 ) then
        write(6,'(A40)') "### ERROR: did not understand input.    "
        return
      else if ( TASK < 2 .or. TASK > 9 ) then
        write(6,'(A40)') "### ERROR: task integer is out of range."
        return
      else
        CALL OPENFILE (50,"new"," output   ",OUTFILE,ERROR)
        if ( ERROR /= 0 ) return
      end if

      write(50,100)"Task to do (from 2 to 9)   Integ. ring size (deg)  "
      write(50,100)"Sum of the radii           Ratio of the radii      "
      write(50,100)"Orbital inclination (deg)  Mass ratio of system    "
      write(50,100)"Eccentricity or ecosw      Peri.longitude or esinw "
      write(50,100)"Gravity darkening (star A) Grav darkening (star B) "
      write(50,100)"Surface brightness ratio   Amount of third light   "
      write(50,100)"LD law type for star A     LD law type for star B  "
      write(50,100)"LD star A (linear coeff)   LD star B (linear coeff)"
      write(50,100)"LD star A (nonlin coeff)   LD star B (nonlin coeff)"
      write(50,100)"Reflection effect star A   Reflection effect star B"
      write(50,100)"Phase of primary eclipse   Light scale factor (mag)"

      if ( TASK >= 3 .and. TASK <= 9 ) then
        write(50,101)"Orbital period of eclipsing binary system (days) "
        write(50,101)"Reference time of primary minimum (HJD)          "
      end if

      if ( TASK == 2 ) then
        write(50,101)"Output file name (continuous character string)   "
      else if ( TASK == 4 ) then
        write(50,101)"Sigma value to reject discrepant observations    "
      else if ( TASK == 5 ) then
        write(50,101)"Number of refits for perturbed initial parameters"
      else if ( TASK == 7 ) then
        write(50,101)"Number of bootstrapping simulations to do        "
      else if ( TASK == 8 ) then
        write(50,101)"Number of Monte Carlo simulations to do          "
      end if

      if ( TASK >= 3 .and. TASK <= 9 ) then
      write(50,103)"Adjust RADII SUM  or  RADII RATIO      (0, 1, 2, 3)"
      write(50,103)"Adjust INCLINATION  or  MASSRATIO      (0, 1, 2, 3)"
      write(50,103)"Adjust ECCENTRICITY  or  OMEGA         (0, 1, 2, 3)"
      write(50,103)"Adjust GRAVDARK1  or  GRAVDARK2        (0, 1, 2, 3)"
      write(50,103)"Adjust SURFACEBRIGHT2  or  THIRDLIGHT  (0, 1, 2, 3)"
      write(50,103)"Adjust LD-lin1  or  LD-lin2            (0, 1, 2, 3)"
      write(50,103)"Adjust LD-nonlin1  or  LD-nonlin2      (0, 1, 2, 3)"
      write(50,103)"Adjust REFLECTION COEFFS 1 and 2       (-1,0,1,2,3)"
      write(50,103)"Adjust PHASESHIFT  or  SCALE FACTOR    (0, 1, 2, 3)"
      write(50,103)"Adjust PERIOD  or  TZERO (min light)   (0, 1, 2, 3)"
      write(50,100)"Name of file containing light curve                "
      write(50,100)"Name of output parameter file                      "
      end if

      if ( TASK == 3 .or. TASK == 4 ) then
        write(50,101)"Name of output light curve file                  "
        write(50,101)"Name of output model light curve fit file        "
      end if

      if ( TASK == 5 .or. TASK == 7 .or. TASK == 8 .or. TASK == 9 ) then
        write(50,101)"Name of output file of individual fits           "
      end if

      write (50,*) " "
      write (50,*) " "
      write (50,102) "# Enter the appropriate numbers on the l",
     &               "eft-hand side of each line of this file."
      write (50,102) "# Most of the lines require two numerica",
     &               "l parameters separated by spaces.       "
      write (50,*) " "
      write (50,102) "# Put a negative number for the mass rat",
     &               "io to force the stars to be spherical.  "
      write (50,102) "# The mass ratio will then be irrelevant",
     &               " (it is only used to get deformations). "
      write (50,*) " "
      write (50,102) "# To fit for r1 and r2 instead of (r1+r2",
     &               ") and k, give a negative value for r1+r2"
      write (50,102) "# Then (r1+r2) will be interpreted to me",
     &               "an r1, and k will be interpreted as -r2."
      write (50,102) "# The adjustment indicators will similar",
     &               "ly refer to r1,r2 rather than (r1+r2),k."
      write (50,*) " "
      write (50,102) "# If e < 10 then e and omega will be ass",
     &               "umed to be e*cos(omega) and e*sin(omega)"
      write (50,102) "# If e >= 10 then e and omega will be as",
     &               "sumed to be (e+10) and omega (degrees). "
      write (50,102) "# The first option is in general better ",
     &               "unless eccentricity is larger or fixed. "
      write (50,*) " "
      write (50,102) "# The possible entries for the type of l",
     &               "imb darkening law are 'lin' (for linear)"
      write (50,102) "# 'log' (logarithmic), 'sqrt' (square-ro",
     &               "ot), 'quad' (quadratic) or 'cub' (cubic)"
      write (50,102) "# Put 'same' for star B to force its coe",
     &               "fficients to be equal those of star A.  "
      write (50,*) " "

      if ( TASK >= 3 .and. TASK <= 9 ) then
        write (50,*) " "
        write (50,102) "# For each adjustable parameter the adju",
     &                 "stment integer can be 0  (parameter will"
        write (50,102) "# be fixed at the input file value),  1 ",
     &                 "(parameter will be freely adjusted),  2 "
        write (50,102) "# (parameter will be fixed for initial f",
     &                 "it but perturbed during later analysis)."
        write (50,102) "# or  3 (adjusted in initial fit but not",
     &                 " perturbed during Monte Carlo analysis)."
        write (50,*) " "
        write (50,102) "# When fitting a light curve  the reflec",
     &                 "tion coefficients can be calculated from"
        write (50,102) "# the system geometry  (put -1 for the a",
     &                 "djustment integers),  held fixed (put 0)"
        write (50,102) "# or freely adjusted to fit the light cu",
     &                 "rve (put 1) - useful for close binaries."
        write (50,*) " "
      end if

      if ( TASK == 7 .or. TASK == 8 ) then
        write (50,*) " "
        write (50,102) "# When doing Monte Carlo or bootstrappin",
     &                 "g the starting parameters for each simu-"
        write (50,102) "# lation are perturbed by a pre-defined ",
     &                 "amount to avoid biassing the results. If"
        write (50,102) "# this is not wanted, put a minus sign i",
     &                 "n front of the number of simulations.   "
        write (50,*) " "
      end if

      write (50,*) " "
      write (50,102) "# TIMES OF MINIMUM LIGHT: add a line bel",
     &               "ow the parameter line to input each one:"
      write (50,102) "#   'TMIN  [cycle]  [time]  [error]'    ",
     &               "                                        "
      write (50,102) "# where [cycle] is cycle number (integer",
     &               " for primary minimum  or integer+0.5 for"
      write (50,102) "# secondary minimum), [time] and [error]",
     &               " are the observed time and uncertainty. "
      write (50,*) " "
      write (50,102) "# LIGHT RATIO: add a line below the para",
     &               "meter line to input each observed one:  "
      write (50,102) "#   'LRAT'  [time]  [light_ratio]  [erro",
     &               "r]                                      "
      write (50,102) "# where [time] is the time(HJD) when the",
     &               " spectroscopic light ratio was measured,"
      write (50,102) "# [light_ratio] is its value and [error]",
     &               " is its measurement uncertainty.        "
      write (50,*) " "
      write (50,102) "# MEASURED THIRD LIGHT VALUE:   include ",
     &               "as observed constraint by adding a line:"
      write (50,102) "# 'THDL'  [value]  [uncertainty]        "
      write (50,102) "# which gives the third light measuremen",
     &               "t and its observational uncertainty.    "
      write (50,*) " "
      write (50,102) "# MEASURED orbital shape parameters (dep",
     &               "ending on the value of eccentricity):   "
      write (50,102) "#  ECSW  [value]  [uncertainty]    (inte",
     &               "rpreted as either e*cos(omega) or e)    "
      write (50,102) "#  ENSW  [value]  [uncertainty]    (inte",
     &               "rpreted as either e*sin(omega) or omega)"
      write (50,*) " "
      write (50,102) "# SINE AND POLYNOMIAL FITTING:   the par",
     &               "ameters of sine curves or polynomials of"
      write (50,102) "# order 5) can be included.  You can hav",
     &               "e up to five sines and five polynomials,"
      write (50,102) "# each acting on a specific parameter. T",
     &               "he information for each one is specified"
      write (50,102) "# by an additional line below the main i",
     &               "nput file parameters.     Line format:  "
      write (50,102) "#   SINE  [par]  [T0]  [P]  [amp]  [vary",
     &               "(T0)]  [vary(P)]  [vary(amp)]           "
      write (50,102) "#   POLY  [par] [pivot] [x] [x^2] [x^3] ",
     &               "[x^4] [x^5] [vary(const)] [vary(x)] ...."
      write (50,102) "# where the required parameters are give",
     &               "n inside square brackets. [T0] is a time"
      write (50,102) "# of zero phase (HJD), [P] is period (da",
     &               "ys), [amp] is amplitude, [const] is a   "
      write (50,102) "# constant and [x^n] are the coefficient",
     &               "s of the polynomial.  Each parameter has"
      write (50,102) "# a [vary()] which is 0, 1, 2 or 3 to in",
     &               "dicate how the parameter is treated.    "
      write (50,102) "# [par] indicates what to apply it to: c",
     &               "hoose between:  J r1 r2 i L3 sf L1 L2   "
      write (50,102) "# Note that the independent parameter is",
     &               " always time (either HJD or phase).     "
      write (50,*) " "
      write (50,102) "# NUMERICAL INTEGRATION:  long exposure ",
     &               "times can be split up into NUMINT points"
      write (50,102) "# occupying a total time interval of NIN",
     &               "TERVAL (seconds) by including this line:"
      write (50,102) "#   NUMI  [NUMINT]  [NINTERVAL]         ",
     &               "                                        "
      write (50,*) " "
      write (50,*) " "


      close (50,status="keep")
      write(6,'(A29,I1,A21,A30)') "An empty input file for task ",TASK,
     &                                   " has been written to ",OUTFILE

100   FORMAT (18X,A51)
101   FORMAT (18X,A49)
102   FORMAT (A40,A40)
103   FORMAT (" 0  0",13X,A51)

      END SUBROUTINE NEWFILE
!=======================================================================
!=======================================================================
      SUBROUTINE INPUT (VERSION,TASK,INFILE,V,VARY,LDTYPE,DATA,DTYPE,
     &                  NDATA,NLR,NMIN,NSIM,SIGMA,NSINE,PSINE,NPOLY,
     &                  PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL,STATUS)
      implicit none
      integer VERSION,TASK                ! IN:  Code version Task numbr
      character INFILE*30                 ! IN:  Name of the  input file
      real*8 V(67)                        ! OUT: Photometric  parameters
      integer VARY(67)                    ! OUT: Par adjustment ingeters
      integer LDTYPE(2)                   ! OUT: Type of LD law to adopt
      real*8 DATA(3,999999)               ! OUT: Data (time, mag, error)
      integer DTYPE(999999)               ! OUT: Type of data (1 2 or 3)
      integer NDATA,NLR,NMIN              ! OUT: Amount of types of data
      integer NSIM                        ! OUT: Number  of  simulations
      real*8 SIGMA                        ! OUT: Sigma  rejection  value
      integer NSINE,NL3,NECW,NESW         ! OUT: Numbrs of sines and L3s
      integer PSINE(5)                    ! OUT: Which par for each sine
      integer NPOLY,PPOLY(5)              ! OUT: Similar for polynomials
      integer NUMINT                      ! OUT: Number of numerical int
      real*8  NINTERVAL                   ! OUT: Time interval of numint
      character OBSFILE*30                ! LOCAL: Input lightcurve file
      character LCFILE*30                 ! LOCAL: Output lightcurv file
      character PARFILE*30                ! LOCAL: Output parameter file
      character FITFILE*30                ! LOCAL: Output model fit file
      integer i,j,k,ERROR,STATUS          ! LOCAL: Counter & error flags
      real*8 LP,LS                        ! LOCAL: Light from  each star
      real*8 ARRAY(999999),SELLECT        ! LOCAL: Help array & FUNCTION
      character DTDATE*8,DTTIME*10        ! LOCAL: runtime time and date
!!      character DTDATE*8,DTTIME*9    ! form needed for g77 compilation
      character LD1*4,LD2*4               ! LOCAL: type of LD law to use
      character CHARHELP200*200           ! LOCAL: variable to read in
      character CHARHELP4*4               ! LOCAL: variable to read in
      real*8 GETMODEL,MAG
      character WHICHPAR*2
      integer IARRAY(6)
      real*8 R1,R2

      ERROR = 0
      STATUS = 0
      CALL OPENFILE (60,"old","input file",INFILE,ERROR)
      if ( ERROR /= 0 ) then
        write(6,*) " "
        STOP
      end if

!-----------------------------------------------------------------------

      read (60,*,iostat=ERROR) TASK,V(18)
      if ( ERROR /= 0 ) then
        write(6,'(A39,A41)') "### ERROR: cannot read first line of in",
     &                     "put parameter file (TASK and INTRING).   "
        STATUS = 1
        write(6,'(A28,I5)') "### Error message returned: ",ERROR
      else if ( TASK < 0 .or. TASK > 9 ) then
        write(6,'(A39,A41)') "### ERROR: task integer does not corres",
     &                     "pond to a value between 1 and 9.         "
        STATUS = 1
      else if ( V(18) < 0.01d0 .or. V(18) > 10.0d0 ) then
        write(6,'(A39,A41)') "### ERROR: integration ring size must b",
     &                     "e between 0.01 and 10.0 degrees.         "
        STATUS = 1
      end if
      if ( STATUS /= 0 ) return

      if ( TASK == 1 ) write(6,'(A25,A55)')"Task 1  outputs limb and ",
     &        "gravity darkening coefficients for given Teff and log g"
      if ( TASK == 2 ) write(6,'(A25,A55)')"Task 2  outputs one model",
     &        " light curve calculated using a set of input parameters"
      if ( TASK == 3 ) write(6,'(A25,A55)')"Task 3  finds the best fi",
     &        "t of the model to observations (internal errors quoted)"
      if ( TASK == 4 ) write(6,'(A25,A55)')"Task 4  finds the best fi",
     &        "t to the observations, sigma clips, and refits the data"
      if ( TASK == 5 ) write(6,'(A25,A55)')"Task 5  finds global best",
     &        " fit, by perturbing parameters and refitting many times"
      if ( TASK == 6 ) write(6,'(A25,A55)')"Task 6  fits observations",
     &        " and finds goodness of fit for several parameter values"
      if ( TASK == 7 ) write(6,'(A25,A55)')"Task 7  finds robust reli",
     &        "able errors by analysing with a bootstrapping algorithm"
      if ( TASK == 8 ) write(6,'(A25,A55)')"Task 8  finds robust erro",
     &        "rs by analysing with a Monte Carlo simulation algorithm"
      if ( TASK == 9 ) write(6,'(A25,A55)')"Task 9  finds robust erro",
     &        "rs using Monte Carlo given significant correlated noise"

      CALL READFF (60,"RADII SUM ",-0.8d0, 0.8d0, V( 2),
     &                "RADIIRATIO", 0.0d2, 1.0d2, V( 3),STATUS)
      CALL READFF (60,"INCLNATION",50.0d0, 1.4d2, V( 6),
     &                "MASS RATIO",-1.0d9, 1.0d3, V(13),STATUS)
      CALL READFF (60,"ECCENTRCTY",-1.0d0,11.0d0, V( 7),
     &                "OMEGA     ",-3.6d2, 3.6d2, V( 8),STATUS)
      CALL READFF (60,"GRAVDARK-1",-1.0d1, 1.0d1, V( 9),
     &                "GRAVDARK-2",-1.0d1, 1.0d1, V(10),STATUS)
      CALL READFF (60,"SURF-BRT-2", 0.0d0, 1.0d2, V( 1),
     &                "THIRDLIGHT", 0.0d0, 1.0d1, V(15),STATUS)

      read (60,*) LD1,LD2
      if ( LD1 == "lin"  ) LDTYPE(1) = 1
      if ( LD1 == "log"  ) LDTYPE(1) = 2
      if ( LD1 == "sqrt" ) LDTYPE(1) = 3
      if ( LD1 == "quad" ) LDTYPE(1) = 4
      if ( LD1 == "cub"  ) LDTYPE(1) = 5
      if ( LD2 == "same" ) LDTYPE(2) = 0
      if ( LD2 == "lin"  ) LDTYPE(2) = 1
      if ( LD2 == "log"  ) LDTYPE(2) = 2
      if ( LD2 == "sqrt" ) LDTYPE(2) = 3
      if ( LD2 == "quad" ) LDTYPE(2) = 4
      if ( LD2 == "cub"  ) LDTYPE(2) = 5
      if ( LDTYPE(1) <= 0 ) then
        write(6,'(A39,A41)') "### ERROR: LD law type for star A should",
     &                   " be 'lin', 'log', 'sqrt', 'quad', or 'cub'. "
        STATUS = 1
      end if
      if ( LDTYPE(2) < 0 )then
        write(6,'(A39,A41)') "### ERROR: LD law for star B must be 'li",
     &                   "n', 'log', 'sqrt', 'quad', 'cub', or 'same'."
        STATUS = 1
      end if
      if ( STATUS /= 0 ) return

      CALL READFF (60,"LD-lin-1  ",-1.0d0, 2.0d0, V( 4),
     &                "LD-lin-2  ",-1.0d0, 2.0d0, V( 5),STATUS)
      if ( LDTYPE(2) == 0 ) V(5) = V(4)
      CALL READFF (60,"LDnonlin-1",-1.0d0, 2.0d0, V(21),
     &                "LDnonlin-2",-1.0d0, 2.0d0, V(22),STATUS)
      if ( LDTYPE(2) == 0 ) V(22) = V(21)
      CALL READFF (60,"REFLECTN-1", 0.0d0, 1.0d0, V(11),
     &                "REFLECTN-2", 0.0d0, 1.0d0, V(12),STATUS)
      CALL READFF (60,"PHASESHIFT",-1.0d0, 1.0d0, V(16),
     &                "SCALEFACTR",-1.0d3, 1.0d3, V(17),STATUS)

      if ( STATUS /= 0 ) return

      if ( TASK >= 3 .and. TASK <= 9 ) then
        CALL READF (60, "PERIOD    ", 0.0d0, 1.0d6, V(19),STATUS)
        CALL READF (60, "TIME-ZERO ",-1.0d4, 3.0d6, V(20),STATUS)
      end if

      if ( TASK == 2 ) then
        CALL READCHAR30 (60,"OUTFILE   ",PARFILE,STATUS)
        close (60)
        CALL OPENFILE (62,"new","lightcurve",PARFILE,STATUS)
        return
      end if

      if ( STATUS /= 0 ) return

      if ( TASK == 5 .or. TASK == 7 .or. TASK == 8 ) then
        read (60,*,iostat=ERROR) NSIM
        if ( ERROR /= 0 ) then
          if (TASK==5) write(6,'(A24,A56)') "### ERROR reading number",
     &       " of perturbed initial parameter sets to refit data with."
          if (TASK==7) write(6,'(A24,A56)') "### ERROR reading the nu",
     &       "mber of bootstrapping simulations to do.                "
          if (TASK==8) write(6,'(A24,A56)') "### ERROR reading the nu",
     &       "mber of Monte Carlo simulations to do.                  "
          STATUS = 1
        end if
        if ( abs(NSIM) < 8 .or. abs(NSIM) > 99999 ) then
          write(6,'(A37,A43)') "### ERROR: number of simulations/sets",
     &                    " to do must be between +/-8 and +/-99999.  "
          STATUS = 1
        end if
      else if ( TASK == 4 ) then
        read (60,*,iostat=ERROR) SIGMA
        if ( ERROR /= 0 ) then
          write(6,'(A37,A43)') "### ERROR reading the sigma number fo",
     &                    "r rejection of discrepant data.            "
          STATUS = 1
        end if
      end if

      if ( STATUS /= 0 ) return

      if ( TASK >= 3 .and. TASK <= 9 ) then
        CALLREAD2(60,"adj(R1+R2)","adj(R2/R1)",VARY( 2),VARY( 3),STATUS)
        CALLREAD2(60,"adj(_INC_)","adj(M2/M1)",VARY( 6),VARY(13),STATUS)
        CALLREAD2(60,"adj(_ECC_)","adj(OMEGA)",VARY( 7),VARY( 8),STATUS)
        CALLREAD2(60,"adj(_GD1_)","adj(_GD2_)",VARY( 9),VARY(10),STATUS)
        CALLREAD2(60,"adj(_SB2_)","adj(_L_3_)",VARY( 1),VARY(15),STATUS)
        CALLREAD2(60,"adj(LD-u1)","adj(LD-u2)",VARY( 4),VARY( 5),STATUS)
        CALLREAD2(60,"adj(LD-V1)","adj(LD-V2)",VARY(21),VARY(22),STATUS)
        CALLREAD2(60,"adj(REFL1)","adj(REFL2)",VARY(11),VARY(12),STATUS)
        CALLREAD2(60,"adj(PSHFT)","adj(SFACT)",VARY(16),VARY(17),STATUS)
        CALLREAD2(60,"adj(PERIOD","adj(TZERO)",VARY(19),VARY(20),STATUS)
        CALL READCHAR30 (60,"OBSFILE   ",OBSFILE,STATUS)
        CALL READCHAR30 (60,"PARAMFILE ",PARFILE,STATUS)
        if ( STATUS /= 0 ) return
      end if

      if ( LDTYPE(1) == 1 ) VARY(21) = 0
      if ( LDTYPE(2) == 1 ) VARY(22) = 0
      if ( LDTYPE(2) == 0 ) VARY(5) = 0
      if ( LDTYPE(2) == 0 ) VARY(22) = 0

      if ( VARY(20) == 1 .and. VARY(16) == 1 ) then
        write(6,'(A39,A41)') ">> TZERO and PSHIFT cannot both be adju",
     &                      "sted so adj(PSHIFT) has been set to zero."
        VARY(16) = 0
      end if

      if(TASK==3.or.TASK==4.or.TASK==5.or.TASK==7.or.TASK==8.or.TASK==9)
     &                  CALL READCHAR30 (60,"LC FILE   ",LCFILE, STATUS)
      if ( TASK == 3 .or. TASK==4)
     &                  CALL READCHAR30 (60,"FITFILE   ",FITFILE,STATUS)


!-----------------------------------------------------------------------

            ! Open the output files and start writing results to them.

      if (TASK >= 3) CALL OPENFILE(62,"new","parameter ",PARFILE,STATUS)
      if (TASK == 3 .or. TASK == 4) then
        CALL OPENFILE (63,"new","lightcurve",LCFILE, STATUS)
        CALL OPENFILE (64,"new","model fit ",FITFILE,STATUS)
      end if
      if (TASK == 5) CALL OPENFILE(63,"new","refit     ",LCFILE, STATUS)
      if (TASK == 7) CALL OPENFILE(63,"new","bootstrap ",LCFILE, STATUS)
      if (TASK == 8) CALL OPENFILE(63,"new","simulation",LCFILE, STATUS)
      if (TASK == 9) CALL OPENFILE(63,"new","simulation",LCFILE, STATUS)
      if ( STATUS /= 0 ) return


      if ( TASK >= 3 ) then
        write (62,'(A38,A42)') "======================================",
     &                     "=========================================="
        write (62,'(A38,A42)') "JKTEBOP  output results               ",
     &                     "  John Southworth (jkt~astro.keele.ac.uk) "
        write (62,'(A38,A42)') "======================================",
     &                     "=========================================="

        CALL date_and_time(DTDATE,DTTIME)
        write (62,'(A32,I2,6X,A24,A2,":",A2,1X,A2,"/",A2,"/",A4,A8)')
     &                       "Version number of JKTEBOP code: ",VERSION,
     &               "Time and date at start: ",DTTIME(1:2),DTTIME(3:4),
     &                               DTDATE(7:8),DTDATE(5:6),DTDATE(1:4)
        write (62,*)   " "

        write (62,'(A31,A30)') "Input parameter file:          ", INFILE
        write (62,'(A31,A30)') "Input light curve file:        ",OBSFILE
        write (62,'(A31,A30)') "Output parameter file:         ",PARFILE
        if ( TASK == 3 .or. TASK == 4 ) then
          write(62,'(A31,A30)')"Output data file:              ", LCFILE
          write(62,'(A31,A30)')"Output light curve fit:        ",FITFILE
        else if ( TASK == 5 ) then
          write(62,'(A31,A30)')"Perturbed refit output file:   ", LCFILE
        else if ( TASK == 7 ) then
          write(62,'(A31,A30)')"Bootstrapping output file:     ", LCFILE
        else if ( TASK == 8 .or. TASK == 9 ) then
          write(62,'(A31,A30)')"Monte Carlo simulation file:   ", LCFILE
        end if
        write (62,*)   " "
      end if

!-----------------------------------------------------------------------
            ! Read in the photometric data

      if ( TASK >= 3 .and. TASK <= 9 ) then
        CALL READDATA (61,OBSFILE,DATA,DTYPE,NDATA,STATUS)
        if ( STATUS /= 0 ) return
      end if

!-----------------------------------------------------------------------
            ! Read in any additional information from the input file.
            ! TMIN: time of minimum light
            ! LRAT: spectroscopic light ratio
            ! THDL: third light measurement
            ! SINE: sine curve parameters
            ! POLY: polynomial parameters
            ! ECSW: either ecosw or eccentricity (depending on V(7))
            ! ESNW: either esinw or w (omega)    (depending on V(7))

      NLR = 0
      NMIN = 0
      do i = 1,9999
        read (60,'(A200)',iostat=ERROR) CHARHELP200
        if ( ERROR /= 0 ) exit
        CHARHELP4 = '#   '

        read (CHARHELP200,*,iostat=ERROR) CHARHELP4
!-----------------------------------------------------------------------
        if ( CHARHELP4 == "TMIN" .or. CHARHELP4 == "tmin" ) then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,(ARRAY(j),j=1,3)
          if ( ERROR /= 0 ) then
            write(6,'(A32,A28,I3)') "### ERROR reading in data for ti",
     &                             "me of minimum light, number ",NMIN+1
            STATUS = 1
          else
            NMIN = NMIN + 1
            NDATA = NDATA + 1
            DATA(1,NDATA) = ARRAY(1)
            DATA(2,NDATA) = ARRAY(2)
            DATA(3,NDATA) = ARRAY(3)
            DTYPE(NDATA) = 3
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "LRAT" .or. CHARHELP4 == "lrat" ) then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,(ARRAY(j),j=1,3)
          if ( ERROR /= 0 ) then
            write(6,'(A32,A18,I3)') "### ERROR reading in data for li",
     &                                        "ght ratio, number ",NLR+1
            STATUS = 1
          else
            NLR = NLR + 1
            NDATA = NDATA + 1
            DATA(1,NDATA) = ARRAY(1)
            DATA(2,NDATA) = ARRAY(2)
            DATA(3,NDATA) = ARRAY(3)
            DTYPE(NDATA) = 2
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "THDL" .or. CHARHELP4 == "thdl") then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,(ARRAY(j),j=1,2)
          if ( ERROR /= 0 ) then
            write(6,'(A32,A24,I3)') "### ERROR reading in data for th",
     &                               "ird light value, number ",NL3+1
            STATUS = 1
          else
            NL3 = NL3 + 1
            NDATA = NDATA + 1
            DATA(1,NDATA) = 0.0d0         ! dummy value as it is unused
            DATA(2,NDATA) = ARRAY(1)
            DATA(3,NDATA) = ARRAY(2)
            DTYPE(NDATA) = 4
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "SINE" .or. CHARHELP4 == "sine" ) then
          read (CHARHELP200,*,iostat=ERROR)  CHARHELP4, WHICHPAR,
     &                                (ARRAY(j),j=1,3),(IARRAY(j),j=1,3)
          if ( ERROR /= 0 ) then
            write(6,'(A47,I3)')
     &         "### ERROR reading in data for sine wave number ",NSINE+1
            STATUS = 1
          else
            NSINE = NSINE + 1
            V(20+3*NSINE) = ARRAY(1)
            V(21+3*NSINE) = ARRAY(2)
            V(22+3*NSINE) = ARRAY(3)
            VARY(20+3*NSINE) = IARRAY(1)
            VARY(21+3*NSINE) = IARRAY(2)
            VARY(22+3*NSINE) = IARRAY(3)
            if ( WHICHPAR == "J"  ) PSINE(NSINE) = 1
            if ( WHICHPAR == "r1" ) PSINE(NSINE) = 2
            if ( WHICHPAR == "r2" ) PSINE(NSINE) = 3
            if ( WHICHPAR == "i"  ) PSINE(NSINE) = 6
            if ( WHICHPAR == "L3" ) PSINE(NSINE) = 15
            if ( WHICHPAR == "sf" ) PSINE(NSINE) = 17
            if ( WHICHPAR == "L1" ) PSINE(NSINE) = -1
            if ( WHICHPAR == "L2" ) PSINE(NSINE) = -2
            if ( PSINE(NSINE) == 0 ) then
              write(6,'(A30,A32,A4)') "### ERROR: valid sine paramete",
     &                   "r has not been specified. It is ",PSINE(NSINE)
              write(6,'(A32,A32)') "### and needs to be 'J', 'r1', 'r",
     &                              "2', 'i', 'L3', 'sf', 'L1' or 'L2' "
              STATUS = 1
            else
            write(6,'(A28,I1,A30,A2,A1)')">> Read parameters for sine ",
     &               NSINE,', to be applied to parameter "',WHICHPAR,'"'
              write(62,'(A25,I1,A30,A2,A1)')"Read parameters for sine ",
     &               NSINE,', to be applied to parameter "',WHICHPAR,'"'
            end if
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "POLY" .or. CHARHELP4 == "poly" ) then
          read (CHARHELP200,*,iostat=ERROR)  CHARHELP4, WHICHPAR,
     &                                (ARRAY(j),j=1,6),(IARRAY(j),j=1,6)
          if ( ERROR /= 0 ) then
            write(6,'(A48,I3)')
     &        "### ERROR reading in data for polynomial number ",NPOLY+1
            write(6,'(A35,A45)') "Here is an example polynomial input",
     &                  " line for a fit to the total system light:   "
            write(6,'(A35,A45)') "poly  sf  55599.39446   0.0  0.0  0",
     &                  ".0  0.0  0.0   0  1  0  0  0  0              "
            write(6,'(A35,A45)') "which I have used before. Remember ",
     &                  "to specify a good pivot point for your data. "
            STATUS = 1
          else
            NPOLY = NPOLY + 1
            j = 32 + 6*NPOLY
            do k = 1,6
              V(j+k-1) = ARRAY(k)
              VARY(j+k-1) = IARRAY(k)
            end do
            if ( WHICHPAR == "J"  ) PPOLY(NPOLY) = 1
            if ( WHICHPAR == "r1" ) PPOLY(NPOLY) = 2
            if ( WHICHPAR == "r2" ) PPOLY(NPOLY) = 3
            if ( WHICHPAR == "i"  ) PPOLY(NPOLY) = 6
            if ( WHICHPAR == "L3" ) PPOLY(NPOLY) = 15
            if ( WHICHPAR == "sf" ) PPOLY(NPOLY) = 17
            if ( WHICHPAR == "L1" ) PPOLY(NPOLY) = -1
            if ( WHICHPAR == "L2" ) PPOLY(NPOLY) = -2
            if ( PPOLY(NPOLY) == 0 ) then
              write(6,'(A30,A38,A4)') "### ERROR: valid polynomial pa",
     &             "rameter has not been specified. It is ",PPOLY(NPOLY)
              write(6,'(A33,A40)') "### and needs to be one of 'J', '",
     &                        "r1', 'r2', 'i', 'L3', 'sf', 'L1' or 'L2'"
              STATUS = 1
            else
              write(6,'(A22,A12,I1,A30,A2,A1)')">> Read parameters for",
     &" polynomial ",NPOLY,', to be applied to parameter "',WHICHPAR,'"'
              write(62,'(A20,A11,I1,A30,A2,A1)') "Read parameters for ",
     & "polynomial ",NPOLY,', to be applied to parameter "',WHICHPAR,'"'
            end if
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "ECSW" .or. CHARHELP4 == "ecsw") then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,(ARRAY(j),j=1,2)
          if ( ERROR /= 0 ) then
            if ( V(7) < 5.0d0 ) then
              write(6,'(A30,A33,I3)') "### ERROR reading in data for ",
     &                        "e*cos(omega) observation, number ",NECW+1
            else
              write(6,'(A30,A33,I3)') "### ERROR reading in data for ",
     &                        "eccentricity observation, number ",NECW+1
            end if
            STATUS = 1
          else
            NECW = NECW + 1
            NDATA = NDATA + 1
            DATA(1,NDATA) = 0.0d0         ! dummy value as it is unused
            DATA(2,NDATA) = ARRAY(1)
            DATA(3,NDATA) = ARRAY(2)
            DTYPE(NDATA) = 5
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "ESNW" .or. CHARHELP4 == "esnw") then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,(ARRAY(j),j=1,2)
          if ( ERROR /= 0 ) then
            if ( V(7) < 5.0d0 ) then
              write(6,'(A30,A33,I3)') "### ERROR reading in data for ",
     &                        "e*sin(omega) observation, number ",NECW+1
            else
              write(6,'(A30,A41,I3)') "### ERROR reading in data for ",
     &                "periastron longitude observation, number ",NECW+1
            end if
            STATUS = 1
          else
            NESW = NESW + 1
            NDATA = NDATA + 1
            DATA(1,NDATA) = 0.0d0         ! dummy value as it is unused
            DATA(2,NDATA) = ARRAY(1)
            DATA(3,NDATA) = ARRAY(2)
            DTYPE(NDATA) = 6
          end if
!-----------------------------------------------------------------------
        else if ( CHARHELP4 == "NUMI" .or. CHARHELP4 == "numi") then
          read (CHARHELP200,*,iostat=ERROR) CHARHELP4,NUMINT,NINTERVAL
          if ( ERROR /= 0 ) then
            write(6,'(A35,A45)') "### ERROR reading in instructions f",
     &                  "or numerical integration.                    "
            STATUS = 1
          end if
!-----------------------------------------------------------------------
        end if
      end do
!-----------------------------------------------------------------------

      if ( NMIN > 0 ) write(6,'(A7,I3,A43)') ">> Read",NMIN,
     &                     " times of minimum light from the input file"
      if ( NMIN > 0 ) write (62,'(A4,I3,A43)') "Read",NMIN,
     &                     " times of minimum light from the input file"

      if ( NLR > 0 )  write(6,'(A7,I3,A47)') ">> Read",NLR,
     &                 " spectroscopic light ratios from the input file"
      if ( NLR > 0 )  write (62,'(A4,I3,A47)') "Read",NLR,
     &                 " spectroscopic light ratios from the input file"

      if ( NL3 > 0 )  write(6,'(A7,I3,A45)') ">> Read",NL3,
     &                   " third light measurements from the input file"
      if ( NL3 > 0 )  write (62,'(A4,I3,A45)') "Read",NL3,
     &                   " third light measurements from the input file"

      if ( NECW > 0 ) then
        if ( V(7) < 5.0d0 ) write(6,'(A7,I3,A46)') ">> Read",NECW,
     &                  " e*cos(omega) measurements from the input file"
        if ( V(7) < 5.0d0 ) write (62,'(A4,I3,A46)') "Read",NECW,
     &                  " e*cos(omega) measurements from the input file"
        if ( V(7) > 5.0d0 ) write(6,'(A7,I3,A46)') ">> Read",NECW,
     &                  " eccentricity measurements from the input file"
        if ( V(7) > 5.0d0 ) write (62,'(A4,I3,A46)') "Read",NECW,
     &                  " eccentricity measurements from the input file"
      end if

      if ( NESW > 0 ) then
        if ( V(7) < 5.0d0 ) write(6,'(A7,I3,A46)') ">> Read",NESW,
     &                  " e*sin(omega) measurements from the input file"
        if ( V(7) < 5.0d0 ) write (62,'(A4,I3,A46)') "Read",NESW,
     &                  " e*sin(omega) measurements from the input file"
        if ( V(7) > 5.0d0 ) write(6,'(A7,I3,A46)') ">> Read",NESW,
     &                  " periast-long measurements from the input file"
        if ( V(7) > 5.0d0 ) write (62,'(A4,I3,A46)') "Read",NESW,
     &                  " periast-long measurements from the input file"
      end if

      if ( NSINE > 0 ) write(6,'(A7,I1,A39)') ">> Read ",NSINE,
     &                         " sine wave datasets from the input file"
      if ( NSINE > 0 ) write (62,'(A4,I1,A39)') "Read ",NSINE,
     &                         " sine wave datasets from the input file"

      if ( NPOLY > 0 ) write(6,'(A7,I1,A40)') ">> Read ",NPOLY,
     &                        " polynomial datasets from the input file"
      if ( NPOLY > 0 ) write (62,'(A4,I1,A40)') "Read ",NPOLY,
     &                        " polynomial datasets from the input file"

      if ( NUMINT > 1 ) write(6,'(A23,A43)') ">> Read instructions fo",
     &                    "r numerical integration from the input file"
      if ( NUMINT > 1 ) write (62,'(A30,I4,A18,F7.2,A8)')
     &                        "Numerical integration invoked:",NUMINT,
     &                        " samples covering ",NINTERVAL," seconds"

      close (unit=60,status="keep")
      if ( STATUS /= 0 ) return

      write (62,*) " "

!-----------------------------------------------------------------------

      if ( DATA(3,1) < 0.0d0 ) then
        if ( NMIN>0 .or. NLR>0 .or. NL3>0 .or. NECW>0 .or. NESW>0 ) then
          write(6,'(A37,A43)') "### Warning: no observational errors ",
     &                    "have been found but additional information "
          write(6,'(A37,A43)') "### (times of minimum light, spectros",
     &                    "copic light ratios, e or omega values) have"
          write(6,'(A37,A43)') "### been entered. Proper uncertaintie",
     &                    "s are therefore needed for all input data. "
        end if
      end if

      if ( V(13) < 0.0d0 ) write (62,'(A19,A61)') "Mass ratio is below",
     &   " zero: the stellar shapes will be forced to be spherical.    "
      if ( V(7) < 5.0d0 )  write (62,'(A19,A61)') "Eccentricity is bel",
     &   "ow 5:  [e,omega] are taken to be [e*cos(omega),e*sin(omega)]."
      if ( V(7) > 5.0d0 )  write (62,'(A19,A61)') "Eccentricity is mor",
     &   "e than 5:  e and omega will be taken to be (e+10) and omega. "
      if ( V(2) < 0.0d0 )  write (62,'(A19,A61)') "r1+r2 is less than ",
     &   "zero so (r1+r2) and k will be interpreted to mean -r1 and r2"

      if ( V(7) < 5.0d0 ) then
        if (V(7)<-1.0d0.or.V(7)>1.0d0.or.V(8)<-1.0d0.or.V(8)>1.0d0) then
          write(6,'(A36,A44)')  "### ERROR: e*sin(omega) or e*cos(ome",
     &                    "ga) are unphysical: must be between -1 and 1"
          STATUS = 1
          return
        end if
      end if
                    ! Find reflection coefficients if they are not fixed

      if ( V(2) >= 0.0d0 ) then
        R1 = V(2) / (1.0d0 + V(3))
        R2 = V(2) / (1.0d0 + (1.0d0/V(3)))
      else
        R1 = abs(V(2))
        R2 = V(3)
      end if
      MAG = GETMODEL (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,V(20),1,LP,LS,
     &                                                 NUMINT,NINTERVAL)
      if ( VARY(11) == -1 )  V(11) = 0.4d0 * (LS/(1.0d0-V(15))) * R1**2
      if ( VARY(12) == -1 )  V(12) = 0.4d0 * (LP/(1.0d0-V(15))) * R2**2

      if ( TASK >= 3 .and. TASK <= 9 ) then
        do i = 11,12
100       FORMAT (A4,I2,A37,I1,A36)
          if ( VARY(i) == -1 ) write (62,100)
     &                 "Adj(",i,") = -1 so reflection effect for star ",
     &                      i-10," is calculated from system geometry. "
          if ( VARY(i) == 0 )  write (62,100)
     &                 "Adj(",i,") = -1 so reflection effect for star ",
     &                      i-10," is fixed at the input file value.   "
          if ( VARY(i) == 1 )  write (62,100)
     &                 "Adj(",i,") = -1 so reflection effect for star ",
     &                      i-10," is freely adjusted to the best fit. "
          if ( VARY(i) == 2 )  write (62,100)
     &                 "Adj(",i,") = -1 so reflection effect for star ",
     &                      i-10," is set to input value but perturbed."
        end do

        if ( VARY(i) == -1 .and. NLR > 0 ) then
          write(6,'(A37,A43)') "## Warning: if there is a spectroscop",
     &                    "ic light ratio it's best to directly adjust"
          write(6,'(A37,A43)') "##  reflection coefficients rather th",
     &                    "an calculate them from the system geometry."
        end if
        write (62,*) " "

                  ! Now to avoid problems with partial derivatives later

        if ( V(7)>5.0d0 .and. VARY(7)+VARY(8)/=0)  V(7)=max(V(7),10.001)
        if ( V(7) < 6.0d0 .and. V(7) >= 1.0d0 )    V(7) = 0.0d0
        if ( V(11) < 0.001d0 .and. VARY(11) == 1 ) V(11) = 0.001d0
        if ( V(12) < 0.001d0 .and. VARY(12) == 1 ) V(12) = 0.001d0

            ! Find a good starting value for the light scale factor

        if ( VARY(17) /= 0 ) then
          do i = 1,NDATA
            ARRAY(i) = DATA(2,i)
          end do
          V(17) = sellect(ARRAY,NDATA,int(0.2d0*NDATA))
        end if
      end if

            ! If LD of star B is forced to be same as star A, then set
            ! the (unused) LD of star B to be fixed to avoid problems.

      if ( LDTYPE(2) == 0 ) then
        V(5) = V(4)
        VARY(5) = 0
      end if

      END SUBROUTINE INPUT
!=======================================================================
!=======================================================================
      SUBROUTINE READFF (UNIT,NAME1,LOWER1,UPPER1,VALUE1,
     &                        NAME2,LOWER2,UPPER2,VALUE2,STATUS)
            ! This reads in two double precision numbers and checks if
            ! they are within given bounds. STATUS is set to 1 if error
            ! occurs and keeps its previous value if no error occurs.
      implicit none
      integer UNIT                        ! IN: unit number to read from
      character NAME1*10,NAME2*10         ! IN: names of the  parameters
      real*8 LOWER1,UPPER1,LOWER2,UPPER2  ! IN: allowed ranges of values
      real*8 VALUE1,VALUE2                ! OUT:values of the parameters
      integer STATUS                      ! OUT:set to 1 if error occurs
      integer ERROR                       ! LOCAL:  error flag for input

      ERROR = 0
      read (UNIT,*,iostat=ERROR) VALUE1,VALUE2
      if ( ERROR /= 0 ) then
        write(6,'(A33,A10,A5,A10)')
     &           "### ERROR reading the parameters ",NAME1," and ",NAME2
        STATUS = 1
      end if

      if ( VALUE1 < LOWER1 .or. VALUE1 > UPPER1 ) then
        write(6,'(A24,A10,A4,F12.5)')
     &                    "### ERROR: the value of ",NAME1," is ",VALUE1
        write(6,'(A21,F12.5,A4,F12.5)')
     &                      "The allowed range is ",LOWER1," to ",UPPER1
        STATUS = 1
      end if

      if ( VALUE2 < LOWER2 .or. VALUE2 > UPPER2 ) then
        write(6,'(A24,A10,A4,F12.5)')
     &                    "### ERROR: the value of ",NAME2," is ",VALUE2
        write(6,'(A21,F12.5,A4,F12.5)')
     &                      "The allowed range is ",LOWER2," to ",UPPER2
        STATUS = 1
      end if

      END SUBROUTINE READFF
!-----------------------------------------------------------------------
      SUBROUTINE READF (UNIT,NAME,LOWER,UPPER,VALUE,STATUS)
            ! This reads one double precision number and checks if it is
            ! within given bounds. STATUS is set to 1 if an error occurs
      implicit none
      integer UNIT                  ! IN: number of unit to read from
      character NAME*10             ! IN: name of parameter
      real*8 LOWER,UPPER            ! IN: allowed range of values
      real*8 VALUE                  ! OUT: value of the parameter
      integer STATUS                ! IN/OUT: set to 1 if error occurs
      integer ERROR                 ! LOCAL: error flag for input

      ERROR = 0
      read (UNIT,*,iostat=ERROR) VALUE
      if ( ERROR /= 0 ) then
        write(6,'(A32,A10)') "### ERROR reading the parameter ",NAME
        STATUS = 1
      end if
      if ( VALUE < LOWER .or. VALUE > UPPER ) then
        write(6,'(A24,A10,A4,F12.5)')
     &                      "### ERROR: the value of ",NAME," is ",VALUE
        write(6,'(A21,F12.5,A4,F12.5)')
     &                        "The allowed range is ",LOWER," to ",UPPER
        STATUS = 1
      end if

      END SUBROUTINE READF
!-----------------------------------------------------------------------
      SUBROUTINE READCHAR30 (UNIT,NAME,VALUE,STATUS)
            ! This reads in a 30-char value. STATUS is set to 1 if an
            ! error occurs and keeps its previous value if no error.
      implicit none
      integer UNIT                  ! IN: number of unit to read from
      character NAME*10             ! IN: name of character to be read
      character VALUE*30            ! OUT: value of the character
      integer STATUS                ! IN/OUT: set to 1 if error occurs
      integer ERROR                 ! LOCAL: error flag for input

      ERROR = 0
      read (UNIT,*,iostat=ERROR) VALUE
      if ( ERROR /= 0 ) then
        write(6,'(A29,A30)') "### ERROR reading the string ",NAME
        STATUS = 1
      end if
      if ( VALUE == "#" ) then
        write(6,'(A39,A41)') "### ERROR: you cannot use '#' as an out",
     &                      "put file name.                           "
        STATUS = 1
      end if

      END SUBROUTINE READCHAR30
!-----------------------------------------------------------------------
      SUBROUTINE READ2 (UNIT,NAME1,NAME2,VALUE1,VALUE2,STATUS)
            ! This reads two integers on one line and tests if both are
            ! between -1 and 3. If not then STATUS is set to 1.
            ! If there is no error then STATUS keeps its previous value.
      implicit none
      integer UNIT                  ! IN: number of unit to read from
      character NAME1*10,NAME2*10   ! IN: names of characters to read
      integer VALUE1,VALUE2         ! OUT: values of the integers
      integer STATUS                ! IN/OUT: set to 1 if error occurs
      integer ERROR                 ! LOCAL: error flag for input

      ERROR = 0
      read (UNIT,*,iostat=ERROR) VALUE1,VALUE2
      if ( ERROR /= 0 ) then
        write(6,'(A31,A10,A5,A10)')
     &            "### ERROR reading the integers ",NAME1," and ",NAME2
        STATUS = 1
      end if
      if ( VALUE1 < -1 .or. VALUE1 > 3 ) then
        write(6,'(A30,A10,A4,I6,A28)') "### ERROR: adjustment integer ",
     &                NAME1," is ",VALUE1," and should be -1, 0, 1 or 2"
        STATUS = 1
      end if
      if ( VALUE2 < -1 .or. VALUE2 > 3 ) then
        write(6,'(A30,A10,A4,I6,A28)') "### ERROR: adjustment integer ",
     &                NAME2," is ",VALUE2," and should be -1, 0, 1 or 2"
        STATUS = 1
      end if

      END SUBROUTINE READ2
!=======================================================================
      SUBROUTINE READDATA (UNIT,OBSFILE,DATA,DTYPE,NDATA,STATUS)
      implicit none
      integer UNIT                  ! IN: Unit number  to read file from
      character*30 OBSFILE          ! IN: Name of input light curve file
      real*8 DATA(3,999999)         ! OUT: Obs'l data (time, mag, error)
      integer DTYPE(999999)         ! OUT: Type of data (here all "1")
      integer NDATA                 ! OUT: Number of datapoints read
      integer STATUS                ! IN/OUT: set to 1 if there is error
      integer i,ERROR,ERRFLAG       ! LOCAL: Loop counter + error flags
      character*200 CHARHELP        ! LOCAL: Helper character string
      real*8 HELP1,HELP2,HELP3      ! LOCAL: Helper variables

      ERROR = 0
      CALL OPENFILE (UNIT,"old","data file ",OBSFILE,ERROR)
      if ( ERROR /= 0 ) then
        STATUS = 1
        return
      end if

      read (UNIT,'(A200)',iostat=ERROR) CHARHELP
      rewind (UNIT)

      read (CHARHELP,*,iostat=ERROR) HELP1,HELP2,HELP3
      if ( ERROR /= 0 ) then
        read (CHARHELP,*,iostat=ERROR) HELP1,HELP2
        if ( ERROR /= 0 ) then
          write(6,'(A48,A30)')
     &        "### ERROR: cannot understand first line of file ",OBSFILE
          STATUS = 1
          return
        end if
        ERRFLAG = 0                  ! There are no observational errors
      else
        ERRFLAG = 1                     ! There are observational errors
      end if

      do i = 1,999999
        if ( ERRFLAG == 0 ) then
          read (UNIT,*,iostat=ERROR) DATA(1,i),DATA(2,i)
          DATA(3,i) = -1.0d0
        else if ( ERRFLAG == 1 ) then
          read (UNIT,*,iostat=ERROR) DATA(1,i),DATA(2,i),DATA(3,i)
        end if
        if ( ERROR /= 0 ) exit
        DTYPE(i) = 1
      end do

      NDATA = i - 1
      if ( NDATA < 5 ) then
        write(6,'(A30)') "### ERROR: too few data to fit"
        STATUS = 1
      end if

      if ( ERRFLAG == 1 ) then
        write (62,'(A5,I6,A39,A30)') "Read ",NDATA,
     &                 " datapoints (with errorbars) from file ",OBSFILE
        write(6,'(A8,I6,A36,A30)') ">> Read ",NDATA,
     &                    " datapoints (with errors) from file ",OBSFILE
      else if  ( ERRFLAG == 0 ) then
        write (62,'(A5,I6,A38,A30)') "Read ",NDATA,
     &                  " datapoints (no error bars) from file ",OBSFILE
        write(6,'(A8,I6,A37,A30)') ">> Read ",NDATA,
     &                   " datapoints (no errorbars) from file ",OBSFILE
      end if


      close (UNIT)

      END SUBROUTINE READDATA
!=======================================================================
      SUBROUTINE OPENFILE (UNIT,STATE,FILE,FILENAME,STATUS)
            ! Opens a file. STATUS = 1 if the action was not
            ! successful and left unchanged if the action was successful
      implicit none
      integer UNIT                  ! IN: unit number to open file to
      character*30 FILENAME         ! IN: name of datafile to open
      character*10 FILE             ! IN: identifier of this datafile
      character*3 STATE             ! IN: "old" or "new"
      integer STATUS                ! OUT: indicates success of opening
      integer ERROR                 ! LOCAL: error flag for opening file

      ERROR = 0
      OPEN (unit=UNIT,file=FILENAME,status=STATE,iostat=ERROR)
      if ( ERROR /= 0 ) then
        write(6,'(A18,A3,1X,A10,A8,A30)')
     &              "### ERROR opening ", STATE,FILE," file:  ",FILENAME
        STATUS = 1
      end if
      if (STATE=="new".and.ERROR==0)  write(6,'(A10,A3,1X,A10,A8,A30)')
     &                       ">> Opened ",STATE,FILE," file:  ",FILENAME

      END SUBROUTINE OPENFILE
!=======================================================================
!=======================================================================
      SUBROUTINE TASK1 ()           ! This task outputs  limb  darkening
            ! coefficients for given Teff and logg and [M/H] and Vmicro.
            ! It simply interfaces with the  JKTLD  code, which performs
            ! bilinear interpolation in Teff,logg for given [M/H],Vmicro
            ! Usage:   jktld  <Teff>  <logg>  <M/H>  <Vmicro>  <outfile>
      implicit none
      real*8 TEFF,LOGG              ! IN: Parameters  to  interpolate to
      real*8 MOH,VMICRO             ! IN: Other parameters forthe tables
      character*20 CTEFF,CLOGG      ! LOCAL: character version of values
      character*20 CMOH,CMICRO      ! LOCAL: character version of values
      character*30 OUTFILE          ! LOCAL: name of output file to make

      write(6,'(A40,$)') "Enter the effective temperature (K)  >> "
      read (*,*) TEFF
      write(6,'(A40,$)') "Enter the surface gravity (log cm/s) >> "
      read (*,*) LOGG
      write(6,'(A40,$)') "Enter the metal abundance  ([M/H])   >> "
      read (*,*) MOH
      write(6,'(A40,$)') "Enter the microturbulence velocity   >> "
      read (*,*) VMICRO
      write(6,'(A40,$)') "Enter the output file name to create >> "
      read (*,*) OUTFILE

      if ( TEFF < 3500.0d0 ) then
        write(6,'(A39,A41)') "### Warning: a Teff below 3500 K is out",
     &                      " of range of most of the LD coeff tables."
      else if ( TEFF < 2000.0d0 ) then
        write(6,'(A39,A41)') "### Warning: a Teff below 2000 K is out",
     &                      " of range of the LD coefficient tables.  "
      else if ( TEFF > 50000.0d0 ) then
        write(6,'(A39,A41)') "### Warning: a Teff above 50000 K is ou",
     &                      "t of range of the LD coefficient tables. "
      end if
      if ( LOGG > 5.0d0 .or. LOGG < 0.0d0 )
     &  write(6,'(A39,A41)') "### Warning: log(g)a outside the range ",
     &                      "0.0 to 5.0 are not covered in the tables."
      if ( MOH /= 0.0d0 .and. VMICRO /= 2.0d0 )
     &  write(6,'(A39,A41)') "### Warning: for [M/H] /= 0 and Vmicro ",
     &                      "/= 2 you will probably get no results.   "

      write (CTEFF,'(F20.8)') TEFF
      write (CLOGG,'(F20.8)') LOGG
      write (CMOH,'(F20.8)') MOH
      write (CMICRO,'(F20.8)') VMICRO

      CALL SYSTEM ( "jktld " // CTEFF // " " //  CLOGG // " " // CMOH
     &                              // " " // CMICRO // " " // OUTFILE )

      END SUBROUTINE TASK1
!=======================================================================
      SUBROUTINE TASK2 (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY)
      implicit none                       ! Produces a model light curve
      real*8 V(67)                        ! IN: light  curve  parameters
      integer LDTYPE(2)                   ! IN: LD law type foreach star
      integer i,ERROR                     ! LOCAL: counters & error flag
      real*8 MAG,LP,LS                    ! LOCAL: EBOP/GETMODEL  output
      real*8 PHASE                        ! LOCAL:  phase for evaluation
      real*8 GETMODEL                     ! FUNCTION: evaluate the model
      integer NSINE                       ! OUT: Numbrs of sines and L3s
      integer PSINE(5)                    ! OUT: Which par for each sine
      integer NPOLY,PPOLY(5)              ! OUT: Similar for polynomials

      real*8 HJD
      real*8 R1,R2
      integer NPHASE                      ! LOCAL: number of phases todo

      V(19) = 1.0d0           ! Set period to 1.0
      V(20) = 0.0d0           ! Set Tzero to 0.0
      LP = 0.0d0
      LS = 0.0d0
                                      ! NSINE=0 and NPOLY=0 and NUMINT=1
      if ( V(2) >= 0.0d0 ) then
        R1 = V(2) / (1.0d0 + V(3))
        R2 = V(2) / (1.0d0 + (1.0d0/V(3)))
      else
        R1 = abs(V(2))
        R2 = V(3)
      end if
      MAG=GETMODEL(V,LDTYPE,0,PSINE,0,PPOLY,V(20),1,LP,LS,1,0.0d0)
      V(11) = 0.4d0 * (LS/(1.0d0-V(15))) * R1**2
      V(12) = 0.4d0 * (LP/(1.0d0-V(15))) * R2**2

      NPHASE = 10001
      if ( R1 < 0.01d0 .or. R2 < 0.01d0 ) NPHASE = 100001

      write(6,'(A40,A40)') ">> The reflection coefficients come from",
     &                      " the system geometry, not the input file"

      write(62,'(A47)')"#  PHASE  MAGNITUDE    L1         L2         L3"

      do i = 1,NPHASE
        PHASE = (i-1) / dble(NPHASE-1)
        HJD = V(20) + PHASE * V(19)
        MAG = GETMODEL(V,LDTYPE,0,PSINE,0,PPOLY,HJD,1,LP,LS,1,0.0d0)
        write (62,'(F8.6,4(1X,F10.6))') PHASE,MAG,LP,LS,V(15)
      end do
      close (62)

      END SUBROUTINE TASK2
!=======================================================================
      SUBROUTINE TASK34 (TASK,V,VARY,LDTYPE,DATA,DTYPE,NDATA,NLR,NMIN,
     &                   SIGMA,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,
     &                   NUMINT,NINTERVAL)
            ! Find the best-fitting light curve parameters for the data.
            ! If TASK=4 it iteratively rejects all datapoints which are
            ! > SIGMA sigma from the best fit and refits the remainder.
      implicit none
      integer TASK                        ! IN: which task to do
      real*8 V(67)                        ! IN: light  curve  parameters
      integer VARY(67)                    ! IN: parameters vary or fixed
      integer LDTYPE(2)                   ! IN: LD law type foreach star
      real*8 DATA(3,999999)               ! IN: time, magnitude, error
      integer DTYPE(999999)               ! IN: type of each datapoint
      integer NDATA,NLR,NMIN              ! IN: number of  types of data
      integer NSINE,NL3,NECW,NESW         ! IN: Numbers of sines and L3s
      integer PSINE(5)                    ! IN: Which par for each sine
      integer NPOLY,PPOLY(5)              ! IN: Similar for polynomials
      integer NUMINT                      ! IN: Number of numerical ints
      real*8  NINTERVAL                   ! IN: Time interval for numint
      real*8 SIGMA                        ! IN:  std. dev. for  clipping
      real*8 CHISQ                        ! SUB: chi-square of model fit
      real*8 VERR(67)                     ! SUB: parameter formal errors
      real*8 OMC(999999)                  ! LOCAL: (O-C) residual values
      real*8 SIG                          ! LOCAL: rms of the O-C values
      real*8 RESIDSQSUM                   ! LOCAL: sum of resid. squares
      integer ITER,IFAIL                  ! LOCAL: iter number & success
      real*8 MAG,LP,LS                    ! LOCAL: EBOP/GETMODEL  output
      integer KEEP(999999),ACOUNT         ! LOCAL: Datapoint bookkeeping
      integer i,j                         ! LOCAL: loop counter
      real*8 GETMODEL                     ! FUNCTION: evaluate the model
      integer NREJ

      do i = 1,999999
        OMC(i) = 0.0d0
        KEEP(i) = 0
      end do
      SIG = 0.0d0
      NREJ = 0

            ! Output initial values and then find the best fit.
            ! If TASK=3 then output various results and finish.

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &  VERR,0,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &VERR,IFAIL,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      if ( TASK == 3 ) then
        CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &    VERR,1,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

        if ( IFAIL == 0 .and. ITER < 200 ) then
          write(6,'(A33,$)') ">> Best fit has been found after "
          if ( ITER < 100 ) write(6,'(I2,$)') ITER
          if ( ITER >= 100 ) write(6,'(I3,$)') ITER
          write(6,'(A12)') " iterations."
        else
          write(6,'(A38)') "### WARNING: a good fit was not found."
        end if

!         write(6,'(A46)')">> Final parameters written to parameter file."
!         write(6,'(A45)')">> Data and residuals written to light curve file."
!         write(6,'(A45)')">> Best-fitting model is written to fit file."

        return
      end if

            ! If TASK=4 then output the best fit and continue.

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &    VERR,2,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
      write(6,'(A40,A40)')  ">> Best fit to all data points has been ",
     &                       "found and written to the parameter file."

            ! Now, for nine iterations, calculate the scatter of the
            ! observations (if there are no observational errors). Then
            ! go through the data checking for ones with O-C values
            ! greater than SIGMA sigmas. Finish if an iteration does not
            ! result in a rejected datapoint. Remove the rejected datap-
            ! oints then refit the data. Note that the light ratios and
            ! minimum times always have uncertainties.

      do j = 1,9

            ! First calculate the rms of the light curve residuals if
            ! there are no observational uncertainties.

        if ( DATA(3,1) < 0.0d0 ) then
          RESIDSQSUM = 0.0d0
          do i = 1,NDATA
            if ( DTYPE(i) == 1 ) then
              MAG = GETMODEL (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                        DATA(1,i),DTYPE(i),LP,LS,NUMINT,NINTERVAL)
              OMC(i) = MAG - DATA(2,i)
              RESIDSQSUM = RESIDSQSUM + OMC(i)**2
            end if
          end do
          SIG = sqrt( RESIDSQSUM / dble(NDATA-NLR-NMIN) )
        else
          do i = 1,NDATA
            if ( DTYPE(i) == 1 ) then
              MAG = GETMODEL (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                        DATA(1,i),DTYPE(i),LP,LS,NUMINT,NINTERVAL)
              OMC(i) = MAG - DATA(2,i)
            end if
          end do
        end if

          ! Now put the array indices of the good datapoints into KEEP

        ACOUNT = 0
        do i = 1,NDATA                      ! if no observational errors
          if ( DTYPE(i) == 1 .and. DATA(3,i) < 0.0d0 ) then
            if ( abs(OMC(i)) <= SIG*SIGMA ) then
              ACOUNT = ACOUNT + 1
              KEEP(ACOUNT) = i
            end if
          else                            ! if ob'l errors were supplied
            if ( abs(OMC(i)/DATA(3,i)) <= SIGMA ) then
              ACOUNT = ACOUNT + 1
              KEEP(ACOUNT) = i
            end if
          end if
        end do

            ! Now keep only those datapoints which are specified by KEEP
            ! Have to recount number of light ratios and minimum times

        do i = 1,ACOUNT
          DATA(1,i) = DATA(1,KEEP(i))
          DATA(2,i) = DATA(2,KEEP(i))
          DATA(3,i) = DATA(3,KEEP(i))
          DTYPE(i) = DTYPE(KEEP(i))
        end do

        NREJ = NDATA - ACOUNT
        write(6,'(A13,I1,A1,I5,A8,I5,A47)') ">> Iteration ",j,":",NREJ,
     &                                       " out of ",NDATA,
     &                 " datapoints have been rejected from the dataset"

            ! Now recount the numbers and types of data and refit them.

        NDATA = ACOUNT
        NLR = 0
        NMIN = 0
        do i = 1,NDATA
          if ( DTYPE(i) == 2 ) NLR = NLR + 1
          if ( DTYPE(i) == 3 ) NMIN = NMIN + 1
        end do

        CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,
     &                CHISQ,VERR,IFAIL,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,
     &                NESW,NUMINT,NINTERVAL)
        if ( NREJ < 1 ) exit
      end do

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &    VERR,1,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
      write(6,'(A32,I3,A45)') ">> Best fit has been found from ",ITER,
     &                 " iterations and written to the parameter file"

      END SUBROUTINE TASK34
!=======================================================================
      SUBROUTINE TASK6 (V,VARY,LDTYPE,DATA,DTYPE,NDATA,NLR,NMIN,NSINE,
     &                 PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
            ! Fits model parameters to an observed light curve. Then
            ! investigates each adjustable parameter by fixing it at a
            ! range of values round the best fit and seeing what the
            ! effect is on the other parameters.
      implicit none
      real*8 V(67)                        ! IN: The light curve params
      integer VARY(67)                    ! IN: Par adjustment integers
      integer LDTYPE(2)                   ! IN: LD law type foreach star
      real*8 DATA(3,999999)               ! IN: Observational data
      integer DTYPE(999999)               ! IN: Type of each datapoint
      integer NDATA,NLR,NMIN              ! IN: Numbers of datapoints
      integer NSINE,NL3,NECW,NESW         ! IN: Numbers of sines and L3s
      integer PSINE(5)                    ! IN: Which par for each sine
      integer NPOLY,PPOLY(5)              ! IN: Similar for polynomials
      integer NUMINT                      ! IN: Number of numerical ints
      real*8  NINTERVAL                   ! IN: Time interval for numint
      integer NUMVARY                     ! LOCAL: Number of vary params
      integer VWHERE(67)                  ! LOCAL: Which parameters vary
      real*8 VSTORE(67)                   ! LOCAL: Store best-fit params
      integer VARYFLAG                    ! LOCAL: Store a VARY integer
      real*8 VERR(67),CHISQ               ! LOCAL: Output from FITEBOP
      integer i,j,k,IFAIL,ITER            ! LOCAL: Loop counters etc

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &    VERR,0,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
      CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &              VERR,IFAIL,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,
     &              NUMINT,NINTERVAL)
      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &    VERR,2,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
      write(6,'(A40,A40)')  ">> Fitting process completed. Best fit f",
     &                       "ound and output to the parameter file.  "

      NUMVARY = 0
      j = 1
      do i = 1,67
        if ( i /= 11 .and. i /= 12 ) then
          if ( VARY(i) /= 0 ) then
            NUMVARY = NUMVARY + 1
            VWHERE(j) = i
            j = j + 1
          end if
        end if
      end do

      do i = 1,67
        VSTORE(i) = V(i)
      end do

            ! Params with VARY=1 are adjustable and those with VARY=2
            ! are fixed when finding the best fit but are perturbed like
            ! adjustable parameters here.
            ! Foreach perturbed parameter store its value and adjustment
            ! integer, then step through various values whilst adjusting
            ! all adjustable parameters to find the best fit.   Start at
            ! the best-fit value  and gradually step away from it in one
            ! direction. Then do the same for the other direction.

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &    VERR,0,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
      CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &              VERR,IFAIL,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,
     &              NUMINT,NINTERVAL)

      do i = 1,NUMVARY
        j = VWHERE(i)
        do k = 1,67
          V(k) = VSTORE(k)
        end do
        j = VWHERE(i)
        if ( j /= 16 .and. j /= 17 ) then
          VARYFLAG = VARY(j)
          VARY(j) = 0
          CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,
     &                 CHISQ,VERR,100+j,NSINE,PSINE,NPOLY,PPOLY,NL3,
     &                 NECW,NESW,NUMINT,NINTERVAL)

          do k = 0,67
            if ( j == 6 ) V(j) = VSTORE(j) - 0.1d0*k
            if ( j /= 6 ) V(j) = VSTORE(j) * (1.0d0 - k/40.0d0)
            CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,
     &                    CHISQ,VERR,IFAIL,NSINE,PSINE,NPOLY,PPOLY,NL3,
     &                    NECW,NESW,NUMINT,NINTERVAL)
            CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,
     &                   CHISQ,VERR,200+j,NSINE,PSINE,NPOLY,PPOLY,NL3,
     &                   NECW,NESW,NUMINT,NINTERVAL)
          end do

          do k = 1,67
            V(k) = VSTORE(k)
          end do

          do k = 1,67
            if ( j == 6 ) V(j) = VSTORE(j) + 0.1d0*k
            if ( j /= 6 ) V(j) = VSTORE(j) * (1.0d0 + k/40.0d0)
            if ( j /= 6 .or. V(6) < 90.1d0 ) then
              CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,
     &                      ITER,CHISQ,VERR,IFAIL,NSINE,PSINE,NPOLY,
     &                      PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
              CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,
     &                     ITER,CHISQ,VERR,200+j,NSINE,PSINE,NPOLY,
     &                     PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
            end if
          end do

          VARY(j) = VARYFLAG
        end if

      end do

      END SUBROUTINE TASK6
!=======================================================================
      SUBROUTINE TASK5789 (TASK,V,VARY,LDTYPE,DATA,DTYPE,NDATA,NLR,NMIN,
     &      NSIM,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
!
!   This big subroutine executes four of the JKTEBOP tasks:
!
! TASK 5:  this has some characteristics of a global search. Once a best
!   fit has been found the params are perturbed and refitted many times.
!   Each result is outputted  and the best one is given as the solution.
!
! TASK 7:  this performs a bootstrapping error analysis. Once a best fit
!   is found many new datasets are made by copying (with replacement) at
!   random from the actual light curve. Each one is fitted to find 1-sig
!   errors. All light ratios and times of minimum are kept each time.
!
! TASK 8:  this does a Monte Carlo simulation error analysis.  Once best
!   fit has been found, the model is evaluated  at the phases of the obs
!   to make a synthetic light curve. Then, for each simulation, Gaussian
!   noise is added and the data are refitted to find the 1-sigma spread.
!
! TASK 9:  this does a Monte Carlo simulation error analysis similar to
!   TASK 8, but instead of adding simulated Gaussian noise to the data,
!   the residuals of the fit are used. For each simulation the set of
!   residuals is moved and applied to the model values for the adjacent
!   datapoints. Residuals which move off the end of the data are wrapped
!   around to the start of the dataset.  Thus the number of simulations
!   done is (NDATA - 1). This has the advantage that correlated noise in
!   the residuals is retained and affects the quality of the fit.
!
! For Tasks 7, 8, 9, the start params are perturbed each time to avoid
!   finding unrealistically small errors in a local minimum. For Tasks
!   7 and 8 this is optional (won't be done if NSIM is less than zero).
!   Can also be turned off for individual params by putting its VARY=3
!
! Monte Carlo simulations are an excellent error indicator, but do not
!   properly account for systematic errors.  The residual permutation
!   (Task 9) is better. Bootstrapping is the best, but is likely to give
!   pesimistic results due to the loss of time sampling in the simulated
!   datasets. Correct solution to problem: model and remove systematics.
!
! For information on bootstrapping, Monte Carlo, and minimisation algor-
! ithms, read Numerical Recipes in Fortran 77 (Press et al 1993) chap.15

      implicit none
      integer TASK                  ! IN: The task to undertake
      real*8 V(67)                  ! IN: The photometric parameters
      integer VARY(67)              ! IN: Parameter adjustment integers
      integer LDTYPE(2)             ! IN: Type of LD law for each star
      real*8 DATA(3,999999)         ! IN: Observational data
      integer DTYPE(999999)         ! IN: Types of observational data
      integer NDATA,NLR,NMIN        ! IN: Numbers of different datatypes
      integer NSINE,NL3,NECW,NESW   ! IN: Numbers of sines and L3 values
      integer PSINE(5)              ! IN: Which parameter for each sine
      integer NPOLY,PPOLY(5)        ! IN: Similar for the polynomials
      integer NSIM                  ! IN: Number of simulations to do
      integer NUMINT                ! IN: Number of numerical integrat's
      real*8  NINTERVAL             ! IN: Time interval numerical integ.
      character PERTURB*1           ! LOCAL: Perturb params ('y' or 'n')
      real*8 VERR(67)               ! SUB:   Param  formal uncertainties
      real*8 DV(67)                 ! LOCAL: Perturbations for paramters
      real*8 VEXTRA(6)              ! SUB:   Extra (dependent) parametrs
      real*8 VALL(67,100000)        ! SUB:   All the best-fitting params
      real*8 VALLEXTRA(6,100000)    ! SUB:   All extra (dependent) pars
      integer NVARY                 ! LOCAL: Number of varying parametrs
      integer ITER                  ! SUB:   Numberof fitting iterations
      real*8 INDATA(3,999999)       ! LOCAL: Synthetic obs  to be fitted
      real*8 RESID(1999998)         ! LOCAL: Residuals of the best fit
      real*8 INV(67)                ! SUB:   Parameters of  synthetic LC
      real*8 CHISQ,ERRSIZE          ! SUB:   Chi-squared of the best fit
      integer NPHOT                 ! LOCAL: Number of photometric datap
      integer SEEDSTART,SEED        ! FUNCTION: Start randomnumber maker
      real*8 RANDOMG                ! FUNCTION: Gaussian  random  number
      real*8 RANDOM                 ! FUNCTION: Flat-distrib  random num
      real*8 GETMODEL               ! FUNCTION: Gets  model  predictions
      real*8 LP,LS,MAG              ! LOCAL: EBOP/GETMODEL output values
      real*8 HELP1                  ! LOCAL: Useful variable storage
      integer i,j,k,m,ERROR         ! LOCAL: Loop counters + error flags
      real*8 ECQPHASES(4)           ! phases of eclipse and quadrature

      NVARY = 0
      do i = 1,67
        if ( VARY(i) == 1 ) NVARY = NVARY + 1
      end do

            ! First get DV values (how much to perturb parameters by).
            ! Set 4x larger than the numerical derivative intervals.
            ! Check whether perturbations should be applied.
            ! Set perturbations a factor 4 of larger for TASK 5.

      CALL GET_DV (V,DV)
      do i = 1,67
        DV(i) = DV(i) * 4.0d0
      end do

      if ( TASK == 9 ) NSIM = NDATA - 1
      if ( NSIM < 0 ) then
        PERTURB = "n"
        NSIM = abs(NSIM)
      else
        PERTURB = "y"
      end if

      if ( TASK == 5 ) then
        do i = 1,67
          DV(i) = DV(i) * 4.0d0
        end do
      end if

            ! Now find the best fit to the light curve and output it.

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &    VERR,0,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &              VERR,ERROR,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,
     &              NUMINT,NINTERVAL)

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &    VERR,2,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

      write(6,'(A40,A40)')  ">> Fitting process completed. Best fit i",
     &                       "s found and outputted to parameter file. "

            ! Store some dependent quantities for output later on.
            ! Write column headings for simulation results output file.

      if ( V(2) >= 0.0d0 ) then
        VEXTRA(1) = V(2) / (1.0d0 + V(3))                   ! r_1
        VEXTRA(2) = V(2) / (1.0d0 + (1.0d0/V(3)))           ! r_2
      else
        VEXTRA(1) = abs(V(2)) + V(3)                        ! r1+r2
        VEXTRA(2) = V(3) / abs(V(2))                        ! k
      end if

      call ECQUADPHASES(V(7),V(8),V(16),ECQPHASES)
      VEXTRA(3) = GETMODEL (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     & V(20)+V(19)*ECQPHASES(2),2,LP,LS,NUMINT,NINTERVAL)   ! L_2 / L_1

      if ( V(7) >= 10.0d0 ) then
        VEXTRA(4) = (V(7)-10.0d0) * cos(V(8)/57.2957795d0)
        VEXTRA(5) = (V(7)-10.0d0) * sin(V(8)/57.2957795d0)
      else
        VEXTRA(4) = sqrt(V(7)**2 + V(8)**2)                 ! e
        VEXTRA(5) = atan2(V(8),V(7))*57.2957795d0           ! omega
        if ( VEXTRA(5) < 0.0d0 )     VEXTRA(5) = VEXTRA(5) + 360.0d0
        if ( VEXTRA(5) >= 360.0d0 )  VEXTRA(5) = VEXTRA(5) - 360.0d0
      end if
      VEXTRA(6) = CHISQ / (NDATA - NVARY)
      if ( DATA(3,1) <= 0.0d0 ) VEXTRA(6) = -999.0d0

      if ( V(2) >= 0.0d0 ) then
      write(63,'(A,$)') "#   N  ITER     SB2                r1+r2    "//
     &     "           k                 LDu1               LDu2     "//
     &     "           i                 ecosw              esinw    "//
     &     "          GD1                GD2                refl1    "//
     &     "          refl2               q                 tidalangl"//
     &     "e         L3                 phasecorr          sfact    "//
     &     "          intring             P                Tzero     "//
     &     "          LDn1               LDn2               r1       "//
     &     "          r2                 L2/L1               e       "//
     &     "         omega               chisq              "
      else
      write(63,'(A,$)') "#   N  ITER     SB2                r1       "//
     &     "           r2                LDu1               LDu2     "//
     &     "           i                 ecosw              esinw    "//
     &     "          GD1                GD2                refl1    "//
     &     "          refl2               q                 tidalangl"//
     &     "e         L3                 phasecorr          sfact    "//
     &     "          intring             P               Tzero      "//
     &     "          LDn1               LDn2               r1+r2    "//
     &     "          k                  L2/L1               e       "//
     &     "         omega               chisq              "
      end if

      do i = 1,5
        if ( NSINE > i-1 )  write (63,'(A5,I1,A18,I1,A18,I1,A13,$)')
     &                            "sine_",i,"_T0          sine_",i,
     &                            "_P           sine_",i,"_amp         "
      end do

      do i = 1,5
        if ( NPOLY > i-1 )  write (63,'(2(A5,I1,A18,I1,A18,I1,A13),$)')
     &          "poly_",i,"_pivot       poly_",i,"_x_coeff     poly_",i,
     &          "_x^2_coeff   ","poly_",i,"_x^3_coeff   poly_",i,
     &          "_x^4_coeff   poly_",i,"_x^5_coeff   "
      end do

      write (63,*) " "


            ! Start the random number generator. Calculate the number of
            ! variable params.  Store original data and best-fit params.

      SEED = SEEDSTART ()

      do i = 1,67
        INV(i) = V(i)
      end do

      do i = 1,NDATA
        INDATA(1,i) = DATA(1,i)
        INDATA(2,i) = DATA(2,i)
        INDATA(3,i) = DATA(3,i)
      end do

      NPHOT = 0
      do i = 1,NDATA
        if ( DTYPE(i) == 1 ) NPHOT = NPHOT + 1
      end do

            ! If the inputted NSIM is below zero no perturbations should
            ! be applied to the initial  light curve parameter estimates
            ! prior to each MC simulation. This is useful if convergence
            ! only happens for a  very narrow range of parameter values.

      if ( TASK == 5 ) write (62,'(A23,A57)') "Task 5: initial paramet",
     &      "ers will be widely perturbed before each fit is done.    "
      if ( TASK == 7 .or. TASK == 8) then
        if (PERTURB == "n") write (62,'(A18,A60)') "Number of simulati",
     &   "ons is below zero: the parameters will not be perturbed.    "
        if (PERTURB == "y") write (62,'(A18,A60)') "Number of simulati",
     &   "ons is above zero: the parameters will be perturbed.        "
        write (62,*) " "
      end if
      if ( TASK == 9 ) write (62,'(A23,A57)') "Task 9: initial paramet",
     &      "ers will be perturbed before each fit is performed.      "

      if ( TASK == 5 ) write(6,'(A46,I6)')
     &             ">> Number of perturbed parameter refits to do:",NSIM
      if ( TASK == 7 ) write(6,'(A45,I6)')
     &              ">> Number of bootstrapping simulations to do:",NSIM
      if ( TASK == 8 ) write(6,'(A43,I6)')
     &                ">> Number of Monte Carlo simulations to do:",NSIM
      if ( TASK == 9 ) write(6,'(A43,I6)')
     &                ">> Number of Monte Carlo simulations to do:",NSIM
      write(6,'(A13,$)') ">> Completed:"

            ! For Monte Carlo simulations we must evaluate the best-fit-
            ! ting model light curve at the observed phases.
            ! TASK 8: if there were no errors in the input light curve,
            ! calculate the rms of the residuals of the fit and use to
            ! scale the Gaussian noise added to the simulated data.
            ! TASK 9: store the residuals. Do this twice over as this
            ! makes it easy later on to slide the residuals along the
            ! data. This is equivalent to allowing the residuals which
            ! drop off the end to be wrapped to the start of the data.

      if ( TASK == 8 .or. TASK == 9 ) then

        HELP1 = 0.0d0
        do i = 1,NDATA
          if ( DTYPE(i) == 1 ) then
            DATA(2,i) = GETMODEL (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                        DATA(1,i),DTYPE(i),LP,LS,NUMINT,NINTERVAL)
            HELP1 = HELP1 + (DATA(2,i) - INDATA(2,i))**2
          end if
        end do

        if ( TASK == 9 ) then
          do i = 1,NDATA
            RESID(i) = INDATA(2,i) - DATA(2,i)
            RESID(i+NPHOT) = RESID(i)
          end do
        end if

        if ( TASK == 8 ) then
          if ( DATA(3,1) <= 0.0d0 ) then
            ERRSIZE = sqrt(HELP1 / (NDATA-NVARY-NLR-NMIN) )
            do i = 1,NDATA
              if ( DTYPE(i) == 1 ) INDATA(3,i) = ERRSIZE
            end do
            write (62,'(A34,A46)') "No observational errors were suppl",
     &                 "ied with the input light curve.  Noise of size"
            write (62,'(F6.3,A12,A62)')abs(ERRSIZE)*1.d3," mmag (stand",
     & "ard error of best fit) will be added to the synthetic datasets"
          else
            write (62,'(A34,A46)') "Observational errors were supplied",
     &                 " with the input light curve.  These have been "
            write (62,'(A34,A46)') "assumed to be correct and used to ",
     &                 "set the size of the simulated Gaussian noise. "
          end if
        end if
      end if

!-----------------------------------------------------------------------

      do i = 1,abs(NSIM)             ! Loop over number  of  simulations
500     continue                     ! Enter here if previous sim failed

            ! TASK 5: use actual dataset for every simulation
        if ( TASK == 5 ) then
          do j = 1,NDATA
            INDATA(2,j) = DATA(2,j)
            INDATA(3,j) = DATA(3,j)
          end do
        end if

            ! TASK 7: randomly sample the observations with replacement
            ! to create a new light curve to fit. Don't do this with the
            ! light ratios or minimum times (at end of DATA array) -
            ! these should be preserved as they are.
        if ( TASK == 7 ) then
          do j = 1,NDATA
            if ( DTYPE(j) == 1 ) then
              do m = 1,100000
                k = 1 + int( random(SEED) * (dble(NDATA)-0.000001d0) )
                if ( DTYPE(k) == 1 ) exit
              end do
              INDATA(1,j) = DATA(1,k)
              INDATA(2,j) = DATA(2,k)
              INDATA(3,j) = DATA(3,k)
            else
              INDATA(1,j) = DATA(1,j)
              INDATA(2,j) = DATA(2,j)
              INDATA(3,j) = DATA(3,j)
            end if
          end do
        end if

            ! TASK 8: create a simulated light curve by adding
            ! Gaussian noise to the best-fit model light curve (which
            ! has been put in the DATA array, replacing the actual obs).
        if ( TASK == 8 ) then
          do j = 1,NDATA
            INDATA(2,j) = DATA(2,j) + randomg(SEED,0.0d0,INDATA(3,j))
          end do
        end if

            ! TASK 9: create simulated light curve (not for the times of
            ! minimum and spectroscopic light ratios, which will be at
            ! the end of the dataset in DATA) by taking the model fit
            ! for each datapoint and adding the residual of the best fit
            ! from the datapoint NSIM distant.
        if ( TASK == 9 ) then
          do j = 1,NPHOT
            INDATA(2,j) = DATA(2,j) + RESID(j+i)
          end do
        end if

            ! Now perturb the initial values of the fitted parameters.
            ! Return to the original best-fit parameters without
            ! perturbation if PERTURB=="n", to avoid one dodgy result
            ! screwing up every later simulation.

        do j = 1,67
          if ( (VARY(j)==1 .and. PERTURB=="y") .or. VARY(j)==2 ) then
!             if ( DV(j) > 0.0d0 ) then
!               INV(j) = V(j) + PFACTOR*(RANDOM(SEED)-0.5)*V(j)*DV(j)
!             else if ( DV(j) <= 0.0d0 ) then
!               INV(j) = V(j) + PFACTOR*(RANDOM(SEED)-0.5)*abs(DV(j))
!             end if
            INV(j) = V(j) + (random(SEED)-0.5) * abs(DV(j))
          else
            INV(j) = V(j)
          end if
        end do

        if ( INV(6) > 90.0d0 ) GOTO 500         ! Avoid inclination > 90

            ! Now fit the simulated light curve, then store the results

        CALL FITEBOP (INDATA,DTYPE,NDATA,NLR,NMIN,INV,VARY,LDTYPE,ITER,
     &                CHISQ,VERR,ERROR,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,
     &                NESW,NUMINT,NINTERVAL)

        if ( ERROR /= 0 )  GOTO 500
        if ( INV(3) <= 0.0d0 ) then
         write(6,'(A32,E12.5)')" Retry: V(3) is less than zero: ",INV(3)
          GOTO 500
        end if
        if ( INV(6) < -360.0d0 .or. INV(6) > 720.0d0 ) then
          write(6,'(A25,E12.5)') " Retry: V(6) is haywire: ",INV(6)
          GOTO 500
        end if
        if ( VEXTRA(6) > 0.0d0 .and. CHISQ > 10.0d0*VEXTRA(6) ) then
          write(6,'(A17,E12.6,A1,E12.6,A20)') " CHISQ mismatch (",CHISQ,
     &                              ",",VEXTRA(6),"): reject iteration."
          GOTO 500
        end if

        do j = 1,67
          VALL(j,i) = INV(j)
        end do

        if ( V(2) >= 0.0d0 ) then
          VALLEXTRA(1,i) = INV(2) / (1.0d0 + INV(3))               ! r_1
          VALLEXTRA(2,i) = INV(2) / (1.0d0 + (1.0d0/INV(3)))       ! r_2
        else
          VALL(3,i) = VALL(3,i)
          VALLEXTRA(1,i) = abs(VALL(2,i)) + VALL(3,i)
          VALLEXTRA(2,i) = VALL(3,i) / abs(VALL(2,i))
        end if

        call ECQUADPHASES(V(7),V(8),V(16),ECQPHASES)
        VALLEXTRA(3,i) = GETMODEL (INV,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &        V(20)+V(19)*ECQPHASES(2),2,LP,LS,NUMINT,NINTERVAL) ! L2/L1
        VALLEXTRA(6,i) = CHISQ! / (NDATA - NVARY)

        if ( V(7) < 10.0d0 ) then
          VALLEXTRA(4,i) = sqrt(INV(7)**2 + INV(8)**2)
          VALLEXTRA(5,i) = atan2(INV(8),INV(7))*57.2957795d0
          if (VALLEXTRA(5,i)<0.0d0)VALLEXTRA(5,i)=VALLEXTRA(5,i)+360.0d0
          if(VALLEXTRA(5,i)>=360.d0)VALLEXTRA(5,i)=VALLEXTRA(5,i)-360.d0
        else
          VALLEXTRA(4,i) = (INV(7)-10.0d0) * cos(INV(8)/57.2957795d0)
          VALLEXTRA(5,i) = (INV(7)-10.0d0) * sin(INV(8)/57.2957795d0)
          if ( VALL(8,i) < 0.0d0 ) VALL(8,i) = VALL(8,i) + 360.0d0
          if ( VALL(8,i) >= 360.0d0 ) VALL(8,i) = VALL(8,i) - 360.0d0
        end if

            ! Write result to the big output file. Output the simulation
            ! number to screen if required. Check that a halt has not
            ! been requested. Then finish the iteration loop.

        write (63,'(I5,1X,I3,67(1X,f18.10))')
     &                 i, ITER, (INV(j),j=1,22), (VALLEXTRA(j,i),j=1,6),
     &                (INV(j),j=23,22+3*NSINE), (INV(j),j=38,37+6*NPOLY)

        if ( i < 10 ) write(6,'(1X,I1,$)') i
        if ( abs((dble(i)/10.0d0)-int(dble(i)/10.0d0)) < 0.0001d0 ) then
          if ( i >= 10 .and. i < 100 )   write(6,'(1X,I2,$)') i
          if ( i >= 100 .and. i < 1000 )  write(6,'(1X,I3,$)') i
          if ( i >= 1000 .and. i < 10000 ) write(6,'(1X,I4,$)') i
          if ( i >= 10000 )                 write(6,'(1X,I5,$)') i
        end if

        CALL STOPCHECK (ERROR)
        if ( ERROR /= 0 ) then
          NSIM = i
          exit
        end if
      end do
      write(6,*) " "

!-----------------------------------------------------------------------

            ! TASK 5: pick out the result with the lowest chi-square
            ! value, put its parameters into V, refit the result, and
            ! print the final quantities to the output file.

      if ( TASK == 5 ) then
        j = 0
        HELP1 = 1.0d8
        do i = 1,NSIM
          if ( VALL(43,i) < HELP1 ) then
            j = i
            HELP1 = VALL(43,i)
          end if
        end do

        do i = 1,67
          V(i) = VALL(i,j)
        end do
        CHISQ = HELP1
        CALL FITEBOP(DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &VERR,ERROR,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)

        write(62,'(A44)') "                                            "
        write(62,'(A44)') "-----------------------------------------   "
        write(62,'(A44)') "Overall best fit from all parameter sets:   "
        write(62,'(A44)') "-----------------------------------------   "
        write(62,'(A44)') "                                            "
        CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &    VERR,3,NSINE,PSINE,NPOLY,PPOLY,NL3,NECW,NESW,NUMINT,NINTERVAL)
      end if

            ! TASKS 7,8,9: call OUTPUTSIM to write out 1-sigma uncertai-
            ! nties in each of the fitted and dependent parameters.

      if ( TASK == 7 .or. TASK == 8 .or. TASK == 9 ) then
        write (62,*) " "
        write (62,'(A38,A56)') "======================================",
     &    "==========================================================="
        if ( TASK == 7 ) write (62,'(A13,I6,A27)')
     &          "Results from ",abs(NSIM)," bootstrapping simulations:"
        if ( TASK == 8 ) write (62,'(A13,I6,A25)')
     &            "Results from ",abs(NSIM)," Monte Carlo simulations:"
        if ( TASK == 9 ) write (62,'(A13,I6,A25)')
     &                 "Results from ",NSIM," residual-shifted fits:  "
        write (62,'(A38,A56)') "======================================",
     &    "==========================================================="
        write (62,*) " "

            ! Now print out the 1-sigma results for each variable param.

        write (62,'(A38,A53)') "              Best fit       one_sigma",
     &       "          median         +68.3%         -68.3%         %"
        CALL OUTPUTSIM (V,VARY,VEXTRA,VALL,VALLEXTRA,NSIM,0.683d0)
        write (62,*) " "

            ! Now also print out the 2-sigma results for reference

        if ( int(dble(NSIM)*(1.d0-0.954d0)*0.5d0) >= 1 ) then
        write (62,'(A38,A53)') "              Best fit       two_sigma",
     &       "          median         +68.3%         -68.3%         %"
          CALL OUTPUTSIM (V,VARY,VEXTRA,VALL,VALLEXTRA,NSIM,0.954d0)
          write (62,*) " "
        end if
      end if

      END SUBROUTINE TASK5789
!=======================================================================
!=======================================================================
      SUBROUTINE STOPCHECK (STATUS)  ! Checks if the file "jktebop.stop"
                                     ! exists and if so puts STATUS=1000
      implicit none
      integer ERROR,STATUS

      STATUS = 0
      ERROR = 100
      open (99,file="jktebop.stop",status="old",iostat=ERROR)
      close (99,status="delete")
      if ( ERROR == 0 ) STATUS = 1000

      END SUBROUTINE STOPCHECK
!=======================================================================
      SUBROUTINE OUTPUTSIM (V,VARY,VEXTRA,VALL,VALLEXTRA,NSIM,CONFINT)
            ! Given a set of Monte Carlo or bootstrapping simulation
            ! results, this subroutine calculates the confindence inter-
            ! val for each of the adjusted and dependent parameters.
      implicit none
      real*8 V(67)                  ! IN: The photometric model paramtrs
      integer VARY(67)              ! IN: Which parametrs fixed/adjusted
      real*8 VEXTRA(6)              ! IN: Extra  (dependent)  parameters
      real*8 VALL(67,100000)        ! IN: All the best-fitting parametrs
      real*8 VALLEXTRA(6,100000)    ! IN: All the extra (dependent) pars
      integer NSIM                  ! IN: The number of simulations done
      real*8 CONFINT                ! IN: Fractional confidence interval
      real*8 ARRAY(NSIM)            ! LOCAL: All results for one paramtr
      real*8 HELP1,HELP2,HELP3      ! LOCAL: Useful storage for variabls
      real*8 CONFERR,PERCENT        ! LOCAL: Error estimates for a param
      character*5 NAME(73)          ! LOCAL: Names of the variables etc.
!       character*5 NAMEXTRA(6)       ! LOCAL: Names of the extra paramtrs
      character*5 MESSAGE           ! LOCAL: Warning message
      integer i,j,k                 ! LOCAL: Loop counters
      real*8 SELLECT                ! FUNCTION: Gives nth array value

      data NAME /       "    J","r1+r2","    k","LD_u1","LD_u2","  inc",
     &  "ecosw","esinw","  GD1","  GD2","refl1","refl2","    q","T_lag",
     &  "  L_3"," d(P)","sfact","iring","P_orb","  T_0","LD_n1","LD_n2",
     &          "sin1T","sin1P","sin1A",
     &          "sin2T","sin2P","sin2A",
     &          "sin3T","sin3P","sin3A",
     &          "sin4T","sin4P","sin4A",
     &          "sin5T","sin5P","sin5A",
     &    "p1_co","p1_x ","p1_x2","p1_x3","p1_x4","p1_x5",
     &    "p2_co","p2_x ","p2_x2","p2_x3","p2_x4","p2_x5",
     &    "p3_co","p3_x ","p3_x2","p3_x3","p3_x4","p3_x5",
     &    "p4_co","p4_x ","p4_x2","p4_x3","p4_x4","p4_x5",
     &    "p5_co","p5_x ","p5_x2","p5_x3","p5_x4","p5_x5",
     &     "  r_1","  r_2","L2/L1","    e","omega","Rchi2" /

      if ( V(7) >= 10.0d0 ) then
        NAME(7) = "    e"
        NAME(8) = "omega"
        NAME(71) = "ecosw"
        NAME(72) = "esinw"
!         NAMEXTRA(4) = "ecosw"
!         NAMEXTRA(5) = "esinw"
      end if

      if ( V(2) < 0.0d0 ) then
        NAME(2) = "   r1"
        NAME(3) = "   r2"
        NAME(68) = "r1+r2"
        NAME(69) = "    k"
!         NAMEXTRA(1) = "r1+r2"
!         NAMEXTRA(2) = "    k"
        do i = 1,NSIM
          VALL(2,i) = abs(VALL(2,i))
        end do
        V(2) = abs(V(2))
      end if

            ! Now for each adjusted parameter (and the other calculated
            ! quantities, put the NSIM parameter evaluations into an
            ! array and select the median and 1_sigma bounds (using the
            ! median +/- 0.5_sigma) and output to the parameter file.


      do j = 1,73

        if ( j <= 67 ) then
          do k = 1,NSIM
            ARRAY(k) = mod(VALL(j,k),1.0d3)
          end do
        else
          do k = 1,NSIM
            ARRAY(k) = mod(VALLEXTRA(j-67,k),1.0d3)
          end do
        end if

        HELP1 = SELLECT(ARRAY,NSIM,int(dble(NSIM)*(1.d0-CONFINT)*0.5d0))
        HELP2 = SELLECT(ARRAY,NSIM,int(dble(NSIM)*0.5d0))
!         if ( j <= 67 ) HELP2 = V(j)
!         if ( j > 67 ) HELP2 = VEXTRA(j-67)
        HELP3 = SELLECT(ARRAY,NSIM,int(dble(NSIM)*(1.d0+CONFINT)*0.5d0))

        HELP1 = abs(HELP2-HELP1)        ! Minus CONFINT sigma error
        HELP3 = abs(HELP3-HELP2)        ! Plus  CONFINT sigma error
        MESSAGE = "     "               ! Warn if errors are asymmetric
        if ( HELP1<0.67d0*HELP3.or. HELP1>1.5d0*HELP3) MESSAGE = "skew?"
        if ( HELP1<0.5d0*HELP3 .or. HELP1>2.d0*HELP3 ) MESSAGE = "SKEW!"

        PERCENT = 50.d0 * abs( (abs(HELP1)+abs(HELP3)) / HELP2 )
        CONFERR = 0.5d0 * abs(HELP3+HELP1)

        if ( j <= 37 .and. VARY(j) /= 0 ) then
          if ( j == 1 .or. j == 2 .or. j == 3 ) then
            write (62,100) j,NAME(j),mod(V(j),1.0d3),CONFERR,
     &                      mod(HELP2,1.0d3),HELP3,HELP1,PERCENT,MESSAGE
          else if ( j /= 14 ) then
            write (62,101) j,NAME(j),mod(V(j),1.0d3),CONFERR,
     &                              mod(HELP2,1.0d3),HELP3,HELP1,MESSAGE
          end if
        end if

        if ( j == 68 .or. j == 69 )                       goto 601
        if ( j == 70 .and. V(1) /= 0.0d0 )                goto 601
        if ( j == 71 .and. (VARY(7)/=0 .or. VARY(8)/=0) ) goto 602
        if ( j == 72 .and. (VARY(7)/=0 .or. VARY(8)/=0) ) goto 602
        if ( j == 73 .and. VEXTRA(6) > 0.0d0 )            goto 602
        goto 603

601     write (62,102) NAME(j),VEXTRA(j-67),CONFERR,HELP2,HELP3,HELP1,
     &                                                   PERCENT,MESSAGE
        goto 603
602     write (62,103) NAME(j),VEXTRA(j-67),CONFERR,HELP2,HELP3,HELP1,
     &                                                           MESSAGE
        goto 603
603     continue

      end do

      if ( NAME(38) == "r1+r2" ) V(2) = abs(V(2))


100   FORMAT (I2,1X,A5,1X,F14.10," +/-",F14.10,2X,
     &               F14.10," +",F14.10," -",F14.10," (",F5.2,"%)  ",A5)
101   FORMAT (I2,1X,A5,1X,F14.10," +/-",F14.10,2X,
     &                    F14.10," +",F14.10," -",F14.10,11X,A5)
102   FORMAT    (3X,A5,1X,F14.10," +/-",F14.10,2X,
     &               F14.10," +",F14.10," -",F14.10," (",F5.2,"%)  ",A5)
103   FORMAT    (3X,A5,1X,F14.10," +/-",F14.10,2X,
     &                    F14.10," +",F14.10," -",F14.10,11X,A5)

      END SUBROUTINE OUTPUTSIM
!=======================================================================
      SUBROUTINE OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,
     &                   CHISQ,VERR,WHAT,NSINE,PSINE,NPOLY,PPOLY,NL3,
     &                   NECW,NESW,NUMINT,NINTERVAL)
            ! This subroutine writes information to the output files.
            ! Precisely what to write is given by parameter WHAT.
            ! If WHAT=0, output initial parameters and useful quantities
            ! If WHAT=1, then output the final fitted parameters, useful
            !   quantities, residual curve, and best-fitting model curve
            ! If WHAT=2, output best fit pars, useful stuff and messages
            ! If WHAT=3, output the best fit parameters and useful stuff
            ! If WHAT=101-122 then this has been invoked from TASK6 to
            !   output the name of the parameter=(WHAT-100)
            ! If WHAT=201-222 then this has come from TASK6 to output
            !   fit results for a certain value of parameter=(WHAT-200)
      implicit none
      integer WHAT                        ! IN:  Indicate what to output
      real*8 DATA(3,999999)               ! IN:  Observational data
      integer DTYPE(999999)               ! IN:  Type of each data point
      real*8 V(67),VERR(67)               ! IN:  Parameters  and  errors
      integer VARY(67)                    ! IN:  Par adjustment integers
      integer NDATA,NLR,NMIN              ! IN:  Numbers of  data points
      integer NSINE,NL3,NECW,NESW         ! IN:  Number of sines and L3s
      integer PSINE(5)                    ! IN:  Which par for each sine
      integer NPOLY,PPOLY(5)              ! IN:  Similar for polynomials
      integer NUMINT                      ! IN: Number of numerical ints
      real*8  NINTERVAL                   ! IN: Time interval for numint
      integer LDTYPE(2)                   ! IN:  LDlaw type foreach star
      integer ITER                        ! IN:  numbr of fit iterations
      real*8 CHISQ                        ! IN:   chi-squared of the fit
      integer NVARY,VWHERE(67)            ! LOCAL: which parameters vary
      character VNAME(67)*19,VNAM(67)*5   ! LOCAL: names  of  parameters
      real*8 ECC,OMEGA,ECOSW,ESINW        ! LOCAL: orbital    parameters
      real*8 R1,R2,R12,RRAT               ! LOCAL: fractional star radii
      real*8 RESIDSUM,RESIDABSSUM         ! LOCAL: sums of the residuals
      real*8 RESIDSQSUM,SE                ! LOCAL: resid, standard error
      real*8 LTOTAL                       ! LOCAL: total light of system
      real*8 A,B,EPSLN1,EPSLN2,EPSLN      ! LOCAL: stellar shape  params
      real*8 MAG,LP,LS,HJD                ! LOCAL: LIGHT subroutine pars
      real*8 PHASE,HELP,HELP2             ! LOCAL: some useful variables
      real*8 OMC(999999)                  ! LOCAL: observed minus calc'd
      real*8 LIMBRIGHT                    ! LOCAL: total limb brightness
      real*8 UNINTMAG                     ! LOCAL: unintegrated magntude
      integer NPHASE                      ! LOCAL: number phases to plot
      character MESSAGE1*19, MESSAGE2*17  ! LOCAL: mesages to go to file
      real*8 ECQPHASES(4)                 ! LOCAL: important orbitphases
      integer i,j,k,ERROR,STATUS          ! LOCAL: counters & errorflags
      real*8 GETPHASE                     ! FUNCTION: calc orbital phase
      real*8 GETMODEL                     ! FUNCTION: calcs model output

      DATA VNAME      /     "Surf. bright. ratio","Sum of the radii   ",
     &"Ratio of the radii ","Limb darkening u_1 ","Limb darkening u_2 ",
     &"Orbit inclination  ","ecc * cos(omega)   ","ecc * sin(omega)   ",
     &"Grav darkening 1   ","Grav darkening 2   ","Reflected light 1  ",
     &"Reflected light 2  ","Mass ratio (q)     ","Tide lead/lag angle",
     &"Third light (L_3)  ","Phase correction   ","Light scale factor ",
     &"Integration ring   ","Orbital period (P) ","Ephemeris timebase ",
     &"Limb darkening n_1 ","Limb darkening n_2 ","Sine 1 timebase    ",
     &"Sine 1 period      ","Sine 1 amplitude   ","Sine 2 timebase    ",
     &"Sine 2 period      ","Sine 2 amplitude   ","Sine 3 timebase    ",
     &"Sine 3 period      ","Sine 3 amplitude   ","Sine 4 timebase    ",
     &"Sine 4 period      ","Sine 4 amplitude   ","Sine 5 timebase    ",
     &"Sine 5 period      ","Sine 5 amplitude   ","Poly 1 pivot       ",
     &"Poly 1 coeff of x  ","Poly 1 coeff of x^2","Poly 1 coeff of x^3",
     &"Poly 1 coeff of x^4","Poly 1 coeff of x^5","Poly 2 pivot       ",
     &"Poly 2 coeff of x  ","Poly 2 coeff of x^2","Poly 2 coeff of x^3",
     &"Poly 2 coeff of x^4","Poly 2 coeff of x^5","Poly 3 pivot       ",
     &"Poly 3 coeff of x  ","Poly 3 coeff of x^2","Poly 3 coeff of x^3",
     &"Poly 3 coeff of x^4","Poly 3 coeff of x^5","Poly 4 pivot       ",
     &"Poly 4 coeff of x  ","Poly 4 coeff of x^2","Poly 4 coeff of x^3",
     &"Poly 4 coeff of x^4","Poly 4 coeff of x^5","Poly 5 pivot       ",
     &"Poly 5 coeff of x  ","Poly 5 coeff of x^2","Poly 5 coeff of x^3",
     &"Poly 5 coeff of x^4","Poly 5 coeff of x^5" /

      data VNAM/"    J","r1+r2","    k","LD_u1","LD_u2","  inc","ecosw",
     &          "esinw","  GD1","  GD2","refl1","refl2","    q","T_lag",
     &          "  L_3"," d(P)","sfact","iring","P_orb","  T_0","LD_n1",
     &  "LD_n2","sin1T","sin1P","sin1A","sin2T","sin2P","sin2A","sin3T",
     &  "sin3P","sin3A","sin4T","sin4P","sin4A","sin5T","sin5P","sin5A",
     &    "p1_co","p1_x ","p1_x2","p1_x3","p1_x4","p1_x5",
     &    "p2_co","p2_x ","p2_x2","p2_x3","p2_x4","p2_x5",
     &    "p3_co","p3_x ","p3_x2","p3_x3","p3_x4","p3_x5",
     &    "p4_co","p4_x ","p4_x2","p4_x3","p4_x4","p4_x5",
     &    "p5_co","p5_x ","p5_x2","p5_x3","p5_x4","p5_x5" /


            ! Find the stellar radii for both options (r1,r2), (r1+r2,k)

      R1 = abs(V(2)) / (1.0d0 + V(3))
      R2 = abs(V(2)) / (1.0d0 + (1.0d0/V(3)))
      R12 = abs(V(2))
      RRAT = V(3)

      if ( V(2) < 0.0d0 ) then
        VNAME(2) = "Fractional radius 1"
        VNAME(3) = "Fractional radius 2"
        VNAM(2) = "  r_1"
        VNAM(3) = "  r_2"
        R1 = abs(V(2))
        R2 = V(3)
        R12 = abs(V(2)) + V(3)
        RRAT = V(3) / abs(V(2))
      end if

            ! Find the light contributions and reflection coefficients
            ! for the two stars. Do this at phase of first quadrature.

      call ECQUADPHASES(V(7),V(8),V(16),ECQPHASES)
      MAG = GETMODEL (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                V(20)+V(19)*ECQPHASES(2),1,LP,LS,NUMINT,NINTERVAL)
!       LTOTAL = LP + LS + V(15)
      if ( VARY(11) == -1 )  V(11) = 0.4d0 * (LS/(1.0d0-V(15))) * R1**2
      if ( VARY(12) == -1 )  V(12) = 0.4d0 * (LP/(1.0d0-V(15))) * R2**2

            ! Find the number of variable parameters

      NVARY = 0
      do i = 1,67
        if ( VARY(i) == 1 .or. VARY(i) == 3 ) NVARY = NVARY + 1
      end do

            ! Output the values of the parameters and various other info

      if ( WHAT == 0 ) then
        write(62,107) "Initial values of the parameters:               "
      else if ( WHAT == 1 .or. WHAT == 2 ) then
        write(62,*)   " "
        write(62,107) "---------------------------------------         "
        write(62,110) ITER," iterations of EBOP completed     "
        write(62,107) "---------------------------------------         "
        write(62,*)   " "
        write(62,105)   "Warning: the uncertainties given below a",
     &                  "re only formal errors and are probably a"
        write(62,105)   "bit optimistic.  Use TASKs 6, 7 and 8 to",
     &                  " get reliable errors for the parameters."
        write(62,105)   "They also assume the noise is only Poiss",
     &                  "onian: use TASK9 if red noise is strong."
        if ( DATA(3,1) < 0.0d0 ) then
          write(62,105) "No observational errors; so parameter er",
     &                  "rors have been scaled by sqrt(chi^2_red)"
        end if
        write(62,*)   " "
        if ( V(6) >= 89.9d0 ) then
          write(62,105) "### WARNING: the light curve solution ha",
     &                  "s an inclination of 89.9 degree or more."
          write(62,105) "The minimisation algorithm cannot go bey",
     &                  "ond 90 degrees, so the solution has been"
          write(62,105) "kept artificially slightly below this va",
     &                  "lue to preserve its numerical stability."
          write (62,*) " "
        end if
        write(62,107) "Final values of the parameters:                 "
      end if

!-----------------------------------------------------------------------

      if ( WHAT == 0 .or. WHAT == 1 .or. WHAT == 2 .or. WHAT == 3 ) then
        if ( V(6) > 90.0d0 )  V(6) = 180.0d0 - V(6)
        if ( V(6) < 0.0d0 ) V(6) = V(6) + 180.0d0
        if ( LDTYPE(2) == 0 ) then
          V(5) = V(4)
          VERR(5) = VERR(4)
        end if

        do i = 1,67           ! 22+3*NSINE
          if ( i == 7 .and. V(7) > 5.0d0 ) then
            MESSAGE1 = "Eccentricity       "
          else if ( i == 8 .and. V(7) > 5.0d0 ) then
            MESSAGE1 = "Periastronlongitude"
          else
            MESSAGE1 = VNAME(i)
          end if

          MESSAGE2 = "                 "
          if ( WHAT == 0 .and. (VARY(i) == 1 .or. VARY(i) == 3) )
     &                                    MESSAGE2 = "  (adjusted)     "
          if ( WHAT == 1 .and. (VARY(i) == 0 .or. VARY(i) == 2) )
     &                                    MESSAGE2 = "   (fixed)       "
          if ( i == 11 .and. VARY(i) == -1) MESSAGE2="  (from geometry)"
          if ( i == 12 .and. VARY(i) == -1) MESSAGE2="  (from geometry)"

          if ( i <= 22+3*NSINE .or. (i>=38 .and. i<=37+6*NPOLY) ) then
            if ( WHAT == 0 ) then
              write(62,113) i,MESSAGE1,V(i),MESSAGE2
            else
              if ( VARY(i) /= 1 .and. VARY(i) /= 3 ) then
                write(62,113) i,MESSAGE1,V(i),MESSAGE2
              else
                write(62,114) i,MESSAGE1,V(i),VERR(i)
              end if
            end if
          end if

        end do

        write (62,*) " "
        if ( LDTYPE(1) == 1 ) write (62,116) "linear          "
        if ( LDTYPE(1) == 2 ) write (62,116) "logarithmic     "
        if ( LDTYPE(1) == 3 ) write (62,116) "square-root     "
        if ( LDTYPE(1) == 4 ) write (62,116) "quadratic       "
        if ( LDTYPE(1) == 5 ) write (62,116) "cubic           "
        if ( LDTYPE(2) == 0 ) write (62,117) "same as primary "
        if ( LDTYPE(2) == 1 ) write (62,117) "linear          "
        if ( LDTYPE(2) == 2 ) write (62,117) "logarithmic     "
        if ( LDTYPE(2) == 3 ) write (62,117) "square-root     "
        if ( LDTYPE(2) == 4 ) write (62,117) "quadratic       "
        if ( LDTYPE(2) == 5 ) write (62,116) "cubic           "

            ! Now check to see if the limb darkening coefficients are
            ! physically possible or reasonable, and print a warning
            ! if they are not. Remember that EBOP can happily function
            ! with physically impossible parameters (e.g. LDlin > 1.0)

        if ( LDTYPE(1) == 1 ) then
          LIMBRIGHT = 1.0d0 - V(4)
        else
          LIMBRIGHT = 1.0d0 - V(4) - V(21)
        end if
        if ( LIMBRIGHT > 1.0d0 ) then
          write (62,*) " "
          write (62,105) "### WARNING: the total limb darkening at",
     &                   " the limb of star A is less than 0.0, so"
          write (62,105) "### the limb darkening coefficient(s) fo",
     &                   "r this star are physically improbable.  "
        else if ( LIMBRIGHT < 0.0d0 ) then
          write (62,*) " "
          write (62,105) "### WARNING: the total limb darkening at",
     &                   " the limb of star A is greater than 1.0,"
          write (62,105) "### so the limb darkening coefficient(s)",
     &                   " for this star are physically impossible"
        end if

        if ( LDTYPE(2) == 1 ) then
          LIMBRIGHT = 1.0d0 - V(5)
        else
          LIMBRIGHT = 1.0d0 - V(5) - V(22)
        end if
        if ( LIMBRIGHT < 0.0d0 ) then
          write (62,105) "### WARNING: the total limb darkening at",
     &                   " the limb of star B is less than 0.0, so"
          write (62,105) "### the limb darkening coefficient(s) fo",
     &                   "r this star are physically impossible.  "
        else if ( LIMBRIGHT > 1.0d0 ) then
          write (62,105) "### WARNING: the total limb darkening at",
     &                   " the limb of star B is more than 1.0, so"
          write (62,105) "### the limb darkening coefficient(s) fo",
     &                   "r this star are physically impossible.  "
        end if

        if ( WHAT == 0 .and. (LDTYPE(1) == 2 .or. LDTYPE(2) == 2 .or.
     &                       LDTYPE(1) == 5 .or. LDTYPE(2) == 5 ) ) then
          write (62,*) " "
          write (62,105) "### WARNING: the logarithmic and cubic l",
     &                   "imb darkening law flux normalisation are"
          write (62,105) "only approximate. Do not use the laws if",
     &                   " you need a precise (<few%) light ratio."
        end if

        write (62,*) " "

        write(62,102)"Phase of primary eclipse:           ",ECQPHASES(1)
        write(62,102)"Phase of first quadrature:          ",ECQPHASES(2)
        write(62,102)"Phase of secondary eclipse:         ",ECQPHASES(3)
        write(62,102)"Phase of second quadrature:         ",ECQPHASES(4)

        write (62,*) " "

            ! Print radii of the stars and check if they are too large
            ! for EBOP to provide a good approximation to their shapes.

        if ( V(2) >= 0.0d0 ) then
          write (62,102) "Fractional primary radius:          ",R1
          write (62,102) "Fractional secondary radius:        ",R2
        else
          write (62,102) "Sum of the fractional radii:        ",R12
          write (62,102) "Ratio of the radii:                 ",RRAT
        end if

        if ( R1 >= 0.6d0 ) then
          write (62,*) " "
          write (62,105) "### WARNING: average radius is greater t",
     &                   "han 0.3 so the radii may be wrong by 5%!"
          write (62,105) "### See North & Zahn (2004, New.Ast.Rev.",
     &                   ", 47, 741) for details. Do not use EBOP!"
          write (62,*) " "
        else if ( R2 >= 0.5d0 )  then
          write (62,*) " "
          write (62,105) "### WARNING: average radius is greater t",
     &                   "han 0.25 so the radii may be wrong by 1%"
          write (62,105) "### See North & Zahn (2004, New Astronom",
     &                   "y Review, 47, 741) for further details. "
          write (62,*) " "
        end if

        write (62,'(A26,F7.4,A3,2X,F17.10)')
     &         "Stellar light ratio (phase",ECQPHASES(2),"): ",LS/LP
        write (62,102) "Primary contribut'n to system light:", LP
        write (62,102) "Secondary contrib'n to system light:", LS
        write (62,*) " "

        if ( V(7) > 5.0d0 ) then
          ECC = V(7) - 10.0d0
          OMEGA = V(8) * 45.0d0 / atan(1.0d0)
          ECOSW = (V(7)-10.0d0) * cos(V(8)/57.2957795d0)
          ESINW = (V(7)-10.0d0) * sin(V(8)/57.2957795d0)
          write (62,102) "Eccentricity * cos(omega):          ",ECOSW
          write (62,102) "Eccentricity * sin(omega):          ",ESINW
        else if ( V(7) /= 0.0d0 .and. V(8) /= 0.0d0 ) then
          ECC = sqrt(V(7)**2 + V(8)**2)
          OMEGA = atan2(V(8),V(7)) * 45.0d0 / atan(1.0d0)
          if ( OMEGA < 0.0d0 ) OMEGA = OMEGA + 360.0d0
          if ( OMEGA > 360.0d0 ) OMEGA = OMEGA - 360.0d0
          write (62,102) "Eccentricity:                       ",ECC
          write (62,102) "Omega (degrees):                    ",OMEGA
        end if

        if ( VARY(11) /= 1 .or. VARY(12) /= -1 ) then
          HELP = 0.4d0 * LS * R1**2
          write (62,102) "Geometric reflection coeff (star A):",HELP
          HELP = 0.4d0 * LP * R2**2
          write (62,102) "Geometric reflection coeff (star B):",HELP
        end if

        CALL BIAX (R1,abs(V(13)),A,B,EPSLN1)
        CALL BIAX (R2,1.0d0/abs(V(13)),A,B,EPSLN2)
        if ( V(13) <= 0.0d0 ) EPSLN1 = 0.0d0
        if ( V(13) <= 0.0d0 ) EPSLN2 = 0.0d0

        write (62,102) "Oblateness of the primary star:     ",EPSLN1
        write (62,102) "Oblateness of the secondary star:   ",EPSLN2
        if ( (WHAT == 1 .or. WHAT == 2) .and. V(13) < -0.000001d0 ) then
          CALL BIAX(R1,abs(V(13)),A,B,EPSLN)
          write(62,102)"Expected oblateness of primary:     ",EPSLN
          CALL BIAX(R2,1.0d0/abs(V(13)),A,B,EPSLN)
          write(62,102)"Expected oblateness of secondary:   ",EPSLN
        end if

        if ( EPSLN1 > 0.04d0 .or. EPSLN2 > 0.04d0 )  then
          write (62,*) " "
          write (62,105) "### WARNING: oblateness is above the rec",
     &                   "ommended maximum value for EBOP of 0.04."
          write (62,105) "See Popper & Etzel (1981, AJ, 86, 102) f",
     &                   "or a justification and further details. "
        end if
        write (62,*) " "
      end if

!-----------------------------------------------------------------------
            ! Now output the observations and the calculated residuals.

      if ( WHAT == 1 .or. WHAT == 2 .or. WHAT == 3 ) then
        if ( WHAT == 1 ) write (63,'(A21,A49)') "#     TIME       MAGN",
     &              "ITUDE    ERROR      PHASE        MODEL      (O-C)"

        RESIDSUM = 0.0d0
        RESIDABSSUM = 0.0d0
        RESIDSQSUM = 0.0d0
        CHISQ = 0.0d0
        do i = 1,NDATA
          if ( DTYPE(i) == 1 ) then
            MAG = GETMODEL(V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,DATA(1,i),
     &                                  DTYPE(i),LP,LS,NUMINT,NINTERVAL)
            OMC(i) = DATA(2,i)-MAG
            if ( WHAT == 1 ) write(63,104)DATA(1,i),DATA(2,i),DATA(3,i),
     &                        GETPHASE(DATA(1,i),V(19),V(20)),MAG,OMC(i)
            RESIDSUM = RESIDSUM + OMC(i)
            RESIDABSSUM = RESIDABSSUM + abs(OMC(i))
            RESIDSQSUM = RESIDSQSUM + OMC(i)**2
            CHISQ = CHISQ + (OMC(i) / DATA(3,i))**2
          end if
        end do

        write(62,102) "Sum of residuals for LC data (mag): ",   RESIDSUM
        write(62,102) "Sum of absolute values of LC resids:",RESIDABSSUM
        write(62,102) "Sum of the squares of the LC resids:", RESIDSQSUM
        write(62,115) "Total datapoints / fitting params:  ",NDATA,NVARY

        if ( WHAT == 1 .or. WHAT == 2 .or. WHAT == 3 ) then
          write (62,102) "Standard error (one LC obs, mmag):  ",
     &              sqrt(RESIDSQSUM / (NDATA-NVARY-NLR-NMIN)) * 1000.0d0
          write (62,102) "rms of residuals (one LC obs, mmag):",
     &                    sqrt(RESIDSQSUM / (NDATA-NLR-NMIN)) * 1000.0d0
          if ( DATA(3,1) > 0.0d0 ) then
            write (62,102)  "Chi-squared from errorbars:         ",CHISQ
            write (62,1022) "Number of degrees of freedom:       ",
     &                                                       NDATA-NVARY
            write (62,102)  "Reduced chi-squared from errorbars: ",
     &                                               CHISQ/(NDATA-NVARY)
          end if
        end if
        write (62,*) " "
      end if

!-----------------------------------------------------------------------
            ! Now output the best fit on a grid of orbital phases

      if ( WHAT == 1 ) then
        NPHASE = 10001
        if ( R1 < 0.01d0 .or. R2 < 0.01d0 ) NPHASE = 100001
        if ( NUMINT == 1 ) write (64,'(A48)')
     &                "# PHASE   MAGNITUDE     L1         L2         L3"
        if ( NUMINT /= 1 ) write (64,'(A63)')
     &  "# PHASE   MAGNITUDE     L1         L2         L3      UNINTMAG"
        do i = 1,NPHASE
          PHASE = (i-1) / dble(NPHASE-1)
          HJD = V(20) + PHASE * V(19)
          MAG = GETMODEL(V,LDTYPE,0,PSINE,0,PPOLY,HJD,1,LP,LS,
     &                                                 NUMINT,NINTERVAL)
          if ( NUMINT /= 1 ) UNINTMAG = GETMODEL(V,LDTYPE,0,PSINE,0,
     &                                    PPOLY,HJD,1,LP,LS,1,NINTERVAL)
          if ( NUMINT == 1 ) write(64,12) PHASE,MAG,LP,LS,V(15)
          if ( NUMINT /= 1 ) write(64,12) PHASE,MAG,LP,LS,V(15),UNINTMAG
12        FORMAT (F8.6,5(1X,F10.6))
        end do
      end if

!-----------------------------------------------------------------------
      ! Now output the results and residuals forthe direct observational
      ! constraints on the times of minimum, spectroscopic light ratios,
      ! third light, ecosw and esinw.

      if ( WHAT == 1 .or. WHAT == 2 .or. WHAT == 3  ) then
        if ( NLR>0 .or. NMIN>0 .or. NL3>0 .or. NECW>0 .or. NESW>0) then
          write (62,'(A36,A36)') "------------------------------------",
     &                           "-------------------------------     "
          write (62,'(A36,A36)') "Results for the additional observed ",
     &                           "quantities                          "
          write (62,'(A36,A36)') "------------------------------------",
     &                           "-------------------------------     "
          write (62,'(A36,A36)') "Type      Time/Cycle     Measurement",
     &                           "     Error      Model     Sigma     "
          do i = 1,NDATA
            if ( DTYPE(i) >= 2 .and. DTYPE(i) <= 6 ) then
              HELP = GETMODEL(V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                        DATA(1,i),DTYPE(i),LP,LS,NUMINT,NINTERVAL)
!               HELP2 = DATA(1,i)
            end if

            if ( DTYPE(i) == 2 ) write (62,120) "LRAT", DATA(1,i),
     &               DATA(2,i),DATA(3,i),HELP,(HELP-DATA(2,i))/DATA(3,i)
            if ( DTYPE(i) == 3 ) write (62,120) "TMIN", DATA(1,i),
     &               DATA(2,i),DATA(3,i),HELP,(HELP-DATA(2,i))/DATA(3,i)
            if ( DTYPE(i) == 4 ) write (62,120) "THDL", DATA(1,i),
     &               DATA(2,i),DATA(3,i),HELP,(HELP-DATA(2,i))/DATA(3,i)

            if ( DTYPE(i) == 5 ) then
              if ( V(7) > 5.0d0 ) write (62,120) "ECC ", DATA(1,i),
     &                            DATA(2,i)-10.d0,DATA(3,i),HELP-10.0d0,
     &                                 (HELP-DATA(2,i))/DATA(3,i)
              if ( V(7) < 5.0d0 ) write (62,120) "ECSW", DATA(1,i),
     &               DATA(2,i),DATA(3,i),HELP,(HELP-DATA(2,i))/DATA(3,i)
            end if

            if ( DTYPE(i) == 6 ) then
              if ( V(7) > 5.0d0 ) write (62,120) "OMGA", DATA(1,i),
     &               DATA(2,i),DATA(3,i),HELP,(HELP-DATA(2,i))/DATA(3,i)
              if ( V(7) < 5.0d0 ) write (62,120) "ESNW", DATA(1,i),
     &               DATA(2,i),DATA(3,i),HELP,(HELP-DATA(2,i))/DATA(3,i)
            end if

!             if ( DTYPE(i) == 5 ) write (62,120) "ECSW", DATA(1,i),
!      &               DATA(2,i),DATA(3,i),HELP,(HELP-DATA(2,i))/DATA(3,i)
!             if ( DTYPE(i) == 6 ) write (62,120) "ESNW", DATA(1,i),
!      &               DATA(2,i),DATA(3,i),HELP,(HELP-DATA(2,i))/DATA(3,i)
          end do
          write (62,'(A36,A36)') "------------------------------------",
     &                           "-------------------------------     "
          write (62,*) " "
        end if

      end if

!--------------------------TASK6----------------------------------------
            ! If WHAT is between 100 and 122 then output the name of the
            ! parameter which is to be investigated next.
            ! If WHAT is between 200 and 222 then output various results
            ! from the last EBOP fit on one line for comparison.

      if ( WHAT > 100 .and. WHAT < 167 ) then
        j = 0
        do i = 1,67
          if ( VARY(i) == 1 ) then
            j = j + 1
            VWHERE(j) = i
          end if
        end do
        write (62,*) " "
        write (62,'(A41,A20)')
     &       "Now investigating the fitting parameter: ",VNAME(WHAT-100)
        write (62,'(A4,3X,37(A5,7X))')    "ITER",VNAM(WHAT-100)," rms ",
     &                   (VNAM(VWHERE(k)),k=1,j)," r_1 "," r_2 ","L2/L1"
        write(6,'(A44,A20)')
     &   ">> Now investigating the fitting parameter:  ",VNAME(WHAT-100)
      end if

!-----------------------------------------------------------------------

      if ( WHAT > 200 .and. WHAT < 267 ) then
        RESIDSQSUM = 0.0d0
        do i = 1,NDATA
          MAG = GETMODEL(V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,DATA(1,i),1,
     &                                           LP,LS,NUMINT,NINTERVAL)
          RESIDSQSUM = RESIDSQSUM + (DATA(2,i)-MAG)**2
        end do
        SE = sqrt(RESIDSQSUM / NDATA) * 1000.0d0
        j = 0
        do i = 1,67
          if ( VARY(i) == 1 ) then
            j = j + 1
            VWHERE(j) = i
          end if
        end do
       call ECQUADPHASES(V(7),V(8),V(16),ECQPHASES)
       write(62,'(I3,23(1X,F10.7))')   ITER, mod(V(WHAT-200),1.0d3), SE,
     &                          (mod(V(VWHERE(k)),1.0d3),k=1,j), R1, R2,
     &                        GETMODEL(V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,
     &                V(20)+V(19)*ECQPHASES(2),2,LP,LS,NUMINT,NINTERVAL)
      end if
!-----------------------------------------------------------------------

! 100   FORMAT (I2,2X,A19,2X,F18.10)
! 101   FORMAT (I2,2X,A19,2X,F18.10,A12)
102   FORMAT (A36,2X,F17.10)
! 1021  FORMAT (A36,2X,E17.10)
1022  FORMAT (A36,2X,I17)
104   FORMAT (F14.6,1X,F11.6,1X,F10.6,3X,F10.8,1X,F10.6,1X,F10.6)
105   FORMAT (A40,A40)
107   FORMAT (A48)
! 108   FORMAT (F14.6,F10.6)
110   FORMAT (I6,A34)
! 111   FORMAT (I2,2X,A18,2X,F18.10,A17)
! 112   FORMAT (A39,F7.4)
113   FORMAT (I2,2X,A19,2X,F18.10,A17)
114   FORMAT (I2,2X,A19,2X,F18.10," +/- ",F14.10)
115   FORMAT (A36,7X,I6,2X,I4)
116   FORMAT ("Limb darkening law for primary star:       ",A16)
117   FORMAT ("Limb darkening law for secondary star:     ",A16)
118   FORMAT (F8.1,1X,F13.5,1X,F8.5,2X,F13.5,1X,F6.2)
119   FORMAT (F13.5,1X,F7.4,1X,F6.4,2X,F7.4,F7.2)
120   FORMAT (A4,3X,F14.6,1X,F14.6,1X,F10.6,1X,F12.6,F7.2)

      END SUBROUTINE OUTPUT
!=======================================================================
      SUBROUTINE GET_DV (V,DV)
            ! This subroutine produces the DV values - the intervals for
            ! evaluating the numerical derivatives for each variable.
      implicit none
      real*8 V(67)                        ! IN:   Photometric parameters
      real*8 DV(67)                       ! OUT:   DV values
      real*8 DVD(67)                      ! LOCAL:  basic DV values
      integer i                           ! LOCAL:  loop counter

      data DVD/ 0.01d0,   0.01d0,  0.01d0,  -0.05d0,   ! J rs k u1
     &         -0.05d0,  -0.1d0,  -0.001d0, -0.001d0,  ! u2 i e w
     &          0.1d0,    0.1d0,   0.01d0,   0.01d0,   ! GD1,2 ref1,2
     &          0.01d0,  -1.0d0,  -0.01d0,  -0.0001d0, ! q tangl L3 dPhi
     &         -0.01d0,   0.5d0,   0.1d-6,             ! sfact intring P
     &         -0.001d0, -0.05d0, -0.05d0,             ! T0 v1 v2
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 1 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 2 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 3 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 4 T0,p,amp
     &       -0.01d0,   0.0001d0,   0.01d0,            ! sine 5 T0,p,amp
     & -0.0001d0,-0.01d0,-0.008d0,-0.005d0,-0.005d0,-0.002d0,   ! poly 1
     & -0.0001d0,-0.01d0,-0.008d0,-0.005d0,-0.005d0,-0.002d0,   ! poly 2
     & -0.0001d0,-0.01d0,-0.008d0,-0.005d0,-0.005d0,-0.002d0,   ! poly 3
     & -0.0001d0,-0.01d0,-0.008d0,-0.005d0,-0.005d0,-0.002d0,   ! poly 4
     & -0.0001d0,-0.01d0,-0.008d0,-0.005d0,-0.005d0,-0.002d0 /  ! poly 5

            ! The DVD values are the basic information from which the DV
            ! values are calculated, and essentially represent the order
            ! of magnitude expected for each parameter (based on experi-
            ! ence) and the way in which the DV values should be done.
            ! If DVD(i) > 0 then is it a multiplicative perturbation.
            ! If DVD(i) < 0 then is it an additive perturbation.

      do i = 1,67

        if ( DVD(i) >= 0.0d0 ) then
          DV(i) = abs(DVD(i)*V(i))
        else
          DV(i) = abs(DVD(i))
        end if

            ! The sine timesbases need special treatment: set them to
            ! one hundredth of the period of the sine wave in question.

        if ( i==23 .or. i==26 .or. i==29 .or. i==32 .or. i==35 ) then
          DV(i) = V(i+1) / 100.0d0
        end if

      end do

      END SUBROUTINE GET_DV
!=======================================================================
      DOUBLEPRECISION FUNCTION GETMODEL (V,LDTYPE,NSINE,PSINE,NPOLY,
     &                          PPOLY,TIME,DTYPE,LA,LB,NUMINT,NINTERVAL)
            ! Output a predicted model value according to the parameters
            ! in array V. Precise meaning of the value depends on DTYPE.
            ! DTYPE=1  it outputs an EBOP magnitude for given time
            ! DTYPE=2  it outputs a light ratio for the given time
            ! DTYPE=3  outputs a time of eclipse for the given =CYCLE=
            ! DTYPE=4  it simply outputs the third light value
      implicit none
      real*8 V(67)                  ! IN: Photometric parameters
      integer LDTYPE(2)             ! IN: LD law type for the two stars
      real*8 TIME                   ! IN: The given TIME, PHASE or CYCLE
      integer DTYPE                 ! IN: 1-6 depending on wanted result
      integer NSINE,PSINE(5)        ! IN: number and parameters of sines
      integer NPOLY,PPOLY(5)        ! IN: number and parameters of polys
      integer NUMINT                ! IN: Number of numerical integratns
      real*8 NINTERVAL              ! IN: Time interval for integrations
      real*8 LA,LB                  ! OUT: Light produced by each star
      real*8 FMAG,LP,LS             ! LOCAL: LIGHT subroutine output
      real*8 ECC,OMEGA,ECOSW,ESINW  ! LOCAL: orbital shape parameters
      real*8 GETMIN,GETPHASE        ! FUNCTIONS
      real*8 FMAGSUM,LASUM,LBSUM
      integer i
      real*8 TIMEIN

      GETMODEL = 0.0d0

            ! DTYPE=1 for light curve datapoint
            ! DTYPE=2 for light ratio
            ! DTYPE=3 for times of minimum light
            ! DTYPE=4 for third light
            ! DTYPE=5 for e*cos(omega)
            ! DTYPE=6 for e*sin(omega)

      if ( DTYPE == 1 ) then
        if ( NUMINT == 1 ) then
         CALL LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIME,FMAG,LP,LS)
          LA = LP
          LB = LS
          GETMODEL = FMAG

        else if ( NUMINT > 1 ) then
          FMAGSUM = 0.0d0
          LASUM = 0.0d0
          LBSUM = 0.0d0
          CALL LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIME,FMAG,LP,LS)
          do i = 1,NUMINT
            TIMEIN  =  TIME  -  NINTERVAL / 86400.0d0 / dble(NUMINT) *
     &        ( dble(NUMINT) - 2.0d0*dble(i) + 1.0d0 ) / 2.0d0
          CALL LIGHT(V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIMEIN,FMAG,LP,LS)
            FMAGSUM = FMAGSUM + FMAG
            LASUM = LASUM + LP
            LBSUM = LBSUM + LS
          end do
          FMAG = FMAGSUM / NUMINT
          LA = LASUM / NUMINT
          LB = LBSUM / NUMINT
          GETMODEL = FMAG
        else
          write(6,*)"NUMINT is less than 1 in function GETMODEL. Abort."
          write(6,*)"NUMINT =    ", NUMINT
          write(6,*)"NINTERVAL = ", NINTERVAL
          stop
        end if

      else if ( DTYPE == 2 ) then
        CALL LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIME,FMAG,LP,LS)
        GETMODEL = LS / LP

      else if ( DTYPE == 3 ) then
        if ( V(7) > 5.0d0 ) then
          ECC = V(7) - 10.0d0
          OMEGA = V(8)
        else
          ECC = sqrt(V(7)**2 + V(8)**2)
          OMEGA = atan2(V(8),V(7)) * 45.0d0 / atan(1.0d0)
          if ( OMEGA < 0.0d0 ) OMEGA = OMEGA + 360.0d0
        end if
        GETMODEL = GETMIN (V(20),V(19),ECC,OMEGA,TIME)

      else if ( DTYPE == 4 ) then
        GETMODEL = V(15)

      else if ( DTYPE == 5 .or. DTYPE == 6 ) then
        if ( V(7) > 5.0d0 ) then
          ECC = V(7)
          OMEGA = V(8)
          ECOSW = (V(7)-10.0d0) * cos(V(8)/57.2957795d0)
          ESINW = (V(7)-10.0d0) * sin(V(8)/57.2957795d0)
          if ( DTYPE == 5 )  GETMODEL = V(7)!ECC
          if ( DTYPE == 6 )  GETMODEL = V(8)!OMEGA
        else
          ECOSW = V(7)
          ESINW = V(8)
          ECC = sqrt(V(7)**2 + V(8)**2)
          OMEGA = atan2(V(8),V(7)) * 45.0d0 / atan(1.0d0)
          if ( OMEGA < 0.0d0 ) OMEGA = OMEGA + 360.0d0
          if ( OMEGA > 360.0d0 ) OMEGA = OMEGA - 360.0d0
          if ( DTYPE == 5 )  GETMODEL = ECOSW
          if ( DTYPE == 6 )  GETMODEL = ESINW
        end if

      else
        GETMODEL = -100.0d0
        write(6,*) "### ERROR: wrong datatype asked for in GETMODEL: ",
     &              DTYPE
        STOP
      end if

      END FUNCTION GETMODEL
!=======================================================================
      DOUBLEPRECISION FUNCTION GETPHASE (HJD,PERIOD,TZERO)
            ! Returns phase from given time and orbital ephemeris
      implicit none
      real*8 HJD,PERIOD,TZERO

      GETPHASE = (HJD - TZERO) / PERIOD
      GETPHASE = GETPHASE - int(GETPHASE)
      if ( GETPHASE < 0.0d0 ) GETPHASE = GETPHASE + 1.0d0

      END FUNCTION GETPHASE
!=======================================================================
      SUBROUTINE ECQUADPHASES (ECCIN,OMEGAIN,PSHIFT,PHASES)
            ! Calculates orbital phases of the eclipses and quadratures.
            ! PHASES(1) and PHASES(3) are the prim and sec eclipse times
            ! PHASES(2) and PHASES(4) are the *photometric( quadratures.
      implicit none
      real*8 ECCIN,OMEGAIN          ! IN: eccentricity and peri.long.
      real*8 PSHIFT                 ! IN: time of primary eclipse
      real*8 PHASES(4)              ! OUT: phases of eclipses and quads
      real*8 PI,DEG2RAD             ! LOCAL: constants
      real*8 ECC,OMEGA              ! LOCAL: values of ecc and peri.long
      real*8 ECOSW,ESINW            ! LOCAL: values of combination terms
      real*8 EFAC,EFAC2             ! LOCAL: useful eccentricity factors
      real*8 TERM1,TERM2            ! LOCAL: calculation helper varibles
      real*8 PHASEDIFF              ! LOCAL: diff of prim and sec minima

      PI = atan(1.0d0) * 4.0d0
      DEG2RAD = 45.0d0 / atan(1.0d0)

            ! Get actual eccentricity and periastron longitude values
            ! from the input, which could be (e+10,w) or (ecosw,esinw)

      if ( ECCIN > 10.0d0 ) then
        ECC = ECCIN - 10.0d0
        OMEGA = OMEGAIN / DEG2RAD
        ECOSW = ECC * cos(OMEGA)
        ESINW = ECC * sin(OMEGA)
      else
        ECC = sqrt(ECCIN**2 + OMEGAIN**2)
        OMEGA = atan2(OMEGAIN,ECCIN)
        ECOSW = ECCIN
        ESINW = OMEGAIN
      end if

      EFAC = sqrt( (1.0d0-ECC) / (1.0d0+ECC) )
      EFAC2 = sqrt(1.0d0 - ECC**2)

!     write(6,*)" "
!     write(6,'(a16,2(f13.8))')"pi,deg2rad      ",pi,deg2rad
!     write(6,'(a16,2(f13.8))')"eccin,omegain   ",eccin,omegain
!     write(6,'(a16,2(f13.8))')"ecc,omega       ",ecc,omega
!     write(6,'(a16,2(f13.8))')"ecosw,esinw     ",ecosw,esinw
!     write(6,'(a16,2(f13.8))')"efac,efac2      ",efac,efac2
!     write(6,'(a16,2(f13.8))')"pshift,omegadeg ",pshift,omega*deg2rad

! The equation for phase difference comes from Hilditch (2001) page 238
! equation 5.66, originally credited to the monograph by  Kopal (1959).

      TERM1 = 2.0d0 * atan( ECOSW / EFAC2 )
      TERM2 = 2.0d0 * ECOSW * EFAC2 / (1.0d0 - ESINW**2)
      PHASEDIFF = ( PI + TERM1 + TERM2 ) / ( 2.0d0 * PI )

      PHASES(1) = PSHIFT
      PHASES(2) = PSHIFT + PHASEDIFF/2.0d0
      PHASES(3) = PSHIFT + PHASEDIFF
      PHASES(4) = PSHIFT + 0.50d0 + PHASEDIFF/2.0d0

!     write(6,*)" "
!     write(6,'(A19,2(F14.8))')"term1,term2        ",term1,term2
!     write(6,'(A19,2(F14.8))')"timediff,phasediff ",phasediff
!     write(6,'(A19,2(F14.8))')"eclipse phases     ",phases(1),phases(3)
!     write(6,'(A19,2(F14.8))')"quadrature phases  ",phases(2),phases(4)

      END SUBROUTINE ECQUADPHASES
!=======================================================================
      DOUBLEPRECISION FUNCTION GETMIN (TZERO,PERIOD,ECCIN,OMEGAIN,CICLE)
            ! Returns time of minimum for given cycle and ephemeris.  If
            ! the orbit's circular then the cycle number is used without
            ! restriction so can refer to any phase. If the orbit is ec-
            ! centric then the cycle number should be integer (indicates
            ! primary minimum) or half-integer (secondary minimum).
      implicit none
      real*8 TZERO,PERIOD           ! IN: reference time, orbital period
      real*8 ECCIN,OMEGAIN          ! IN: orbital (e,w) or (ecosw,esinw)
      real*8 CICLE                  ! IN: cycle number of minimum to use
      real*8 CICLEFRAC              ! LOCAL: fraction part of cycle nmbr
      real*8 ECC,OMEGA              ! LOCAL: eccentricity and peri.long.
      real*8 ECOSW,ESINW            ! LOCAL: eccentr'y combination terms
      real*8 PSEP,PHASES(4)         ! LOCAL: phase sep and useful phases
      real*8 PI,DEG2RAD             ! LOCAL: useful variables

      PI = atan(1.0d0) * 4.0d0
      DEG2RAD = 45.0d0 / atan(1.0d0)

            ! First must deal with the possibility that e and omega are
            ! actually e*cos(omega) and e*sin(omega)

      if ( ECCIN > 10.0d0 ) then
        ECC = ECCIN - 10.0d0
        OMEGA = OMEGAIN / DEG2RAD
        ECOSW = ECC * cos(OMEGA)
        ESINW = ECC * sin(OMEGA)
      else
        ECC = sqrt(ECCIN**2 + OMEGAIN**2)
        OMEGA = atan2(OMEGAIN,ECCIN)
        ECOSW = ECCIN
        ESINW = OMEGAIN
      end if

            ! If orbit is circular then simply use the orbital ephemeris
            ! If orbit is eccentric then call ECQUADPHASES to calculate
            ! the phase difference between primary and secondary minima.

      if ( abs(ECC) < 1.0d-7 ) then
        GETMIN = TZERO  +  PERIOD * CICLE
      else
        CICLEFRAC = mod(CICLE,1.0d0)

        if ( ( CICLEFRAC >= 0.0d0 .and. CICLEFRAC < 0.001d0 ) .or.
     &       ( CICLEFRAC > 0.999d0 .and. CICLEFRAC <= 1.0d0 ) ) then
          GETMIN = TZERO + PERIOD * CICLE
        else if ((CICLEFRAC > -0.501d0 .and. CICLEFRAC < -0.499d0) .or.
     &           (CICLEFRAC > 0.499d0 .and. CICLEFRAC < 0.501d0) ) then
          CALL ECQUADPHASES (ECCIN,OMEGAIN,0.0d0,PHASES)
          PSEP = PHASES(3) - PHASES(1)
          if ( PSEP < 0.0d0 ) PSEP = PSEP + 1.0d0
          GETMIN = TZERO + PERIOD * (CICLE-0.50d0+PSEP)
        else
          write(6,'(A37,A43)') "### ERROR: found a cycle number which",
     &                    " is not integer or half-integer. Abort.    "
          write(6,'(A18,F20.10)') "### Cycle number =",CICLE
          write(6,*) " "
          stop
        end if

      end if

      END FUNCTION GETMIN
!=======================================================================
!=======================================================================
      SUBROUTINE FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,
     &                    CHISQ,VERR,ERROR,NSINE,PSINE,NPOLY,PPOLY,NL3,
     &                    NECW,NESW,NUMINT,NINTERVAL)
            ! This subroutine calls the  MRQMIN  algorithm (Press et al,
            ! 1992, Numerical recipes in FORTRAN 77, p.678)  which finds
            ! the best-fitting EBOP model for the data using the
            ! Levenberg-Marquardt optimisation method.
            ! Unfortunately I have had to use a COMMON block here to
            ! avoid passing parameters through MRQMIN and related code.
      implicit none
      real*8 DATA(3,999999)               ! IN: Observational data
      integer DTYPE(999999)               ! IN: Types of datapoints
      integer NDATA,NLR,NMIN              ! IN: Numbers of datapoints
      real*8 V(67)                        ! IN: EBOP parameter values
      integer VARY(67)                    ! IN: Whether parameters vary
      integer LDTYPE(2)                   ! IN: Type of LD law for stars
      integer NSINE,NL3,NECW,NESW         ! IN: Numbers of sines and L3s
      integer PSINE(5)                    ! IN: Which par for each sine
      integer NPOLY,PPOLY(5)              ! IN: Similar for polynomials
      integer NUMINT                      ! IN: Number of numerical ints
      real*8 NINTERVAL                    ! IN: Time interval for numint
      integer ITER                        ! OUT: Number of iterations
      real*8 CHISQ                        ! OUT: Reduced chi-squared
      real*8 VERR(67)                     ! OUT: Formal parameter errors
      integer ERROR                       ! OUT: Whether fit successful
      real*8 X(999999),Y(999999),SIG(999999) ! LOCAL: data for MRQMIN
      real*8 LAMBDA                       ! LOCAL: Marquardt lambda
      real*8 COVAR(67,67),ALPHA(67,67)    ! LOCAL: Covariance etc matrix
      real*8 OCHISQ                       ! LOCAL: Previous chi-squared
      integer i,j,NVARY                   ! LOCAL: helpful integers
      integer CLDTYPE(2),CDTYPE(999999)   ! COMMON: for EBOP subroutine
      integer CNSINE,CPSINE(5)
      integer CNPOLY,CPPOLY(5)
      real*8 CNINTERVAL
      integer CNUMINT

      common / FOREBOP / CLDTYPE,CNSINE,CPSINE,CNPOLY,CPPOLY,CDTYPE,
     $                                                CNUMINT,CNINTERVAL
      ERROR = 0
!       CHISQ = 0.0d0

      CLDTYPE(1) = LDTYPE(1)
      CLDTYPE(2) = LDTYPE(2)
      do i = 1,NDATA
        CDTYPE(i) = DTYPE(i)
      end do
      CNSINE = NSINE
      CNPOLY = NPOLY
      do i = 1,5
        CPSINE(i) = PSINE(i)
        CPPOLY(i) = PPOLY(i)
      end do
      CNUMINT = NUMINT
      CNINTERVAL = NINTERVAL

      do i = 1,NDATA
        X(i)   = DATA(1,i)
        Y(i)   = DATA(2,i)
        SIG(i) = abs(DATA(3,i))
      end do

      NVARY = 0
      do i = 1,67
        if ( VARY(i) == 1 .or. VARY(i) == 3 ) NVARY = NVARY + 1
      end do

            ! Now find the best fit using MRQMIN. This requires an init-
            ! ial call with LAMBDA less than zero to initialise things.
            ! Then iterate to the best fit and assume this has been got
            ! once LAMBDA > 10^7. If no observational errors have been
            ! supplied then calculate them and send in once more. And
            ! finally set LAMBDA = 0.0 to get useful parameters out.

      LAMBDA = -1.0d0

      CALL MRQMIN (X,Y,SIG,NDATA,V,VARY,67,COVAR,ALPHA,67,CHISQ,
     &                                                     LAMBDA,ERROR)
      OCHISQ = 1.0d10

      do i = 1,50
        CALL MRQMIN (X,Y,SIG,NDATA,V,VARY,67,COVAR,ALPHA,67,CHISQ,
     &                                                     LAMBDA,ERROR)
        if ( LAMBDA >= 1.0d10 .and. abs(OCHISQ/CHISQ) < 1.001 ) exit
        if ( ERROR /= 0 ) exit
        OCHISQ = CHISQ
      end do
      ITER = i

      LAMBDA = 0.0d0
      CALL MRQMIN (X,Y,SIG,NDATA,V,VARY,67,COVAR,ALPHA,67,CHISQ,
     &                                                     LAMBDA,ERROR)

            ! Now record the formal errors outputted by MRQMIN

      do i = 1,67
        VERR(i) = sqrt(COVAR(i,i))
        if (DATA(3,1) < 0.0d0) VERR(i)=VERR(i)*sqrt(CHISQ/(NDATA-NVARY))
      end do

      if ( V(6) > 90.0d0 ) V(6) = 180.0d0 - V(6)

      if ( LDTYPE(2) == 0 ) then
        V(5) = V(4)
        VERR(5) = VERR(4)
        V(22) = V(21)
        VERR(22) = VERR(21)
      end if

      if ( DATA(3,1) > 0.0d0 )  CHISQ = CHISQ / (NDATA-NVARY)
      if ( DATA(3,1) <= 0.0d0 ) CHISQ = -1.0d0

      END SUBROUTINE FITEBOP
!=======================================================================
      SUBROUTINE EBOP (WHICHP,X,V,Y,DYDA,NCOEFFS,VARY,DODYDA)
            ! This evaluates the model value (Y) for one datapoint (X)
            ! INDEX is the datapoint number and DTYPE(INDEX) is its type
            ! Optionally (if DODYDA = 'y') it also calculates the numer-
            ! ical derivatives of X with respect to the variable params.
      implicit none
      integer WHICHP                ! IN: Which data point being studied
      real*8 X                      ! IN: Time to  calculate results for
      real*8 V(67)                  ! IN: The   photometric   parameters
      integer NCOEFFS               ! IN: Total number of adjustd params
      integer VARY(67)              ! IN:  Which params  being  adjusted
      character*1 DODYDA            ! IN: 'y' or 'n' todothe derivatives
      real*8 Y                      ! OUT:  Output result for input time
      real*8 DYDA(67)               ! OUT:   The numerical differentials
      real*8 LP,LS                  ! LOCAL: Light produced by each star
      real*8 OUT1,OUT2              ! LOCAL: Help in finding derivatives
      real*8 STORE                  ! LOCAL: Help in finding derivatives
      real*8 DV(67)                 ! LOCAL: Amount to perturb params by
      integer i,j,k,ERROR           ! LOCAL: Loop counters and errorflag
      real*8 GETMODEL               ! FUNCTION: Returns model evaluation
      integer LDTYPE(2)             ! IN/COMMON: LD law types  for stars
      integer DTYPE(999999)         ! IN/COMMON: number of LC datapoints
      integer NSINE,PSINE(5)
      integer NPOLY,PPOLY(5)
      real*8 R1,R2
      integer NUMINT                      ! IN: Number of numerical ints
      real*8 NINTERVAL                    ! IN: Time interval for numint

      common / FOREBOP / LDTYPE,NSINE,PSINE,NPOLY,PPOLY,DTYPE,
     &                                                  NUMINT,NINTERVAL

      CALL GET_DV(V,DV)

            ! First get the model prediction for this datapoint. And use
            ! this call to GETMODEL to get the reflection coefficients

      if ( V(2) >= 0.0d0 ) then
        R1 = V(2) / (1.0d0 + V(3))
        R2 = V(2) / (1.0d0 + (1.0d0/V(3)))
      else
        R1 = abs(V(2))
        R2 = V(3)
      end if

      Y = GETMODEL (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,X,DTYPE(WHICHP),
     &                                           LP,LS,NUMINT,NINTERVAL)
      if ( VARY(11) == -1 )  V(11) = 0.4d0 * (LS/(1.0d0-V(15))) * R1**2
      if ( VARY(12) == -1 )  V(12) = 0.4d0 * (LP/(1.0d0-V(15))) * R2**2

            ! Now for each adjustable parameter work out the adjustment
            ! to make to its value for calculating partial derivatives.
            ! NOTE: for the sine reference time the size ofthe numerical
            ! interval depends on the period of the sine: this must be
            ! treated separately to avoid numerical intervals which are
            ! either too large or too small.

      if ( DODYDA == 'y' ) then

        do i = 1,NCOEFFS
          if ( VARY(i) == 1 .or. VARY(i) == 3 ) then

            STORE = V(i)
            V(i) = STORE + DV(i)
            OUT1 = GETMODEL (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,X,
     &                             DTYPE(WHICHP),LP,LS,NUMINT,NINTERVAL)
            V(i) = STORE - DV(i)
            OUT2 = GETMODEL (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,X,
     &                             DTYPE(WHICHP),LP,LS,NUMINT,NINTERVAL)
            V(i) = STORE
            DYDA(i) = (OUT1 - OUT2) / (2.0d0 * DV(i))
          else
            DYDA(i) = 0.0d0
          end if
        end do

        if ( V(6) > 89.9d0) then
          DYDA(6) = DYDA(6) + (V(6)-89.9d0)**2
        end if

      else
        do i = 1,NCOEFFS
          DYDA(i) = 0.0d0
        end do
      end if


      END SUBROUTINE EBOP
!=======================================================================
!=======================================================================
!=======================================================================
      SUBROUTINE BIAX (R,Q,A,B,EPS)
            ! EBOP subroutine to calculate biaxial ellipsoid dimensions
            ! and oblateness for each star after Chandrasekhar (1933).
      implicit none
      real*8 R,Q,A,B,EPS

      if ( Q <= 0.0d0 )  then
        A = R
        B = R
        EPS = 0.0d0
      else
        A = R * ( 1.0d0 + (1.0d0 + 7.0d0*Q)/6.0d0 * R**3.0d0)
        B = R * ( 1.0d0 + (1.0d0 - 2.0d0*Q)/6.0d0 * R**3.0d0)
        EPS = (A - B) / A
        B = ( (1.0d0 - EPS) * R**3.0d0) ** (1.0d0/3.0d0)
        A = B / (1.0d0 - EPS)
      end if

      END SUBROUTINE BIAX
!=======================================================================
      SUBROUTINE GETEW (ECOSW,ESINW,E,W)
            ! EBOP subroutine to calculate e and w from e(cos)w e(sin)w
      implicit none
      real*8 ECOSW,ESINW,E,W

      if ( ECOSW == 0.0d0  .and.  ESINW == 0.0d0 ) then
        E = 0.0d0
        W = 0.0d0
      else
        W = atan2( ESINW,ECOSW )
        E = sqrt( ESINW*ESINW + ECOSW*ECOSW )
        W = W * 180.0d0 / 3.1415926536d0
      end if

      END SUBROUTINE GETEW
!=======================================================================
      SUBROUTINE LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,HJD,FMAG,LP,LS)
      implicit real*8 (a-h,o-z)
      real*8 V(67),HJD,GETPHASE
      real*8 LP,LS,LECL,LE
      real*8 LD1U,LD2U            ! linear LD coeff for each star
      real*8 LD1Q,LD2Q            ! quadratic LD coeff for each star
      real*8 LD1S,LD2S            ! square-root LD coeff for each star
      real*8 LD1L,LD2L            ! logarithmic LD coeff for each star
      real*8 LD1C,LD2C            ! cubic LD coeff for each star
      real*8 LDU,LDQ,LDS,LDL,LDC  ! LD coeffs for the star in question
      integer LDTYPE(2)           ! LD law type for both stars
      integer GIMENEZ             ! 1 to use original FMAX calculations
                                  ! 2 to use Gimenez' modified calcs
                                  ! 3 to use Gimenez' nonlinear LD calcs
      integer NSINE,PSINE(5)
      integer NPOLY,PPOLY(5)
      real*8 PHASE,SINT,SINP,SINA,SINTERM,LPMULT,LSMULT

!       data GIMENEZ / 3 /
!       data PI,TWOPI,RAD / 3.1415926536E0,6.28318531E0,0.0174532925E0 /
!       data LPMULT,LSMULT / 1.0 , 1.0 /

      GIMENEZ = 3
      LPMULT = 1.0d0
      LSMULT = 1.0d0
      PI = 3.1415926536d0
      TWOPI = 6.28318531d0
      RAD = 0.0174532925d0


C
C        DETERMINE PRIMARY AND SECONDARY BIAXIAL DIMENSIONS
C        USE SPHERICAL RADII FOR THE COMPUTATION OF ECLIPSE FUNCTIONS
C        USE OBLATENESSES FOR THE COMPUTATION OF THE OUTSIDE ECLIPSE
C        PHOTOMETRIC VARIATIONS WITH LIMB AND GRAVITY DARKENING
C

      if ( V(2) >= 0.0d0 ) then
        RP   = V(2) / ( 1.0d0 + V(3) )
        RS   = V(2) / ( 1.0d0 + (1.0d0/V(3)) )
      else
        RP   = abs( V(2) )
        RS   = V(3)
      end if

      if ( V(7) > 5.0d0 ) then
        ECOSW = (V(7)-10.0d0) * cos(V(8)/57.2957795d0)
        ESINW = (V(7)-10.0d0) * sin(V(8)/57.2957795d0)
      else
        ECOSW  = V( 7)
        ESINW  = V( 8)
      end if


      BS     = V(1)
      FI     = V(6)
      YP     = V(9)
      YS     = V(10)
      SP     = V(11)
      SS     = V(12)
      Q      = V(13)
      TANGL  = V(14)
      EL     = V(15)
      DPH    = 1.0d0 - V(16)
      SFACT  = V(17)
      DGAM   = V(18)

      LD1U = V(4)               ! linear terms
      LD2U = V(5)
      LD1L = 0.0d0              ! log terms
      LD2L = 0.0d0
      LD1S = 0.0d0              ! sqrt terms
      LD2S = 0.0d0
      LD1Q = 0.0d0              ! quadratic terms
      LD2Q = 0.0d0
      LD1C = 0.0d0              ! cubic terms
      LD2C = 0.0d0
      if ( LDTYPE(1) == 2 ) LD1L = V(21)
      if ( LDTYPE(1) == 3 ) LD1S = V(21)
      if ( LDTYPE(1) == 4 ) LD1Q = V(21)
      if ( LDTYPE(1) == 5 ) LD1C = V(21)
      if ( LDTYPE(2) == 2 ) LD2L = V(22)
      if ( LDTYPE(2) == 3 ) LD2S = V(22)
      if ( LDTYPE(2) == 4 ) LD2Q = V(22)
      if ( LDTYPE(2) == 5 ) LD2C = V(22)
      if ( LDTYPE(2) == 0 ) then
        LD2U = LD1U
        LD2L = LD1L
        LD2S = LD1S
        LD2Q = LD1Q
        LD2C = LD1C
      end if

      if ( NSINE > 0 ) then
        do i = 1,NSINE
          SINT = V(20+i*3)      ! sine reference time
          SINP = V(21+i*3)      ! sine period
          SINA = V(22+i*3)      ! sine amplitude
!           TWOPI =  atan(1.0d0) * 8.0d0
          SINTERM = SINA * sin( TWOPI * (HJD-SINT) / SINP )
          if ( PSINE(i) ==  1 )  BS = BS * (1.0d0+SINTERM)
          if ( PSINE(i) ==  2 )  RP = RP * (1.0d0+SINTERM)
          if ( PSINE(i) ==  3 )  RS = RS * (1.0d0+SINTERM)
          if ( PSINE(i) ==  6 )  FI = FI + SINTERM
          if ( PSINE(i) == 15 )  EL = EL * (1.0d0+SINTERM)
          if ( PSINE(i) == 17 )  SFACT = SFACT + SINTERM
          if ( PSINE(i) == -1 )  LPMULT = LPMULT * (1.0d0+SINTERM)
          if ( PSINE(i) == -2 )  LSMULT = LSMULT * (1.0d0+SINTERM)
        end do
      end if

      if ( NPOLY > 0 ) then
        do i = 1,NPOLY
          j = 32 + (i*6)
          SINTERM = V(j+1)*(HJD-V(j)) + V(j+2)*((HJD-V(j))**2) + V(j+3)*
     &  ((HJD-V(j))**3) + V(j+4)*((HJD-V(j))**4)+ V(j+5)*((HJD-V(j))**5)
          if ( PPOLY(i) ==  1 )  BS = BS + SINTERM   ! Yes this is POLY-
          if ( PPOLY(i) ==  2 )  RP = RP + SINTERM   ! TERM but's called
          if ( PPOLY(i) ==  3 )  RS = RS + SINTERM   ! SINETERM to avoid
          if ( PPOLY(i) ==  6 )  FI = FI + SINTERM    ! proliferation of
          if ( PPOLY(i) == 15 )  EL = EL + SINTERM           ! variables
          if ( PPOLY(i) == 17 )  SFACT = SFACT + SINTERM
          if ( PPOLY(i) == -1 )  LPMULT = LPMULT + SINTERM
          if ( PPOLY(i) == -2 )  LSMULT = LSMULT + SINTERM
        end do
      end if

      PHASE = GETPHASE(HJD,V(19),V(20))

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !         RS=RP*RATIO
      if ( Q <= 0.0d0 ) then
        CALL BIAX (RP,0.0d0,RPA,RPB,EP)
        CALL BIAX (RS,0.0d0,RSA,RSB,ES)
      else
        CALL BIAX (RP,Q,RPA,RPB,EP)
        CALL BIAX (RS,1.0d0/Q,RSA,RSB,ES)
      end if

C
C        CORRECT THE OBSERVED PHASE FOR ANY EPOCH ERROR IN EPHEMERIS
C
      THETA=PHASE+DPH
C
      SINI  = SIN(FI*RAD)
      SINI2 = SINI*SINI
      COSI2 = 1.0d0  - SINI2
C
C        TRANSLATE TIDAL LEAD/LAG ANGLE TO RADIANS
      TANGR=TANGL*RAD
C
C     EQUATION 9
C        CONVERT PHASE TO RADIANS
      FMN=THETA*TWOPI
C
C        GET CURRENT VALUES OF E, AND W
      CALL GETEW (ECOSW,ESINW,E,W)
C
C        TEST FOR CIRCULAR ORBIT
      IF (E)   17,20,17
   20 COSVW=COS(FMN)
      SINVW=SIN(FMN)
      RV=1.0d0
      GO TO 25
C
C        SOLUTION OF KEPLER'S EQUATION BY DIFFERENTIAL CORRECTIONS
C        (NON-ZERO ECCENTRICITY ONLY . . . )
C
C     EQUATION 6
C
   17 OMEGA = 450.0d0  - W
   23 IF (OMEGA - 360.0d0)         22,21,21
   21 OMEGA = OMEGA - 360.0d0
      GO TO 23
   22 OMEGA = OMEGA*RAD
C        SINE AND COSINE OF OMEGA
      COSW=COS(OMEGA)
      SINW=SIN(OMEGA)
C
C        COMPUTE MEAN ANOMALY CORRECTION TO PHASE
C        CORRESPONDING TO V=OMEGA=90-W
C        AT WHICH PHASE COS(V-OMEGA)=1
      E0=ATAN2(SQRT(1.0d0-E*E)*SINW,COSW+E)
C
C        MEAN ANOMALY OF MID-PRIMARY ECLIPSE
      FMA0=E0-E*SIN(E0)
C
C        MEAN ANOMALY
      FMA=FMN+FMA0
C     FIRST APPROXIMATION OF ECCENTRIC ANOMALY
      EA=FMA+E*SIN(FMA)
C
      DO 10 J=1,15
C        EVALUATE SINE AND COSINE OF ECCENTRIC ANOMALY
      SINE=SIN(EA)
      COSE=COS(EA)
      DENOM=1.0d0-E*COSE
      DISC=FMA-EA+E*SINE
      EA=EA+DISC/DENOM
C        TEST FOR CONVERGENCE
      IF (ABS(DISC) - 2.0d-5)     15,15,10
   10 CONTINUE
C
C
C        EVALUATE SINE AND COSINE OF TRUE ANOMALY
   15 COSV=(COSE-E)/DENOM
      SINV=SINE*SQRT(1.0d0-E*E)/DENOM
C
C        RADIUS VECTOR
      RV = (1.0d0-E*E)/(1.0d0+E*COSV)
C
C        THE PHOTOMETRIC PHASE ARGUMENT IN TERMS OF ORBIT PARAMETERS
C        VW = V-OMEGA
      COSVW=COSV*COSW+SINV*SINW
      SINVW=SINV*COSW-COSV*SINW
C
   25 COS2=COSVW*COSVW
      SIN2=1.0d0-COS2
C
      CSVWT=COS(TANGR)*COSVW-SIN(TANGR)*SINVW
C
C
C        PHOTOMETRIC EFFECTS
C
C

      FMAXP = 0.0d0
      FMAXS = 0.0d0
      DELTP = 0.0d0
      DELTS = 0.0d0
      SHORT = 0.0d0

!-----------------------------------------------------------------------
! Alvaro Gimenez and J Diaz-Cordoves have corrected the treatment of LD
! and stellar shapes.  This treatment can be used by putting GIMENEZ=2
! Their treatment for nonlinear LD can be used by putting GIMENEZ=3
!-----------------------------------------------------------------------
! This whole thing affects only the brightness normalisation of the two
! eclipsing stars: any problems here affect the radiative parameters
! but not the geometric parameters (radii, inclination etc).
!-----------------------------------------------------------------------

      if ( GIMENEZ==1 ) then                          ! LINEAR LD ONLY

!       FMAXP=((1.0E0-UP)+0.666666667E0*UP*(1.0E0+0.2E0*EP)) ! Original
!      1      *(1.0E0+3.0E0*YP*EP)/(1.0E0-EP)                ! lines
!       FMAXS=((1.0E0-US)+0.666666667E0*US*(1.0E0+0.2E0*ES)) ! if the
!      1      *(1.0E0+3.0E0*YS*ES)/(1.0E0-ES)                ! stars
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP    ! are
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES    ! oblate
!       SHORT=SINI2*CSVWT*CSVWT

!    26 FMAXP=1.0E0-UP/3.0E0                                 ! Original
!       FMAXS=1.0E0-US/3.0E0                                 ! lines if
!       DELTP=0.0E0                                          ! the stars
!       DELTS=0.0E0                                          ! are
!       SHORT=0.0                                            ! spherical

        if ( Q >= 0.0d0 ) then
          FMAXP=((1.0d0-LD1U)+0.666666667d0*LD1U*(1.0d0+0.2d0*EP))
     1        *(1.0d0+3.0d0*YP*EP)/(1.0d0-EP)
          FMAXS=((1.0d0-LD2U)+0.666666667d0*LD2U*(1.0d0+0.2d0*ES))
     1        *(1.0d0+3.0d0*YS*ES)/(1.0d0-ES)
          DELTP=(15.0d0+LD1U)/(15.0d0-5.0d0*LD1U)*(1.0d0+YP)*EP
          DELTS=(15.0d0+LD2U)/(15.0d0-5.0d0*LD2U)*(1.0d0+YS)*ES
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0d0-LD1U/3.0d0
          FMAXS=1.0d0-LD2U/3.0d0
          DELTP=0.0d0
          DELTS=0.0d0
          SHORT=0.0
        end if
!-----------------------------------------------------------------------

      else if ( GIMENEZ==2 ) then                     ! LINEAR LD ONLY

!       FMAXP=(1.0E0-UP*(1.0E0-2.0E0/5.0E0*EP)/3.0E0+YP*EP   ! Original
!      1      *(3.0E0-13.0E0/15.0E0*UP))/(1.0E0-EP)          ! lines
!       FMAXS=(1.0E0-US*(1.0E0-2.0E0/5.0E0*ES)/3.0E0+YS*ES   ! if the
!      1      *(3.0E0-13.0E0/15.0E0*US))/(1.0E0-ES)          ! stars
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP    ! are
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES    ! oblate
!       SHORT=SINI2*CSVWT*CSVWT

!    26 FMAXP=1.0E0-UP/3.0E0                                 ! Original
!       FMAXS=1.0E0-US/3.0E0                                 ! lines if
!       DELTP=0.0E0                                          ! the stars
!       DELTS=0.0E0                                          ! are
!       SHORT=0.0                                            ! spherical

        if ( Q >= 0.0d0 ) then
          FMAXP=(1.0d0-LD1U*(1.0d0-2.0d0/5.0d0*EP)/3.0d0+YP*EP
     1          *(3.0d0-13.0d0/15.0d0*LD1U))/(1.0d0-EP)
          FMAXS=(1.0d0-LD2U*(1.0d0-2.0d0/5.0d0*ES)/3.0d0+YS*ES
     1          *(3.0d0-13.0d0/15.0d0*LD2U))/(1.0d0-ES)
          DELTP=(15.0d0+LD1U)/(15.0d0-5.0d0*LD1U)*(1.0d0+YP)*EP
          DELTS=(15.0d0+LD2U)/(15.0d0-5.0d0*LD2U)*(1.0d0+YS)*ES
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0d0-LD1U/3.0d0
          FMAXS=1.0d0-LD2U/3.0d0
          DELTP=0.0d0
          DELTS=0.0d0
          SHORT=0.0d0
        end if
!-----------------------------------------------------------------------
! And this is Gimenez's code for including nonlinear LD. He includes
! the linear (UP), quadratic (UP, U2P) and square-root (UP, U3P) laws.
!-----------------------------------------------------------------------

      else if ( GIMENEZ==3 ) then

!      FMAXP=1.0E0-UP*(1.0E0-2.0E0*EP/5.0E0)/3.0E0-
!     1      U2P*(1.0E0-3.0E0*EP/5.0E0)/6.0E0-
!     1      U3P*(1.0E0-4.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP
!     1      *(1.5E0-13.0E0*UP/30.0E0-U2P/5.0E0-23.0E0*U3P/90.0E0)
!      FMAXP=FMAXP/(1.0E0-EP)
!      FMINP=1.0E0-UP*(1.0E0+4.0E0*EP/5.0E0)/3.0E0-
!     1      U2P*(1.0E0+6.0E0*EP/5.0E0)/6.0E0-
!     1      U3P*(1.0E0+8.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP
!     1      *(1.0E0-7.0E0*UP/15.0E0-4.0E0*U2P/15.0E0-13.0E0*U3P/45.0E0)
!      FMINS=1.0E0-US*(1.0E0+4.0E0*ES/5.0E0)/3.0E0-
!     1      U2S*(1.0E0+6.0E0*ES/5.0E0)/6.0E0-
!     1      U3S*(1.0E0+8.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES
!     1      *(1.0E0-7.0E0*US/15.0E0-4.0E0*U2S/15.0E0-13.0E0*U3S/45.0E0)
!      FMAXS=1.0E0-US*(1.0E0-2.0E0*ES/5.0E0)/3.0E0-
!     1      U2S*(1.0E0-3.0E0*ES/5.0E0)/6.0E0-
!     1      U3S*(1.0E0-4.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES
!     1      *(1.5E0-13.0E0*US/30.0E0-U2S/5.0E0-23.0E0*U3S/90.0E0)
!      FMAXS=FMAXS/(1.0E0-ES)
!      DELTP=1.0E0-FMINP/FMAXP
!      DELTS=1.0E0-FMINS/FMAXS
!      SHORT=SINI2*CSVWT*CSVWT

!   26 FMAXP=1.0E0-UP/3.0E0-U2P/6.0E0-U3P/5.0E0
!      FMAXS=1.0E0-US/3.0E0-U2S/6.0E0-U3S/5.0E0
!      DELTP=0.0E0
!      DELTS=0.0E0
!      SHORT=0.0

        if ( Q >= 0.0d0 .or. LDTYPE(1)==1 .or. LDTYPE(1)==5 .or.
     &                         LDTYPE(2)==1 .or. LDTYPE(2)==5 ) then
          FMAXP=1.0d0-LD1U*(1.0d0-2.0d0*EP/5.0d0)/3.0d0-
     &          LD1Q*(1.0d0-3.0d0*EP/5.0d0)/6.0d0-
     &          LD1S*(1.0d0-4.0d0*EP/9.0d0)/5.0d0+2.0d0*YP*EP
     &         *(1.5E0-13.0d0*LD1U/30.0d0-LD1Q/5.0d0-23.0d0*LD1S/90.0d0)
          FMAXP=FMAXP/(1.0d0-EP)
          FMINP=1.0d0-LD1U*(1.0d0+4.0d0*EP/5.0d0)/3.0d0-
     &          LD1Q*(1.0d0+6.0d0*EP/5.0d0)/6.0d0-
     &          LD1S*(1.0d0+8.0d0*EP/9.0d0)/5.0d0+2.0d0*YP*EP
     &   *(1.0d0-7.0d0*LD1U/15.0d0-4.0d0*LD1Q/15.0d0-13.0d0*LD1S/45.0d0)
          FMINS=1.0d0-LD2U*(1.0d0+4.0d0*ES/5.0d0)/3.0d0-
     &          LD2Q*(1.0d0+6.0d0*ES/5.0d0)/6.0d0-
     &          LD2S*(1.0d0+8.0d0*ES/9.0d0)/5.0d0+2.0d0*YS*ES
     &   *(1.0d0-7.0d0*LD2U/15.0d0-4.0d0*LD2Q/15.0d0-13.0d0*LD2S/45.0d0)
          FMAXS=1.0d0-LD2U*(1.0d0-2.0d0*ES/5.0d0)/3.0d0-
     &          LD2Q*(1.0d0-3.0d0*ES/5.0d0)/6.0d0-
     &          LD2S*(1.0d0-4.0d0*ES/9.0d0)/5.0d0+2.0d0*YS*ES
     &         *(1.5E0-13.0d0*LD2U/30.0d0-LD2Q/5.0d0-23.0d0*LD2S/90.0d0)
          FMAXS=FMAXS/(1.0d0-ES)
          DELTP=1.0d0-FMINP/FMAXP
          DELTS=1.0d0-FMINS/FMAXS
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0-LD1U/3.0-LD1Q/6.0-LD1S/5.0+LD1L*2.0/9.0-LD1C/10.0
          FMAXS=1.0-LD2U/3.0-LD2Q/6.0-LD2S/5.0+LD2L*2.0/9.0-LD2C/10.0
          DELTP=0.0d0
          DELTS=0.0d0
          SHORT=0.0d0
        end if
!----------------------------------------------------------------------
      end if
!----------------------------------------------------------------------
! Complete original code before the above messing:
! C
! C
! C        PHOTOMETRIC EFFECTS
! C
! C
! C        TEST FOR SIMPLE CASE OF TWO SPHERICAL STARS
!       IF (EP .EQ. 0.  .AND.  ES .EQ. 0.)   GO TO 26
! C
! C        EITHER OR BOTH STARS ARE OBLATE
! C
!       FMAXP=((1.0E0-UP)+0.666666667E0*UP*(1.0E0+0.2E0*EP))
!      1      *(1.0E0+3.0E0*YP*EP)/(1.0E0-EP)
!       FMAXS=((1.0E0-US)+0.666666667E0*US*(1.0E0+0.2E0*ES))
!      1      *(1.0E0+3.0E0*YS*ES)/(1.0E0-ES)
! C        CHANGE IN INTENSITY RATIO DUE TO OBLATENESS RELATED VARIABLES
! C        FROM QUADRATURE TO MINIMUM
! C        FACE ON TO END ON
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES
! C        FORE-SHORTENING FUNCTION OF OBLATENESS
!       SHORT=SINI2*CSVWT*CSVWT
!       GO TO 27
! C
! C        BOTH STARS ARE SPHERICAL
! C
!    26 FMAXP=1.0E0-UP/3.0E0
!       FMAXS=1.0E0-US/3.0E0
!       DELTP=0.0E0
!       DELTS=0.0E0
!       SHORT=0.0
!----------------------------------------------------------------------

C
C        UN-NORMALIZED BRIGHTNESS OF STELLAR COMPONENTS AT QUADRATURE
   27 OP=PI*RPB*RPB*FMAXP
      OS=PI*RSB*RSB*FMAXS*BS
C        THE NORMALIZING FACTOR
      OTOT=OP+OS
C        BRIGHTNESS CONTRIBUTION FROM EACH COMPONENT
      LP=OP/OTOT*(1.0d0-DELTP*SHORT)
      LS=OS/OTOT*(1.0d0-DELTS*SHORT)
C
C        REFLECTION AND RERADIATION EQUATION
      IF (SP .EQ. 0.0d0  .AND.  SS .EQ. 0.0d0)   GO TO 28
      HEAT=SINI*COSVW
      HEAT2=0.5d0+0.5d0*HEAT*HEAT
      DLP=SP*(HEAT2+HEAT)
      DLS=SS*(HEAT2-HEAT)
      GO TO 29
   28 DLP=0.0d0
      DLS=0.0d0
C
C        WHICH ECLIPSE COULD THIS BE
   29 IF (COSVW)         40,40,30
C
C     PRIMARY ECLIPSE
C
   30 R1 = RP
      R2 = RS
!----------------------------------------------------------------------!
! JKT mod (10/8/2006): the line these replaced was      UU = UP        !
!----------------------------------------------------------------------!
      LDU = LD1U                                                       !
      LDL = LD1L                                                       !
      LDS = LD1S                                                       !
      LDQ = LD1Q                                                       !
      LDC = LD1C                                                       !
!----------------------------------------------------------------------!
      LE=LP
      DLE=DLP
      GO TO 60
C
C
C     SECONDARY ECLIPSE
C
   40 R1 = RS
      R2 = RP
!-----------------------------------------------------------------------
! JKT mod (10/8/2006): the line these replaced was      UU = US        !
!----------------------------------------------------------------------!
      LDU = LD2U                                                       !
      LDL = LD2L                                                       !
      LDS = LD2S                                                       !
      LDQ = LD2Q                                                       !
      LDC = LD2C                                                       !
!----------------------------------------------------------------------!
      LE=LS
      DLE=DLS
C
   60 SUM = 0.0d0
      ALAST = 0.0d0
      AREA=0.0d0
C
C     EQUATION  5
C
      DD = SINVW*SINVW + COSVW*COSVW*COSI2
      IF (DD .LE. 1.0d-6)  DD=0.0d0
      DD = DD*RV*RV
      D = SQRT(ABS(DD))
      R22 = R2*R2
C
C     EQUATION 17
C
      GAMN = 90.01d0*RAD
      DGAMA = DGAM*RAD
      DGM = DGAMA/2.0d0
      RK = 0.0d0
      GAM = 0.0d0
   50 GAM = GAM + DGAMA
C        HAS LIMIT OF INTEGRATION BEEN REACHED
      IF (GAM - GAMN)              48,48,49
C
   48 RR = R1*SIN(GAM)
      R12 = RR*RR
C
      AA = 0.0d0
C        ARE THE PROJECTED DISKS CONCENTRIC
      IF (D)                       405,406,405
  406 IF (RR - R2)                 230,230,403
  403 IF (RK - R2)                 404, 49, 49
  404 AA = PI*R22
      GO TO 215
C        TEST FOR NO ECLIPSE
  405 IF (D-R1-R2)                 240,216,216
  216 SUM = 0.0d0
      GO TO 49
C        DECIDE WHICH AREA EQUATIONS FOR NON-CONCENTRIC ECLIPSE
  240 IF (D-RR-R2)                 245,215,215
  245 IF (D-R2+RR)                 230,230,250
  250 IF (R1-R2)                   255,255,280
  255 IF (DD-R22+R12)              205,210,210
  280 IF (D-RR+R2)                 290,260,260
  260 IF (RR-R2)                   255,255,265
  265 IF (DD-R12+R22)              270,210,210
C
C     EQUATION 12
C
  270 S1 = ABS((R12 - R22 - DD)*0.5d0/D)
      A1 = ABS(R2-S1)
      B2 = ABS(RR-S1-D  )
      AA=PI*R22-(R22*ACOS((R2-A1)/R2)
     1   - (R2-A1)*SQRT(2.0d0*R2*A1-A1*A1))
     2   +R12*ACOS((RR-B2)/RR)-(RR-B2)*SQRT(2.0d0*RR*B2-B2*B2)
      GO TO 215
C
  290 IF (R1 - R2 - D)             260,260,295
  295 IF (RK - R2 - D)             300,215,215
  300 RR = R2 + D
      R12 = RR*RR
      GAMN = 0.0d0
      GO TO 260
C
  230 AA = PI*R12
      GO TO 215
C
C     EQUATION 10
C
  205 S = ABS((R12 - R22 + DD)*0.5d0/D)
      A = ABS(RR-S)
      B1 = ABS(R2-S-D)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0d0*RR*A - A*A)
      AB1 = R22*ACOS((R2-B1)/R2) - (R2-B1)*SQRT(2.0d0*R2*B1-B1*B1)
      AA = PI*R12 - A1 + AB1
      GO TO 215
C
C     EQUATION 1
C
  210 S = ABS((R12 - R22 + DD)*0.5d0/D)
      A = ABS(RR-S)
      B = ABS(S-D+R2)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0d0*RR*A - A*A)
      AA1 = R22*ACOS((R2-B)/R2) - (R2-B)*SQRT(2.0d0*R2*B - B*B)
      AA = A1 + AA1
C
  215 DAREA = AA - ALAST
!----------------------------------------------------------------------!
! JKT modification (10/9/2006). The removed line was:                  !
!     SUM = SUM + DAREA*(1.0d0  - UU + UU*COS(GAM-DGM))                !
!----------------------------------------------------------------------!
      COSGAM = cos(GAM-DGM)                                            !
      SUM = SUM + DAREA*(1.0d0 - LDU*(1.0d0-COSGAM)                    !
     &          - LDL*COSGAM*log(COSGAM) - LDS*(1.0d0-sqrt(COSGAM))    !
     &         - LDQ*(1.0d0-COSGAM)**2 - LDC*(1.0d0-COSGAM)**3)        !
!----------------------------------------------------------------------!
      ALAST = AA
      AREA = AREA + DAREA
C
      RK = RR
      GO TO 50
C
C        LIGHT LOSS FROM ECLIPSE
C
   49 ADISK = PI*R1*R1
!----------------------------------------------------------------------!
! JKT modification (10/9/2006).  See 1992A+A...259..227D for more info.!
! The removed line was:           ALPHA = SUM/(ADISK*(1.0-UU/3.0))     !
!----------------------------------------------------------------------!
      ALPHA = 1.0d0 - LDU/3.0d0 + LDL*2.0d0/9.0d0 -                    !
     &          LDS/5.0d0 - LDQ/6.0d0 - LDC/10.0d0                     !
      ALPHA = SUM/(ADISK*ALPHA)                                        !
!----------------------------------------------------------------------!
      LECL = ALPHA*LE
      AREA = AREA/ADISK
      REFL=DLP+DLS-AREA*DLE
C
C        THEORETICAL INTENSITY WITH THIRD LIGHT AND QUADRATURE
C        SCALE FACTOR APPLIED
C
!----------------------------------------------------------------------!
! This is the original line from EBOP:
!----------------------------------------------------------------------!
!      FLITE = ((LP+LS-LECL+REFL)*(1.0d0-EL)+EL)*SFACT
!----------------------------------------------------------------------!

      LP = LP * LPMULT               ! sine/poly applied to star 1 light
      LS = LS * LSMULT               ! sine/poly applied to star 2 light
      FLITE = ((LP+LS-LECL+REFL)*(1.0d0-EL)+EL)
      FMAG = -2.5d0 * log10(FLITE) + SFACT

      LP = LP * (1.0d0-EL)           ! account for third light *AFTER*
      LS = LS * (1.0d0-EL)           ! FLITE and FMAG have been found

      END
!=======================================================================
!=======================================================================
!=================     NUMERICAL RECIPES SUBROUTINES     ===============
!=======================================================================
!=======================================================================
      SUBROUTINE MRQMIN (x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     &                                                     alamda,ifail)
      implicit none
      integer NDATA,MA,NCA          ! IN: NDATA, numcoeffs, maxcoeffs
      real*8 X(ndata),Y(ndata)      ! IN: data to be fitted
      real*8 SIG(ndata)             ! IN: data errors in y
      real*8 A(ma)                  ! IN: coefficients
      integer IA(ma)                ! IN: adjust (1) or fix (0) coeffs
      real*8 COVAR(nca,nca)         ! OUT: curvature matrix
      real*8 ALPHA(nca,nca)         ! OUT: covariance matrix
      real*8 ALAMDA                 ! IN/OUT: Marquardt lambda factor
      real*8 CHISQ                  ! OUT: chi-squared of the fit

      integer MMAX,j,k,l,m,mfit,ifail
      parameter (MMAX = 67)
      real*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit

!       write(6,*)"a"
!       write(6,*)a
!       write(6,*)"ia"
!       write(6,*)ia
!       write(6,*)ndata,x(1),x(ndata)
!       write(6,*)ndata,y(1),y(ndata)
!       write(6,*)ndata,sig(1),sig(ndata)
!       write(6,*)"nca,chisq,alamda,ifail"
!       write(6,*)nca,chisq,alamda,ifail

      if(alamda.lt.0.0d0)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).eq.1) mfit=mfit+1
11      continue
        alamda=0.0010d0
        OCHISQ = 1.0d10
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,OCHISQ,chisq)
        OCHISQ = CHISQ
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif

      j=0
      do 14 l=1,ma
        if(ia(l).eq.1) then
          j=j+1
          k=0
          do 13 m=1,ma
            if(ia(m).eq.1) then
              k=k+1
              covar(j,k)=alpha(j,k)
            endif
13        continue
          covar(j,j)=alpha(j,j)*(1.0d0+alamda)
          da(j)=beta(j)
        endif
14    continue
      call gaussj(covar,mfit,nca,da,1,1,ifail)
      if ( ifail /= 0 ) return

      if(alamda.eq.0.0d0)then
        call covsrt(covar,nca,ma,ia,mfit)
        return
      endif

      j=0
      do 15 l=1,ma
        if(ia(l).eq.1) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,OCHISQ,chisq)
      if(chisq.lt.ochisq)then
        alamda=0.1d0*alamda
        ochisq=chisq
        j=0
        do 17 l=1,ma
          if(ia(l).eq.1) then
            j=j+1
            k=0
            do 16 m=1,ma
              if(ia(m).eq.1) then
                k=k+1
                alpha(j,k)=covar(j,k)
              endif
16          continue
            beta(j)=da(j)
            a(l)=atry(l)
          endif
17      continue
      else
        alamda=10.0d0*alamda
        chisq=ochisq
      endif
      return
      END
!-----------------------------------------------------------------------
      SUBROUTINE mrqcof (x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,OCHISQ,
     &                                                            chisq)
      implicit none
      integer ma,nalp,NDATA,ia(ma),MMAX,mfit,i,j,k,l,m
      parameter ( MMAX = 67 )
      real*8 OCHISQ,CHISQ
      real*8 a(ma),alpha(nalp,nalp),beta(ma)
      real*8 wt,dyda(MMAX)
      real*8 X(NDATA),Y(NDATA),SIG(NDATA)
      real*8 YMOD(NDATA),DY(NDATA),SIG2(NDATA)

            ! Modified by JKT, 16/05/2007
            ! Now if CHISQ > OCHISQ it doesn't waste computing time in
            ! calculating the derivatives (as these are not used)

      mfit=0
      do j = 1,ma
        if ( ia(j)==1 ) mfit = mfit+1
      end do

      do j = 1,mfit
        do k = 1,j
          alpha(j,k) = 0.0d0
        end do
        beta(j) = 0.0d0
      end do

      CHISQ = 0.0d0
      do i = 1,NDATA

        CALL EBOP (i,x(i),a,ymod(i),dyda,ma,ia,'n')
!       CALL EBOP (INDEX,X,V,Y,DYDA,NCOEFFS,VARY,DODYDA)

        SIG2(i) = 1.0d0 / (SIG(i)*SIG(i))
        DY(i) = Y(i) - YMOD(i)
        CHISQ = CHISQ + DY(i)*DY(i)*SIG2(i)
      end do
      if ( CHISQ > OCHISQ ) return

      do i = 1,NDATA
        CALL EBOP (i,x(i),a,ymod(i),dyda,ma,ia,'y')

        j=0
        do l=1,ma
          if(ia(l).eq.1) then
            j=j+1
            wt=dyda(l)*sig2(i)
            k=0
            do m=1,l
              if(ia(m).eq.1) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(m)
              endif
            end do
            beta(j)=beta(j)+dy(i)*wt
          endif
        end do
      end do

      do j = 2,mfit
        do k = 1,j-1
          alpha(k,j) = alpha(j,k)
        end do
      end do

      END
!-----------------------------------------------------------------------
      SUBROUTINE GAUSSJ (a,n,np,b,m,mp,ifail)
      implicit none
      integer m,mp,n,np,NMAX,ifail
      real*8 a(np,np),b(np,mp),big,dum,pivinv
      parameter ( NMAX = 67 )
      integer i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      ifail=0
      irow=0
      icol=0
!       write(6,*)a
!       write(6,*)" "
      do j=1,n
        ipiv(j)=0
      end do
      do i=1,n
        big=0.0d0
        do j=1,n
          if(ipiv(j).ne.1)then
            do k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(6,'(A25,A27,I3,A1)') "### PROBLEM A: singular m",
     &                               "atrix in gaussj (parameter:",k,")"
                ifail = 1
                return
              endif
            end do
          endif
        end do
!         print'(i5,$)',icol
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
          end do
          do l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
          end do
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.0d0) then
          write(6,'(A31,A21,I3,A1)') "### PROBLEM B: singular matrix ",
     &                                  "in gaussj (parameter:",ICOL,")"
          ifail = 1
          return
        end if
        pivinv=1.0d0/a(icol,icol)
        a(icol,icol)=1.0d0
        do l=1,n
          a(icol,l)=a(icol,l)*pivinv
        end do
        do l=1,m
          b(icol,l)=b(icol,l)*pivinv
        end do
        do ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.0d0
            do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
            end do
            do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
            end do
          endif
        end do
      end do
      do l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
          end do
        endif
      end do
      return
      END
!-----------------------------------------------------------------------
      SUBROUTINE COVSRT (covar,npc,ma,ia,mfit)
      implicit none
      integer ma,mfit,npc,ia(ma),i,j,k
      real*8 covar(npc,npc),swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.0d0
          covar(j,i)=0.0d0
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).eq.1)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END
!=======================================================================
!=======================================================================
!======================================================================
      INTEGER FUNCTION SEEDSTART ()
            ! This uses the Fortran-intrinsic function SYSTEM_CLOCK to
            ! generate an integer SEED dependent on real time.
            ! SYSTEM_CLOCK returns the present time, the number of
            ! times a second this changes (100 on my computer) and the
            ! maximum value of present time (after which present time
            ! starts again from zero) (2147483647 on my computer).
            ! SEED is outputted as a four-figure integer.
      implicit none
      integer a,b                   ! unused constant values

      CALL SYSTEM_CLOCK (SEEDSTART,a,b)
                        ! SEEDSTART becomes the current clock count

      END FUNCTION SEEDSTART
!=====================================================================
      DOUBLEPRECISION FUNCTION RANDOM (SEED)
            ! Minimal random number generator from Park and Miller.
            ! Return a random deviate between 0.0 and 1.0 with a flat
            ! ditribution. Its period is (2**31 - 1) = 2.15e9
            ! The constants are the best found by Park and Miller.
            ! Set and reset SEED to any integer value to inititalize
            ! the sequence and do not change it thereafter - the code
            ! updates SEED itself.
            ! Does not work if SEED=0, unless it bit-swaps.

      implicit none
      integer SEED                  ! All-important seed value
      integer IA,IM,IR,IQ,MASK      ! Various constants
      real*8 AM                     ! Another constant
      integer K                     ! Yet another constant

      IA = 16807
      IM = 2147483647               ! 2**31 ie largest integer the ..
      AM = 1.0d0 / IM               ! .. computer can handle.
      IQ = 127773
      IR = 2836
      MASK = 123459876

      SEED = IEOR(SEED,MASK)              ! IEOR is integer exclusive-or
      K = SEED / IQ                       !  bit-swapping, the non-nume-
      SEED = IA * (SEED - K*IQ) - K*IR    !  ric operation needed by all
      if (SEED < 0.0d0) SEED = SEED + IM  !  random number generstors.
      RANDOM = AM * SEED
      SEED = IEOR(SEED,MASK)              ! Seed is updated here

      END FUNCTION RANDOM
!=======================================================================
      DOUBLEPRECISION FUNCTION RANDOMG (SEED,MEAN,SD)
            ! This produces a random number with a Gaussian distribution
            ! It uses RANDOM to generate a random number distribution
            ! with zero mean and unit variance then scales the result
            ! with SD and offsets it with MEAN to produce the effect.
      implicit none
      integer SEED                  ! Seeding integer
      real*8 MEAN,SD                ! Desired mean and S.D. (variance)
      integer ISET                  ! LOCAL: integer variable
      real*8 FAC,GSET,RSQ,V1,V2     ! LOCAL: real variables
      real*8 RANDOM            ! FUNCTION: flat-distrib random generator
      SAVE ISET,GSET
      integer i

      if ( SEED < 0 ) ISET = 0      ! Reinitialise
      if ( ISET == 0 ) then
        do i = 1,100000
          V1 = 2.0d0 * RANDOM(SEED) - 1.0d0
          V2 = 2.0d0 * RANDOM(SEED) - 1.0d0
          RSQ = V1**2 + V2**2
          if ( RSQ /= 0.0d0 .and. RSQ <= 1.0d0 ) exit
        end do
        FAC = SQRT( -2.0d0 * log(RSQ) / RSQ )
        GSET = V1 * FAC
        RANDOMG = V2 * FAC
        ISET = 1
      else
        RANDOMG = GSET
        ISET = 0
      end if
      RANDOMG = ( RANDOMG * SD ) + MEAN

      END FUNCTION RANDOMG
!=====================================================================
      DOUBLEPRECISION FUNCTION SELLECT (ARRAY,NUM,K)
            ! Returns the Kth smallest value in ARRAY(1:NUM).
            ! ARRAY is simply sorted during this procedure.
      implicit none
      integer NUM                   ! IN: size of the array
      real*8 ARRAY(NUM)             ! IN: Input array
      integer K                     ! OUT: value to find
      real*8 STORE                  ! LOCAL: storage variable
      integer i,j,TAG               ! LOCAL: loop counters and a flag

      TAG = 0
      do i = 1,NUM
        STORE = 10.0**10.0
        do j = i,NUM
            if ( ARRAY(j) < STORE )  then
            STORE = ARRAY(j)
            TAG = j
          end if
        end do
        ARRAY(TAG) = ARRAY(i)
        ARRAY(i) = STORE
      end do
      SELLECT = ARRAY(K)

      END FUNCTION SELLECT
!=======================================================================
      DOUBLEPRECISION FUNCTION SIGMA (ARRAY,NUM)
            ! Returns the standard deviation of array ARRAY
      implicit none
      integer NUM                   ! IN: size of the array
      real*8 ARRAY(NUM)             ! IN: array
      real*8 MEAN,VAR,SUM,SUMSQ     ! LOCAL: properties of ARRAY values
      integer i                     ! LOCAL: loop counter

      SIGMA = 0.0d0
      SUM = 0.0d0
      SUMSQ = 0.0d0
      do i = 1,NUM
        SUM = SUM + ARRAY(i)
        SUMSQ = SUMSQ + ARRAY(i)**2
      end do

      MEAN = SUM / NUM
      VAR = (SUMSQ / NUM) - MEAN**2
      if ( VAR > 0.0d0 )  SIGMA = sqrt( VAR )

      END FUNCTION SIGMA
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
