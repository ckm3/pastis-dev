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