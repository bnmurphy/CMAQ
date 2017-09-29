
!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!


C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/yoj/arc/CCTM/src/cloud/cloud_acm_ae5/scavwdep.F,v 1.7 2011/10/21 16:10:28 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SCAVWDEP ( JDATE, JTIME, WTBAR, WCBAR, TBARC, PBARC,
     &                      CTHK1, AIRM, PRATE1, TAUCLD, POLC, CEND,
     &                      REMOV, REMOVAC, ALFA0NUC, ALFA2NUC, ALFA3NUC,
     &                      ALFA0AIT, ALFA2AIT, ALFA3AIT )
C-----------------------------------------------------------------------
C  Description: Compute simplistic incloud scavenging and wet removal
 
C  Revision History:
C      No   Date   Who  What
C      -- -------- ---  -----------------------------------------
C       0 01/15/98 sjr  created program
C       1 07/27/99 sjr  added scavenging coefficient scaling factors
C       3 Dec 00   Jeff move CGRID_MAP into f90 module
C       4 Dec 02   sjr  revised calls to HLCONST
C       5 Jun 10   J.Young: convert for Namelist redesign
C       6 Mar 11   S.Roselle: replaced I/O API include files with UTILIO_DEFN;
C                       removed deprecated TRIMLEN
C       6 Mar 11   B.Hutzell: added TRACER string for AE_SPCs to omit while
C                       computing Aitken mode mass
C       7 Jul 15   B.Murphy: added scavenging calculation for nucleation
C                       mode aerosols in addition to just Aitken
 
C  Called by:  RADMCLD and RESCLD
 
C  Calls the following subroutines:  GETALPHA
 
C  Calls the following functions:  HLCONST
 
C  Arguments    Type      I/O       Description
C  ---------   -------  ------  --------------------------------
C    JDATE     integer   input  current model julian date (yyyyddd)
C    JTIME     integer   input  current model time (hhmmss)
C    WTBAR      real     input  avg total water content (kg/m3)
C    WCBAR      real     input  avg liquid water content (kg/m3)
C    TBARC      real     input  avg cloud temperature (K)
C    PBARC      real     input  avg cloud pressure (Pa)
C    CTHK1      real     input  cloud thickness (m)
C    AIRM       real     input  total air mass (moles/m2) in cloudy air
C    PRATE1     real     input  precipitation rate (mm/hr)
C    TAUCLD     real     input  cloud lifetime (s)
C    POLC       real     input  ave vert conc incloud (moles sp/ mole air)
C    CEND       real    output  ending incloud conc (moles/mole)
C    REMOV      real    output  moles/m2 or mm*mol/lit scavenged
C    REMOVAC    real    output  variable storing H+ deposition
C    ALFA0NUC   real    output  scav coef for nucleaiton aerosol number
C    ALFA2NUC   real    output  scav coef for nucleation mode sfc area
C    ALFA3NUC   real    output  scav coef for nucleation aerosol mass
C    ALFA0AIT   real    output  scav coef for aitken aerosol number
C    ALFA2AIT   real    output  scav coef for aitken mode sfc area
C    ALFA3AIT   real    output  scav coef for aitken aerosol mass
 
C-----------------------------------------------------------------------

      USE CGRID_SPCS                    ! CGRID mechanism species
      USE UTILIO_DEFN

      IMPLICIT NONE

C...........Includes:

      INCLUDE SUBST_CONST               ! constants

      CHARACTER( 120 ) :: XMSG = ' '    ! Exit status message

C...........Parameters:

      REAL, PARAMETER :: H2ODENS = 1000.0  ! density of water at 20 C
                                           ! and 1 ATM (kg/m3)

C Number of species in CGRID
      INTEGER, SAVE :: MXSPCS

C Number of species scavenged
      INTEGER, SAVE :: N_CGRID_SCAV

      REAL, PARAMETER :: TWOTHIRDS = 2.0 / 3.0

      REAL, PARAMETER :: KGPG = 1.0E-03   ! kilograms per gram

C...........Arguments:

      INTEGER, INTENT(  IN ) :: JDATE   ! current model date, coded YYYYDDD
      INTEGER, INTENT(  IN ) :: JTIME   ! current model time, coded HHMMSS
      REAL,    INTENT(  IN ) :: WTBAR   ! total wat cont (kg/m2) int. thru cloud depth
      REAL,    INTENT(  IN ) :: WCBAR   ! liq water content of cloud (kg/m3)
      REAL,    INTENT(  IN ) :: TBARC   ! mean cloud temp (K)
      REAL,    INTENT(  IN ) :: PBARC   ! mean cloud pressure (Pa)
      REAL,    INTENT(  IN ) :: CTHK1   ! aq chem calc cloud thickness
      REAL,    INTENT(  IN ) :: AIRM    ! total air mass (moles/m2) in cloudy air
      REAL,    INTENT(  IN ) :: PRATE1  ! storm rainfall rate (mm/hr)
      REAL,    INTENT(  IN ) :: TAUCLD  ! cloud lifetime
      REAL,    INTENT( OUT ) :: ALFA0NUC ! scav coef for aitken aerosol number
      REAL,    INTENT( OUT ) :: ALFA2NUC ! scav coef for aitken aerosol sfc area
      REAL,    INTENT( OUT ) :: ALFA3NUC ! scav coef for aitken aerosol mass
      REAL,    INTENT( OUT ) :: ALFA0AIT ! scav coef for aitken aerosol number
      REAL,    INTENT( OUT ) :: ALFA2AIT ! scav coef for aitken aerosol sfc area
      REAL,    INTENT( OUT ) :: ALFA3AIT ! scav coef for aitken aerosol mass
      REAL,    INTENT( OUT ) :: REMOVAC ! variable storing H+ deposition
      REAL,    INTENT(  IN ) :: POLC ( : )  ! avg vert conc incloud (moles/mole)
      REAL,    INTENT( OUT ) :: CEND ( : )  ! ending incloud conc (moles/mole)
      REAL,    INTENT( OUT ) :: REMOV( : )  ! moles/m2 or mm*mol/lit scavenged

C...........Local Variables:

      LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! flag for first pass thru
      CHARACTER( 16 ), SAVE :: PNAME = 'SCAVWDEP'  ! program name

      CHARACTER( 16 ), ALLOCATABLE, SAVE :: CGRID_SCAV( : )  ! CGRID species scavenged
      INTEGER, ALLOCATABLE, SAVE :: CGRID_SCAV_MAP( : ) ! CGRID map to scavenged spc
      REAL, ALLOCATABLE,    SAVE :: CGRID_SCAV_FAC( : )  ! CGRID scav coef factors
      INTEGER, ALLOCATABLE, SAVE :: L_NUMNUC( : ) ! pointers to nucleation aerosol #
      INTEGER, ALLOCATABLE, SAVE :: L_MASNUC( : ) ! pointers to nucleation aerosols
      INTEGER, ALLOCATABLE, SAVE :: L_SRFNUC( : ) ! pntrs to nucleation aerosol surface area
      INTEGER, ALLOCATABLE, SAVE :: L_NUMAIT( : ) ! pointers to aitken aerosol #
      INTEGER, ALLOCATABLE, SAVE :: L_MASAIT( : ) ! pointers to aitken aerosols
      INTEGER, ALLOCATABLE, SAVE :: L_SRFAIT( : ) ! pntrs to aitken aerosol surface area
      INTEGER ASTAT

      INTEGER       I
      INTEGER, SAVE :: N_NUMNUC            ! # nucleation aerosol number species
      INTEGER, SAVE :: N_MASNUC            ! # nucleation aerosol mass species
      INTEGER, SAVE :: N_SRFNUC            ! # nucleation aerosol sfc area species
      INTEGER, SAVE :: N_NUMAIT            ! # aitken aerosol number species
      INTEGER, SAVE :: N_MASAIT            ! # aitken aerosol mass species
      INTEGER, SAVE :: N_SRFAIT            ! # aitken aerosol sfc area species
      INTEGER       PNTR                ! relative pointer variable
      INTEGER       SPC                 ! liquid species loop counter
      INTEGER       VAR                 ! variable loop counter

      REAL          ALFA                ! scavenging coefficient (1/s)
      REAL          KH                  ! Henry's law constant (mol/l/atm)
      REAL          NUMNUC              ! Nucleation mode aerosol # (#/m3)
      REAL          MASNUC              ! Total nucleation mode mass (ug/m3)
      REAL          SRFNUC              ! Nucleation mode surface area (m2/m3)
      REAL          NUMAIT              ! Aitken mode aerosol # (#/m3)
      REAL          MASAIT              ! Total Aitken mode mass (ug/m3)
      REAL          SRFAIT
      REAL          ONE_OVER_TWASH      ! 1 / TWASH
      REAL,    SAVE :: HPLUS = 1.0E-4   ! typical value hydrogen ion concentration [mol/l]
      REAL          RHOAIR              ! air density in kg/m3
      REAL          RTCH                ! chemical gas const times temp
      REAL          TWASH               ! washout time for clouds (sec) with low liq wat content
      REAL          TWF                 ! washout scaling factor (mol/l/atm)

C...........External Functions:

      INTEGER, EXTERNAL :: INDEXN
      REAL,    EXTERNAL :: HLCONST

C-----------------------------------------------------------------------

C...INITIALIZATION SCAVWDEP module:
C...  event-statistics variables.  Open output files.

      IF ( FIRSTIME ) THEN

        FIRSTIME = .FALSE.

        N_CGRID_SCAV = N_GC_SCAV + N_AE_SCAV + N_NR_SCAV + N_TR_SCAV

C...first check to make sure that some species in CGRID were specified
C...  for scavenging, otherwise notify the user and return

        IF ( N_CGRID_SCAV .LE. 0 ) THEN
          XMSG = 'No species were specified for scavenging by cloud ' //
     &           'or rain water...SCAVENGING WILL NOT BE PERFORMED!'
          CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
          RETURN
        END IF

        ALLOCATE ( CGRID_SCAV    ( N_CGRID_SCAV ),
     &             CGRID_SCAV_MAP( N_CGRID_SCAV ),
     &             CGRID_SCAV_FAC( N_CGRID_SCAV ), STAT = ASTAT )
        IF ( ASTAT .NE. 0 ) THEN
          XMSG = 'Failure allocating CGRID_SCAV, CGRID_SCAV_MAP or CGRID_SCAV_FAC'
          CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
        END IF

        MXSPCS = N_GC_SPCD + N_AE_SPC + N_NR_SPC + N_TR_SPC

        ALLOCATE ( L_NUMNUC( MXSPCS ),
     &             L_MASNUC( MXSPCS ),
     &             L_SRFNUC( MXSPCS ),
     &             L_NUMAIT( MXSPCS ),
     &             L_MASAIT( MXSPCS ),
     &             L_SRFAIT( MXSPCS ), STAT = ASTAT )
        IF ( ASTAT .NE. 0 ) THEN
          XMSG = 'Failure allocating L_NUMAIT, L_MASAIT or L_SRFAIT'
          CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
        END IF

C...prepare indices for scavenged species

C...  load the CGRID to scavenged species pointers for the gases

        SPC = 0
        DO VAR=1, N_GC_SCAV
          SPC = SPC + 1
          CGRID_SCAV( SPC ) = GC_SCAV( VAR )
          CGRID_SCAV_MAP( SPC ) = GC_SCAV_MAP( VAR ) + GC_STRT - 1
          CGRID_SCAV_FAC( SPC ) = GC_SCAV_FAC( VAR )
        END DO

C...  load the CGRID to scavenged species pointers for the aerosols

        DO VAR=1, N_AE_SCAV
          SPC = SPC + 1
          CGRID_SCAV( SPC ) = AE_SCAV( VAR )
          CGRID_SCAV_MAP( SPC ) = AE_SCAV_MAP( VAR ) + AE_STRT - 1
          CGRID_SCAV_FAC( SPC ) = AE_SCAV_FAC( VAR )
        END DO

C...  load the CGRID to scavenged species pointers for the non-reactives

        DO VAR=1, N_NR_SCAV
          SPC = SPC + 1
          CGRID_SCAV( SPC ) = NR_SCAV( VAR )
          CGRID_SCAV_MAP( SPC ) = NR_SCAV_MAP( VAR ) + NR_STRT - 1
          CGRID_SCAV_FAC( SPC ) = NR_SCAV_FAC( VAR )
        END DO

C...  load the CGRID to scavenged species pointers for the tracers

        DO VAR=1, N_TR_SCAV
          SPC = SPC + 1
          CGRID_SCAV( SPC ) = TR_SCAV( VAR )
          CGRID_SCAV_MAP( SPC ) = TR_SCAV_MAP( VAR ) + TR_STRT - 1
          CGRID_SCAV_FAC( SPC ) = TR_SCAV_FAC( VAR )
        END DO

C...create the pointers from CGRID to the species needed by AQCHEM

        N_NUMNUC = INDEXN( 'NUM_NUCLN ', N_CGRID_SCAV, CGRID_SCAV,
     &                     L_NUMNUC )
        N_MASNUC = INDEXN( 'NUCLN ', N_CGRID_SCAV, CGRID_SCAV, 
     &                     L_MASNUC )
        N_SRFNUC = INDEXN( 'SRF_NUCLN ', N_CGRID_SCAV, CGRID_SCAV,
     &                     L_SRFNUC )
        N_NUMAIT = INDEXN( 'NUM_AITKEN ', N_CGRID_SCAV, CGRID_SCAV,
     &                     L_NUMAIT )
        N_MASAIT = INDEXN( 'AITKEN ', N_CGRID_SCAV, CGRID_SCAV, 
     &                     L_MASAIT )
        N_SRFAIT = INDEXN( 'SRF_AITKEN ', N_CGRID_SCAV, CGRID_SCAV,
     &                     L_SRFAIT )

      END IF

C...for subsequent calls, check to make sure some species were
C...  specified, otherwise there is no need to perform scavenging

      IF ( N_CGRID_SCAV .LE. 0 ) THEN
        RETURN
      END IF

      RTCH = ( MOLVOL / STDTEMP ) * TBARC
      TWASH = WTBAR * 1000.0 * CTHK1 * 3600.0
     &      / ( H2ODENS * AMAX1( 1.0E-20, PRATE1 ) )
ccc          TWASH = AMAX1( TWASH, TAUCLD / 60.0 )   ! <different units?sec&min
      TWASH = AMAX1( TWASH, TAUCLD )
      ONE_OVER_TWASH = 1.0 / TWASH
      TWF = H2ODENS / ( WTBAR * RTCH )
      REMOVAC = 0.0

      RHOAIR = ( AIRM / CTHK1 ) * MWAIR * KGPG

C...initialize scavenging variables
      ALFA0NUC = 0.0
      ALFA2NUC = 0.0
      ALFA3NUC = 0.0
      ALFA0AIT = 0.0
      ALFA2AIT = 0.0
      ALFA3AIT = 0.0

C...compute total Nucleation mode number (#/m3)

      NUMNUC = 0.0

      DO I = 1, N_NUMNUC
        PNTR = CGRID_SCAV_MAP( L_NUMNUC( I ) )
        NUMNUC = NUMNUC + ( POLC( PNTR ) * AIRM / CTHK1 )
      END DO

C...compute total Nucleation mode mass (ug/m3)

      MASNUC = 0.0

      DO I = 1, N_MASNUC
        PNTR = CGRID_SCAV_MAP( L_MASNUC( I ) )
        IF (( INDEX( CGRID_SCAV( L_MASNUC( I ) ), 'NUM' ) .EQ. 0 ) .AND.
     &      ( INDEX( CGRID_SCAV( L_MASNUC( I ) ), 'SRF' ) .EQ. 0 ) .AND.
     &      ( INDEX( CGRID_SCAV( L_MASNUC( I ) ), 'H2O' ) .EQ. 0 ) .AND.
     &      ( INDEX( CGRID_SCAV( L_MASNUC( I ) ), 'TRACER' ) .EQ. 0 ) ) THEN
          MASNUC = MASNUC + ( POLC( PNTR ) * AIRM / CTHK1
     &           * AE_MOLWT( PNTR - AE_STRT + 1 ) * 1.0E6 )
        END IF
      END DO

C...compute total Nucleation mode surface area (m2/m3)

      SRFNUC = 0.0

      DO I = 1, N_SRFNUC
        PNTR = CGRID_SCAV_MAP( L_SRFNUC( I ) )
        SRFNUC = SRFNUC + ( POLC( PNTR ) * AIRM / CTHK1 )
      END DO
 
C...compute the scavenging coefficients for aitken mode aerosol mass
C...  and number
C...  NOTE:  for now, scavenging coefficients are computed for only
C...         the liquid water content, not on the total water content
C...         therefore, no ice phase scavenging is considered at this
C...         time, but it should be added in the future!

      CALL GETALPHA ( NUMNUC, MASNUC, SRFNUC, WCBAR, TBARC, PBARC,
     &                RHOAIR, ALFA0NUC, ALFA2NUC, ALFA3NUC )
 

C...compute total Aitken mode number (#/m3)

      NUMAIT = 0.0

      DO I = 1, N_NUMAIT
        PNTR = CGRID_SCAV_MAP( L_NUMAIT( I ) )
        NUMAIT = NUMAIT + ( POLC( PNTR ) * AIRM / CTHK1 )
      END DO

C...compute total Aitken mode mass (ug/m3)

      MASAIT = 0.0

      DO I = 1, N_MASAIT
        PNTR = CGRID_SCAV_MAP( L_MASAIT( I ) )
        IF (( INDEX( CGRID_SCAV( L_MASAIT( I ) ), 'NUM' ) .EQ. 0 ) .AND.
     &      ( INDEX( CGRID_SCAV( L_MASAIT( I ) ), 'SRF' ) .EQ. 0 ) .AND.
     &      ( INDEX( CGRID_SCAV( L_MASAIT( I ) ), 'H2O' ) .EQ. 0 ) .AND.
     &      ( INDEX( CGRID_SCAV( L_MASAIT( I ) ), 'TRACER' ) .EQ. 0 ) ) THEN
          MASAIT = MASAIT + ( POLC( PNTR ) * AIRM / CTHK1
     &           * AE_MOLWT( PNTR - AE_STRT + 1 ) * 1.0E6 )
        END IF
      END DO

C...compute total Aitken mode surface area (m2/m3)

      SRFAIT = 0.0

      DO I = 1, N_SRFAIT
        PNTR = CGRID_SCAV_MAP( L_SRFAIT( I ) )
        SRFAIT = SRFAIT + ( POLC( PNTR ) * AIRM / CTHK1 )
      END DO

C...compute the scavenging coefficients for aitken mode aerosol mass
C...  and number
C...  NOTE:  for now, scavenging coefficients are computed for only
C...         the liquid water content, not on the total water content
C...         therefore, no ice phase scavenging is considered at this
C...         time, but it should be added in the future!

      CALL GETALPHA ( NUMAIT, MASAIT, SRFAIT, WCBAR, TBARC, PBARC,
     &                RHOAIR, ALFA0AIT, ALFA2AIT, ALFA3AIT )

C...gas scavenging and wet deposition

      SPC = 0

      DO VAR = 1, N_GC_SCAV
        SPC = SPC + 1
        PNTR = CGRID_SCAV_MAP( SPC )
        KH = HLCONST( CGRID_SCAV( SPC ), TBARC, .TRUE., HPLUS )
        IF ( KH .GT. 0.0 ) THEN
          ALFA = CGRID_SCAV_FAC( SPC ) * ONE_OVER_TWASH / ( 1.0 + TWF / KH )
          CEND ( PNTR ) = POLC( PNTR ) * EXP( -ALFA * TAUCLD )
          REMOV( PNTR ) = ( POLC( PNTR ) - CEND( PNTR ) ) * AIRM
        ELSE
          ALFA = 0.0
          CEND ( PNTR ) = POLC( PNTR )
          REMOV( PNTR ) = 0.0
        END IF
      END DO

C...aerosol scavenging and wet deposition

      DO VAR = 1, N_AE_SCAV
        SPC = SPC + 1

        PNTR = CGRID_SCAV_MAP( SPC )

        IF ( INDEX( CGRID_SCAV( SPC ), 'NUCLN' ) .GT. 0 ) THEN
          IF ( INDEX( CGRID_SCAV( SPC ), 'NUM' ) .GT. 0 ) THEN
            ALFA = CGRID_SCAV_FAC( SPC ) * ALFA0NUC
          ELSE IF ( INDEX( CGRID_SCAV( SPC ), 'SRF' ) .GT. 0 ) THEN
            ALFA = CGRID_SCAV_FAC( SPC ) * ALFA2NUC
          ELSE
            ALFA = CGRID_SCAV_FAC( SPC ) * ALFA3NUC
          END IF 
        ELSEIF ( INDEX( CGRID_SCAV( SPC ), 'AITKEN' ) .GT. 0 ) THEN
          IF ( INDEX( CGRID_SCAV( SPC ), 'NUM' ) .GT. 0 ) THEN
            ALFA = CGRID_SCAV_FAC( SPC ) * ALFA0AIT
          ELSE IF ( INDEX( CGRID_SCAV( SPC ), 'SRF' ) .GT. 0 ) THEN
            ALFA = CGRID_SCAV_FAC( SPC ) * ALFA2AIT
          ELSE
            ALFA = CGRID_SCAV_FAC( SPC ) * ALFA3AIT
          END IF
        ELSE
csjr          IF ( INDEX( CGRID_SCAV( SPC ), 'SRF' ) .EQ. 0 ) THEN
            ALFA = CGRID_SCAV_FAC( SPC ) * ONE_OVER_TWASH
csjr          ELSE
csjr            ALFA = CGRID_SCAV_FAC( SPC ) * TWOTHIRDS * ONE_OVER_TWASH
csjr          END IF
        END IF

        CEND ( PNTR ) = POLC( PNTR ) * EXP( -ALFA * TAUCLD )
        REMOV( PNTR ) = ( POLC( PNTR ) - CEND( PNTR ) ) * AIRM

      END DO

C...non-reactive scavenging and wet deposition

      DO VAR = 1, N_NR_SCAV
        SPC = SPC + 1
        PNTR = CGRID_SCAV_MAP( SPC )
        KH = HLCONST( CGRID_SCAV( SPC ), TBARC, .TRUE., HPLUS )
        IF ( KH .GT. 0.0 ) THEN
          ALFA = CGRID_SCAV_FAC( SPC ) * ONE_OVER_TWASH / ( 1.0 + TWF / KH )
          CEND ( PNTR ) = POLC( PNTR ) * EXP( -ALFA * TAUCLD )
          REMOV( PNTR ) = ( POLC( PNTR ) - CEND( PNTR ) ) * AIRM
        ELSE
          ALFA = 0.0
          CEND ( PNTR ) = POLC( PNTR )
          REMOV( PNTR ) = 0.0
        END IF
      END DO

C...tracer scavenging and wet deposition

      DO VAR = 1, N_TR_SCAV
        SPC = SPC + 1
        PNTR = CGRID_SCAV_MAP( SPC )
        KH = HLCONST( CGRID_SCAV( SPC ), TBARC, .TRUE., HPLUS )
        IF ( KH .GT. 0.0 ) THEN
          ALFA = CGRID_SCAV_FAC( SPC ) * ONE_OVER_TWASH / ( 1.0 + TWF / KH )
          CEND ( PNTR ) = POLC( PNTR ) * EXP( -ALFA * TAUCLD )
          REMOV( PNTR ) = ( POLC( PNTR ) - CEND( PNTR ) ) * AIRM
        ELSE
          ALFA = 0.0
          CEND ( PNTR ) = POLC( PNTR )
          REMOV( PNTR ) = 0.0
        END IF
      END DO

      RETURN
      END
