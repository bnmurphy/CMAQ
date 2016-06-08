
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
C $Header: /project/yoj/arc/CCTM/src/aero/aero5/aero_depv.F,v 1.12 2012/01/19 13:12:14 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)aero_depv.F       1.3 /project/mod3/CMAQ/src/ae_depv/aero_depv/SCCS/s.aero_depv.F 18 Jun 1997 12:55:48

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE AERO_DEPV ( CGRID, JDATE, JTIME, TSTEP, MOSAIC, VDEP_AE )

C-----------------------------------------------------------------------
C aerosol dry deposition routine
C   written 4/9/97 by Dr. Francis S. Binkowski
C   uses code from modpar and vdvg from the aerosol module.
C   This routine uses a single block to hold information
C   for the lowest layer.
C NOTES: This version assumes that RA is available on the met file.
c        Array structure for vector optimization
C 26 Apr 97 Jeff - many mods
C 13 Dec 97 Jeff - expect uncoupled CGRID, concs as micro-g/m**3, #/m**3
C
C 1/11/99 David Wong at LM - change NUMCELLS to CELLNUM in the loop index
C FSB 3/17/99 changed to accommodate surface area/second moment and
C    encapsulated the actual drydep calculation into a subroutine which
C    is attached to this code
C Jeff - Dec 00 - move CGRID_MAP into f90 module
C FSB 12/11/2000. Logic added to allow deposition of particles at their
C     "wet" diameters; that is, accounting for the water on the particles.
C     This is done by adjusting the third and second moments for the
C     presence of water assuming that the geometric standard deviations
C     are not changed by this process. This appears to be a very good
C     assumption.
C 30 Aug 01 J.Young: Dyn alloc; Use HGRD_DEFN
C    Jan 03 J.Young: Change CGRID dimensions, eliminate re-allocations
C  6 Mar 03 J.Young: eliminate a lot of allocate/deallocates
C  7 Aug 03 S.Roselle: updated code for loading the min aero conc array
C 17 Dec 03 S.Roselle: Adjust 2nd and 3rd moments to include SOA,
C     without affecting the geometric standard deviations.
C 31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C     domain specifications in one module
C 07 Jun 05 P.Bhave: Added code to handle new species in the AE4
C     mechanism: ANAI, ANAJ, ANAK, ACLI, ACLJ, ACLK, ASO4K, AH2OK,
C     and ANO3K; look for ASEAS only when using AE3 mechanism
C 30 Jan 08 S.Napelenok & P.Bhave: Added code to handle new SOA species
C     in AE5; defined DRY aerosol to include nonvolatile SOA spcs
C 14 Apr 08 J.Kelly: Added code to handle new species ANH4K and SRFCOR.
C     Also added code to handle variable coarse mode standard deviation
C     in AE5 (no longer fixed at 2.2).
C 08 Sep 08 P.Bhave: Backward compatibility with AE4 mechanisms
C     standardized names of all coarse-mode variables
C 19 Apr 10 S.Howard: aero re-engineering for modularity
C 23 Apr 10 J.Young: replace chem mechanism include files with namelists
C 10 Mar 11 S.Howard: Renamed met_data to aeromet_data
C 25 Mar 11 S.Roselle: Replaced I/O API include files with UTILIO_DEFN
C 20 May 11 D.Schwede: Modified for mosaic
C 31 Aug 11 J.Bash: Moved shared mosaic variables to MOSAIC_MOD
C 27 Sep 11 David Wong: replaced all run time dynamic arrays with allocatable
C                       arrays to avoid run time memory issue  
C 08 Jun 12 J.Young: remove full character blank padding for GNU Fortran (GCC) 4.1.2
C 07 Jul 14 B.Hutzell: replaced mechanism include file(s) with fortran module
C 07 Nov 14 J.Bash: Updated for the ASX_DATA_MOD shared data module. 
C    May 16 B. Murphy, H. Pye: Updated treatment of aerosol moments
C-----------------------------------------------------------------------

      USE GRID_CONF           ! horizontal & vertical domain specifications
      USE RXNS_DATA           ! chemical mechanism data
      USE CGRID_SPCS          ! CGRID mechanism species
      USE UTILIO_DEFN      
      USE AERO_DATA           ! aero variable data
      USE AEROMET_DATA        ! Includes CONST.EXT
      USE Mosaic_Mod, Only: ADEPVJ  ! Shared mosaic variables
      Use LSM_Mod, Only: N_LUFRAC
      USE ASX_DATA_MOD, Only: Met_Data

      IMPLICIT NONE

C Includes:

      INCLUDE SUBST_FILES_ID  ! file name parameters

C Arguments
      REAL,    POINTER       :: CGRID( :,:,:,: )
      INTEGER, INTENT( IN )  :: JDATE               ! current model date , coded YYYYDDD
      INTEGER, INTENT( IN )  :: JTIME               ! current model time , coded HHMMSS
      INTEGER, INTENT( IN )  :: TSTEP               ! model time step, coded HHMMSS
      LOGICAL, INTENT( IN )  :: MOSAIC              ! use mosaic option
      REAL,    INTENT( OUT ) :: VDEP_AE( :,:,: )    ! surrogate deposition velocities [ m s**-1 ]

C Parameters
      REAL,    PARAMETER :: T0 = 288.15      ! [ K ] ! starting standard surface temp.
      REAL,    PARAMETER :: TWO3 = 2.0 / 3.0
      INTEGER, PARAMETER :: N_AE_DEP_SPC = 12 ! no. of surrogates for aerosol dry dep velocities

C set up species dimension and indices for deposition velocity internal array VDEP
      INTEGER, PARAMETER :: VDNNUC = 1,  ! NPF mode number
     &                      VDNATK = 2,  ! Aitken mode number
     &                      VDNACC = 3,  ! Accumulation mode number
     &                      VDNCOR = 4,  ! Coarse mode number
     &                      VDMNUC = 5,  ! NPF mode mass
     &                      VDMATK = 6,  ! Aitken mode mass
     &                      VDMACC = 7,  ! Accumulation mode mass
     &                      VDMCOR = 8,  ! Coarse mode mass
     &                      VDSNUC = 9,  ! NPF mode surface area
     &                      VDSATK = 10, ! Aitken mode surface area
     &                      VDSACC = 11, ! Accumulation mode surface area
     &                      VDSCOR = 12  ! Coarse mode surface area

C Local variables:

      CHARACTER( 16 ) :: VDAE_NAME( N_AE_DEP_SPC )! dep vel surrogate name table
      DATA         VDAE_NAME( 1 ) / 'VNUMNUC' /
      DATA         VDAE_NAME( 2 ) / 'VNUMATKN' /
      DATA         VDAE_NAME( 3 ) / 'VNUMACC ' /
      DATA         VDAE_NAME( 4 ) / 'VNUMCOR ' /
      DATA         VDAE_NAME( 5 ) / 'VMASSH  ' /
      DATA         VDAE_NAME( 6 ) / 'VMASSI  ' /
      DATA         VDAE_NAME( 7 ) / 'VMASSJ  ' /
      DATA         VDAE_NAME( 8 ) / 'VMASSC  ' /
      DATA         VDAE_NAME( 9 ) / 'VSRFNUC' /
      DATA         VDAE_NAME(10 ) / 'VSRFATKN' /
      DATA         VDAE_NAME(11 ) / 'VSRFACC ' /
      DATA         VDAE_NAME(12 ) / 'VSRFCOR ' /

      INTEGER, ALLOCATABLE, SAVE :: DEPV_SUR( : )   ! pointer to surrogate

C Meteorological variables

      CHARACTER( 16 ), SAVE :: AE_VRSN ! Aerosol version name

      INTEGER, SAVE :: NCELLS              ! number of cells per layer

      REAL, ALLOCATABLE, SAVE  :: XXLSGNU( :,: )    ! log of standard deviation
      REAL, ALLOCATABLE, SAVE  :: XXLSGAT( :,: )
      REAL, ALLOCATABLE, SAVE  :: XXLSGAC( :,: )
      REAL, ALLOCATABLE, SAVE  :: XXLSGCO( :,: )
 
      REAL, ALLOCATABLE, SAVE  :: DGNUC( :,: )      ! geometric mean diameter
      REAL, ALLOCATABLE, SAVE  :: DGATK( :,: )
      REAL, ALLOCATABLE, SAVE  :: DGACC( :,: )
      REAL, ALLOCATABLE, SAVE  :: DGCOR( :,: )

      REAL, ALLOCATABLE, SAVE  :: PDENSNU( :,: )    ! particle density         
      REAL, ALLOCATABLE, SAVE  :: PDENSAT( :,: )
      REAL, ALLOCATABLE, SAVE  :: PDENSAC( :,: )
      REAL, ALLOCATABLE, SAVE  :: PDENSCO( :,: )

      REAL, ALLOCATABLE, SAVE  :: XLM( :,: )        ! mean free path [ m ]
      REAL, ALLOCATABLE, SAVE  :: AMU( :,: )        ! dynamic viscosity [ kg m**-1 s**-1 ]

      REAL, ALLOCATABLE, SAVE :: VDEP( :,:,: )    ! deposition  velocity [ m/s ]
      REAL, ALLOCATABLE, SAVE :: VDEPJ( :,:,:,: ) ! deposition  velocity [ m/s ]

      REAL M3_WET, M3SUBT, M3_DRY
      REAL M2_WET, M2_DRY

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV                   ! unit number for the log file
      CHARACTER( 16 ), SAVE :: PNAME = 'AERO_DEPV'
      CHARACTER( 16 ) :: VNAME            ! varable name
      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER  C, R, V, N, J         ! loop counters
      INTEGER  SPC, S                ! species loop counter
      INTEGER  ALLOCSTAT

      INTERFACE
        SUBROUTINE GETDEP_V ( XLM, AMU, DGNUC, DGATK, DGACC, DGCOR,
     &                        XXLSGNU, XXLSGAT, XXLSGAC, XXLSGCO,
     &                        PDENSNU, PDENSAT, PDENSAC, PDENSCO,
     &                        VDEP, VDEPJ, MOSAIC )
          REAL, INTENT( IN ) :: XLM( :,: )      ! atmospheric mean free path [ m ]
          REAL, INTENT( IN ) :: AMU( :,: )      ! atmospheric dynamic viscosity [ kg/(m s) ]
          REAL, INTENT( IN ) :: DGNUC( :,: )    ! nuclei mode geometric mean diameter  [ m ]
          REAL, INTENT( IN ) :: DGATK( :,: )    ! aitken mode geometric mean diameter  [ m ]
          REAL, INTENT( IN ) :: DGACC( :,: )    ! accumulation geometric mean diameter [ m ]
          REAL, INTENT( IN ) :: DGCOR( :,: )    ! coarse mode geometric mean diameter  [ m ]
          REAL, INTENT( IN ) :: XXLSGNU( :,: )  ! Nucleation mode
          REAL, INTENT( IN ) :: XXLSGAT( :,: )  ! Aitken mode
          REAL, INTENT( IN ) :: XXLSGAC( :,: )  ! accumulation mode
          REAL, INTENT( IN ) :: XXLSGCO( :,: )  ! coarse mode
          REAL, INTENT( IN ) :: PDENSNU( :,: )  ! average particle density in nuclei mode
          REAL, INTENT( IN ) :: PDENSAT( :,: )  ! average particle density in nuclei mode
          REAL, INTENT( IN ) :: PDENSAC( :,: )  ! average particle density in accumulation mode
          REAL, INTENT( IN ) :: PDENSCO( :,: )  ! average particle density in coarse mode
          REAL, INTENT( OUT ) :: VDEP( :,:,: ) ! deposition  velocity [ m/s ]
          REAL, INTENT( OUT ) :: VDEPJ( :,:,:,: ) ! deposition  velocity [ m/s ] for each land use category
          LOGICAL, INTENT( IN ) :: MOSAIC
        END SUBROUTINE GETDEP_V
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()

         NCELLS = NCOLS * NROWS

C  Allocate arrays
         ALLOCATE( XXLSGNU( NCOLS,NROWS ),
     &             XXLSGAT( NCOLS,NROWS ),
     &             XXLSGAC( NCOLS,NROWS ),
     &             XXLSGCO( NCOLS,NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating XXLSGNU, XXLSGAT, XXLSGAC or XXLSGCO'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         ALLOCATE( DGNUC( NCOLS,NROWS ),
     &             DGATK( NCOLS,NROWS ),
     &             DGACC( NCOLS,NROWS ),
     &             DGCOR( NCOLS,NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DGNUC, DGATK, DGACC or DGCOR'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         ALLOCATE( PDENSNU( NCOLS,NROWS ),
     &             PDENSAT( NCOLS,NROWS ),
     &             PDENSAC( NCOLS,NROWS ),
     &             PDENSCO( NCOLS,NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating PDENSBU, PDENSAT, PDENSAC or PDENSCO'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         ALLOCATE( XLM( NCOLS,NROWS ),
     &             AMU( NCOLS,NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating XLM or AMU'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         ALLOCATE ( VDEP( NCOLS,NROWS,N_AE_DEP_SPC ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating VDEP'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

            ALLOCATE ( VDEPJ( N_LUFRAC,NCOLS,NROWS,N_AE_DEP_SPC ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating VDEPJ'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         IF ( MOSAIC ) THEN
            ALLOCATE ( ADEPVJ( N_LUFRAC,N_AE_DEPV,NCOLS,NROWS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating ADEV'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         END IF

         ALLOCATE ( DEPV_SUR( N_AE_DEPV ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DEPV_SUR'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C Set the dep vel surrogate pointers
         DO V = 1, N_AE_DEPV
            N = INDEX1( AE_DEPV( V ), N_AE_DEP_SPC, VDAE_NAME )
            IF ( N .NE. 0 ) THEN
               DEPV_SUR( V ) = N
            ELSE
               XMSG = 'Could not find ' // AE_DEPV( V ) // ' in aerosol' //
     &                ' surrogate table. >>> Dep vel set to zero <<< '
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
               DEPV_SUR( V ) = 0
            END IF
         END DO

      END IF    ! FIRSTIME      

      IF ( N_AE_SPC .LE. 0 ) RETURN

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Put the grid cell physical data in the block arrays
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DO R = 1, MY_NROWS
      DO C = 1, MY_NCOLS

C ***    Set meteorological data for the grid cell.
         AIRDENS = Met_Data%DENS1( C,R )
         AIRTEMP = Met_Data%TEMP2( C,R )
         AIRPRES = Met_Data%PRSFC( C,R )

C ***    extract grid cell concentrations of aero species from CGRID
C        into aerospc_conc in aero_data module
C        Also determines second moment from surface area and adds wet
C        species
         CALL EXTRACT_AERO( CGRID( C,R,1,: ), .TRUE. )

C ***    Calculate geometric mean diameters and standard deviations of the
C        "wet" size distribution
         CALL GETPAR( .FALSE. )     

C        Save getpar values to arrays
         XXLSGNU( C,R ) = AEROMODE_LNSG( 1 )
         XXLSGAT( C,R ) = AEROMODE_LNSG( 2 )
         XXLSGAC( C,R ) = AEROMODE_LNSG( 3 )
         XXLSGCO( C,R ) = AEROMODE_LNSG( 4 )

         DGNUC( C,R )   = AEROMODE_DIAM( 1 )
         DGATK( C,R )   = AEROMODE_DIAM( 2 )
         DGACC( C,R )   = AEROMODE_DIAM( 3 )
         DGCOR( C,R )   = AEROMODE_DIAM( 4 )

         PDENSNU( C,R ) = AEROMODE_DENS( 1 )
         PDENSAT( C,R ) = AEROMODE_DENS( 2 )
         PDENSAC( C,R ) = AEROMODE_DENS( 3 )
         PDENSCO( C,R ) = AEROMODE_DENS( 4 )
 
C        Calculate mean free path [ m ]:
         XLM( C,R ) = 6.6328E-8 * STDATMPA * AIRTEMP / ( T0 * AIRPRES )

C ***    Calcualte dynamic viscosity [ kg m**-1 s**-1 ]:
         AMU( C,R ) = 1.458E-6 * AIRTEMP * SQRT( AIRTEMP )
     &              / ( AIRTEMP + 110.4 )

      END DO ! Column LOOP
      END DO   ! Row LOOP

C *** get dry deposition velocities:

      CALL GETDEP_V ( XLM, AMU, DGNUC, DGATK, DGACC, DGCOR,
     &                XXLSGNU, XXLSGAT, XXLSGAC, XXLSGCO,
     &                PDENSNU, PDENSAT, PDENSAC, PDENSCO,
     &                VDEP, VDEPJ, MOSAIC )

C Return dry deposition velocities for aerosols (first layer only).

      DO R = 1, NROWS
         DO C = 1, NCOLS
            DO V = 1, N_AE_DEPV
               IF ( DEPV_SUR( V ) .GT. 0 ) THEN
                  VDEP_AE( V,C,R ) = VDEP( C,R,DEPV_SUR( V ) )
               ELSE
                  VDEP_AE( V,C,R ) = 0.0
               END IF
            END DO
         END DO
      END DO

      IF ( MOSAIC ) THEN
        DO R = 1, NROWS
           DO C = 1, NCOLS
              DO V = 1, N_AE_DEPV
                DO J = 1, N_LUFRAC
                  IF ( DEPV_SUR( V ) .GT. 0 ) THEN
                     ADEPVJ( J,V,C,R ) = VDEPJ( J,C,R,DEPV_SUR( V ) )
                  ELSE
                     ADEPVJ( J,V,C,R ) = 0.0
                  END IF
                END DO
              END DO
           END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE AERO_DEPV

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE GETDEP_V ( XLM, AMU, DGNUC, DGATK, DGACC, DGCOR,
     &                      XXLSGNU, XXLSGAT, XXLSGAC, XXLSGCO,
     &                      PDENSNU, PDENSAT, PDENSAC, PDENSCO,
     &                      VDEP, VDEPJ, MOSAIC )

C *** Calculate deposition velocity for Aitken, accumulation, and
C     coarse modes.
C     Reference:
C     Binkowski F. S., and U. Shankar, The regional particulate
C     model 1. Model description and preliminary results.
C     J. Geophys. Res., 100, D12, 26191-26209, 1995.
 
C    May 05 D.Schwede: added impaction term to coarse mode dry deposition
C 25 May 05 J.Pleim:  Updated dry dep velocity calculation for aerosols
C                     to Venkatram and Pleim (1999)
C 20 Jul 05 J.Pleim:  Changed impaction term using modal integration of
C                     Stokes**2 / 400 (Giorgi, 1986, JGR)
C 14 Apr 08 J.Kelly:  Added code to calculate deposition velocity of
C                     coarse surface area and to account for variable
C                     standard deviation of the coarse mode.
C 08 Sep 08 P.Bhave:  Backward compatibility with AE4 mechanisms
C                     standardized names of all coarse-mode variables
C 18 May 15 B.Murphy: Added Nucleation mode aerosols
C-----------------------------------------------------------------------

      USE LSM_MOD, Only: n_lufrac      ! to get n_lufrac
      USE AEROMET_DATA   ! Includes CONST.EXT
      USE ASX_DATA_MOD, Only: Met_Data, Mosaic_Data, Grid_Data

      IMPLICIT NONE

C *** input arguments

C     atmospheric properties
      REAL, INTENT( IN ) :: XLM( :,: )      ! atmospheric mean free path [ m ]
      REAL, INTENT( IN ) :: AMU( :,: )      ! atmospheric dynamic viscosity [ kg/(m s) ]

C     aerosol properties:

C     modal diameters: [ m ]
      REAL, INTENT( IN ) :: DGNUC( :,: )    ! nuclei mode geometric mean diameter  [ m ]
      REAL, INTENT( IN ) :: DGATK( :,: )    ! aitken mode geometric mean diameter  [ m ]
      REAL, INTENT( IN ) :: DGACC( :,: )    ! accumulation geometric mean diameter [ m ]
      REAL, INTENT( IN ) :: DGCOR( :,: )    ! coarse mode geometric mean diameter  [ m ]

C     log of modal geometric standard deviations
      REAL, INTENT( IN ) :: XXLSGNU( :,: )  ! Nucleation mode
      REAL, INTENT( IN ) :: XXLSGAT( :,: )  ! Aitken mode
      REAL, INTENT( IN ) :: XXLSGAC( :,: )  ! accumulation mode
      REAL, INTENT( IN ) :: XXLSGCO( :,: )  ! coarse mode

C     average modal particle densities  [ kg/m**3 ]
      REAL, INTENT( IN ) :: PDENSNU( :,: )  ! average particle density in nuclei mode
      REAL, INTENT( IN ) :: PDENSAT( :,: )  ! average particle density in aitken mode
      REAL, INTENT( IN ) :: PDENSAC( :,: )  ! average particle density in accumulation mode
      REAL, INTENT( IN ) :: PDENSCO( :,: )  ! average particle density in coarse mode

C     Mosaic variables
      LOGICAL, INTENT( IN ) :: MOSAIC

C *** output arguments

      REAL, INTENT( OUT ) :: VDEP( :,:,: ) ! deposition  velocity [ m/s ]
      REAL, INTENT( OUT ) :: VDEPJ( :,:,:,: ) ! deposition  velocity [ m/s ] for each land use category

C *** array indices hardcoded to match SUBROUTINE AERO_DEPV
      INTEGER, PARAMETER :: VDNNUC = 1,  ! Nucleation mode number
     &                      VDNATK = 2,  ! Aitken mode number
     &                      VDNACC = 3,  ! accumulation mode number
     &                      VDNCOR = 4,  ! coarse mode number
     &                      VDMNUC = 5,  ! Nucleation mode mass
     &                      VDMATK = 6,  ! Aitken mode mass
     &                      VDMACC = 7,  ! accumulation mode mass
     &                      VDMCOR = 8,  ! coarse mode mass
     &                      VDSNUC = 9,  ! Nucleation mode surface area
     &                      VDSATK =10,  ! Aitken mode surface area
     &                      VDSACC =11,  ! accumulation mode surface area
     &                      VDSCOR =12   ! coarse mode surface area

C modal Knudsen numbers
      REAL KNNUC   ! Nucleation mode Knudsen number
      REAL KNATK   ! Aitken mode Knudsen number
      REAL KNACC   ! accumulation "
      REAL KNCOR   ! coarse mode

C modal particle diffusivities for number, 2nd, and 3rd moment, or mass:
      REAL DCHAT0NU, DCHAT2NU, DCHAT3NU
      REAL DCHAT0AT, DCHAT2AT, DCHAT3AT
      REAL DCHAT0AC, DCHAT2AC, DCHAT3AC
      REAL DCHAT0CO, DCHAT2CO, DCHAT3CO

C modal sedimentation velocities for number, 2nd, and 3rd moment, or mass:
      REAL VGHAT0NU, VGHAT2NU, VGHAT3NU
      REAL VGHAT0AT, VGHAT2AT, VGHAT3AT
      REAL VGHAT0AC, VGHAT2AC, VGHAT3AC
      REAL VGHAT0CO, VGHAT2CO, VGHAT3CO

      INTEGER NCELL, J, C, R

      REAL DCONST1, DCONST1NU, DCONST1AT, DCONST1AC, DCONST1CO
      REAL DCONST2, DCONST3NU, DCONST3AT, DCONST3AC, DCONST3CO
      REAL SC0NU, SC0AT, SC0AC, SC0CO     ! Schmidt numbers for number
      REAL SC2NU, SC2AT, SC2AC, SC2CO     ! Schmidt numbers for 2ND MOMENT
      REAL SC3NU, SC3AT, SC3AC, SC3CO     ! Schmidt numbers for 3rd moment
      REAL STOKENU, STOKEAT, STOKEAC, STOKECO ! Stokes numbers for each mode
      REAL RD0NU, RD0AT, RD0AC, RD0CO     ! canopy resistance for number
      REAL RD2NU, RD2AT, RD2AC, RD2CO     ! canopy resistance for 2nd moment
      REAL RD3NU, RD3AT, RD3AC, RD3CO     ! canopy resisteance for 3rd moment
      REAL UTSCALE              ! scratch function of USTAR and WSTAR
      REAL NU                   ! kinematic viscosity [ m**2 s**-1 ]
      REAL USTFAC               ! scratch function of USTAR, NU, and GRAV
      REAL TWOXLM               ! 2 X atmospheric mean free path

      REAL, PARAMETER :: BHAT    = 1.246 ! Constant from Cunningham slip correction
      REAL, PARAMETER :: THREEPI = 3.0 * PI
      REAL, PARAMETER :: TWO3    = 2.0 / 3.0

C Scalar variables for VARIABLE standard deviations.

      REAL    L2SGNU   ! see usage
      REAL    L2SGAT
      REAL    L2SGAC   
      REAL    L2SGCO

      REAL    ENU1             ! Nucleation mode exp( log^2( sigmag )/8 )
      REAL    EAT1             ! Aitken mode exp( log^2( sigmag )/8 )
      REAL    EAC1             ! accumulation mode exp( log^2( sigmag )/8 )
      REAL    ECO1             ! coarse mode exp( log^2( sigmag )/8 )

      REAL    ESNU04           ! Nucleation   " **4
      REAL    ESAT04           ! Aitken       " 
      REAL    ESAC04           ! accumulation "
      REAL    ESCO04           ! coarse       "

      REAL    ESNU08           ! Nucleation   " **8
      REAL    ESAT08           ! Aitken       " 
      REAL    ESAC08           ! accumulation "
      REAL    ESCO08           ! coarse       "

      REAL    ESNU12           ! Nucleation   " **12
      REAL    ESAT12           ! Aitken       " 
      REAL    ESAC12           ! accumulation "     
      REAL    ESCO12           ! coarse       "     

      REAL    ESNU16           ! Nucleation   " **16
      REAL    ESAT16           ! Aitken       " 
      REAL    ESAC16           ! accumulation "
      REAL    ESCO16           ! coarse       "

      REAL    ESNU20           ! Nucleation   " **20
      REAL    ESAT20           ! Aitken       " 
      REAL    ESAC20           ! accumulation "
      REAL    ESCO20           ! coarse       "

      REAL    ESNU28           ! Nucleation   " **28
      REAL    ESAT28           ! Aitken       " 
      REAL    ESAC28           ! accumulation "
      REAL    ESCO28           ! coarse       "

      REAL    ESNU32           ! Nucleation   " **32
      REAL    ESAT32           ! Aitken       " 
      REAL    ESAC32           ! accumulation "
      REAL    ESCO32           ! coarse       "

      REAL    ESNU36           ! Nucleation   " **36
      REAL    ESAT36           ! Aitken       " 
      REAL    ESAC36           ! accumulation "
      REAL    ESCO36           ! coarse       "

      REAL    ESNU48           ! Nucleation   " **48
      REAL    ESAT48           ! Aitken       " 
      REAL    ESAC48           ! accumulation "     
      REAL    ESCO48           ! coarse       "     

      REAL    ESNU64           ! Nucleation   " **64
      REAL    ESAT64           ! Aitken       " 
      REAL    ESAC64           ! accumulation "
      REAL    ESCO64           ! coarse       "

      REAL    ESNU128          ! Nucleation   " **128
      REAL    ESAT128          ! Aitken       " 
      REAL    ESAC128          ! accumulation "
      REAL    ESCO128          ! coarse       "

      REAL    ESNU160          ! Nucleation   " **160
      REAL    ESAT160          ! Aitken       " 
      REAL    ESAC160          ! accumulation "
      REAL    ESCO160          ! coarse       "

      REAL    ESNUM12          ! Nucleation   " **(-12)
      REAL    ESATM12          ! Aitken       " 
      REAL    ESACM12          ! accumulation "        
      REAL    ESCOM12          ! coarse       "        

      REAL    ESNUM16          ! Nucleation   " **(-16)
      REAL    ESATM16          ! Aitken       " 
      REAL    ESACM16          ! accumulation "        
      REAL    ESCOM16          ! coarse       "        

      REAL    ESNUM20          ! Nucleation   " **(-20)
      REAL    ESATM20          ! Aitken       " 
      REAL    ESACM20          ! accumulation "
      REAL    ESCOM20          ! coarse       "

      REAL    ESNUM32          ! Nucleation   " **(-32)
      REAL    ESATM32          ! Aitken       " 
      REAL    ESACM32          ! accumulation "
      REAL    ESCOM32          ! coarse       "

      REAL    EIM              ! Impaction efficiency

C-----------------------------------------------------------------------

      VDEP  = 0.0   ! array assignment
      IF ( MOSAIC ) THEN
         VDEPJ = 0.0   ! array assignment
      END IF

      DO R = 1, SIZE( Met_Data%TEMP2, 2 )
         DO C = 1, SIZE( Met_Data%TEMP2, 1 )

C *** Calculate Knudsen numbers

            TWOXLM = XLM( C,R ) + XLM( C,R )
            KNNUC = TWOXLM / DGNUC( C,R )
            KNATK = TWOXLM / DGATK( C,R )
            KNACC = TWOXLM / DGACC( C,R )
            KNCOR = TWOXLM / DGCOR( C,R )

C *** Calculate functions of variable standard deviation.

            L2SGNU = XXLSGNU( C,R ) * XXLSGNU( C,R )
            L2SGAT = XXLSGAT( C,R ) * XXLSGAT( C,R )
            L2SGAC = XXLSGAC( C,R ) * XXLSGAC( C,R )
            L2SGCO = XXLSGCO( C,R ) * XXLSGCO( C,R )

            ENU1   = EXP( 0.125 * L2SGNU )
            EAT1   = EXP( 0.125 * L2SGAT )
            EAC1   = EXP( 0.125 * L2SGAC )
            ECO1   = EXP( 0.125 * L2SGCO )

            ESNU04  = ENU1 ** 4
            ESAT04  = EAT1 ** 4
            ESAC04  = EAC1 ** 4
            ESCO04  = ECO1 ** 4

            ESNU08  = ESNU04 * ESNU04
            ESAT08  = ESAT04 * ESAT04
            ESAC08  = ESAC04 * ESAC04
            ESCO08  = ESCO04 * ESCO04

            ESNU12  = ESNU04 * ESNU08
            ESAT12  = ESAT04 * ESAT08
            ESAC12  = ESAC04 * ESAC08
            ESCO12  = ESCO04 * ESCO08

            ESNU16  = ESNU08 * ESNU08
            ESAT16  = ESAT08 * ESAT08
            ESAC16  = ESAC08 * ESAC08
            ESCO16  = ESCO08 * ESCO08

            ESNU20  = ESNU16 * ESNU04
            ESAT20  = ESAT16 * ESAT04
            ESAC20  = ESAC16 * ESAC04
            ESCO20  = ESCO16 * ESCO04

            ESNU28  = ESNU20 * ESNU08
            ESAT28  = ESAT20 * ESAT08
            ESAC28  = ESAC20 * ESAC08
            ESCO28  = ESCO20 * ESCO08

            ESNU32  = ESNU16 * ESNU16
            ESAT32  = ESAT16 * ESAT16
            ESAC32  = ESAC16 * ESAC16
            ESCO32  = ESCO16 * ESCO16

            ESNU36  = ESNU16 * ESNU20
            ESAT36  = ESAT16 * ESAT20
            ESAC36  = ESAC16 * ESAC20
            ESCO36  = ESCO16 * ESCO20

            ESNU48  = ESNU36 * ESNU12
            ESAT48  = ESAT36 * ESAT12
            ESAC48  = ESAC36 * ESAC12
            ESCO48  = ESCO36 * ESCO12

            ESNU64  = ESNU32 * ESNU32
            ESAT64  = ESAT32 * ESAT32
            ESAC64  = ESAC32 * ESAC32
            ESCO64  = ESCO32 * ESCO32

            ESNU128 = ESNU64 * ESNU64
            ESAT128 = ESAT64 * ESAT64
            ESAC128 = ESAC64 * ESAC64
            ESCO128 = ESCO64 * ESCO64

            ESNU160 = ESNU128* ESNU32
            ESAT160 = ESAT128* ESAT32
            ESAC160 = ESAC128* ESAC32
            ESCO160 = ESCO128* ESCO32

C *** calculate inverses:

            ESNUM12 = 1.0 / ESNU12
            ESATM12 = 1.0 / ESAT12
            ESACM12 = 1.0 / ESAC12
            ESCOM12 = 1.0 / ESCO12

            ESNUM16 = 1.0 / ESNU16
            ESATM16 = 1.0 / ESAT16
            ESACM16 = 1.0 / ESAC16
            ESCOM16 = 1.0 / ESCO16

            ESNUM20 = 1.0 / ESNU20
            ESATM20 = 1.0 / ESAT20
            ESACM20 = 1.0 / ESAC20
            ESCOM20 = 1.0 / ESCO20

            ESNUM32 = 1.0 / ESNU32
            ESATM32 = 1.0 / ESAT32
            ESACM32 = 1.0 / ESAC32
            ESCOM32 = 1.0 / ESCO32

            DCONST1   = BOLTZMANN * Met_Data%TEMP2( C,R ) / ( THREEPI * AMU( C,R ) )
            DCONST1NU = DCONST1 / DGNUC( C,R )
            DCONST1AT = DCONST1 / DGATK( C,R )
            DCONST1AC = DCONST1 / DGACC( C,R )
            DCONST1CO = DCONST1 / DGCOR( C,R )
            
            DCONST2   = GRAV / ( 18.0 * AMU( C,R ) )
            DCONST3NU = DCONST2 * PDENSNU( C,R ) * DGNUC( C,R ) * DGNUC( C,R )
            DCONST3AT = DCONST2 * PDENSAT( C,R ) * DGATK( C,R ) * DGATK( C,R )
            DCONST3AC = DCONST2 * PDENSAC( C,R ) * DGACC( C,R ) * DGACC( C,R )
            DCONST3CO = DCONST2 * PDENSCO( C,R ) * DGCOR( C,R ) * DGCOR( C,R )

C h-mode
            DCHAT0NU  = DCONST1NU * ( ESNU04  + BHAT * KNNUC * ESNU16 )
            DCHAT2NU  = DCONST1NU * ( ESNUM12 + BHAT * KNNUC * ESNUM16 )
            DCHAT3NU  = DCONST1NU * ( ESNUM20 + BHAT * KNNUC * ESNUM32 )
            VGHAT0NU  = DCONST3NU * ( ESNU16  + BHAT * KNNUC * ESNU04 )
            VGHAT2NU  = DCONST3NU * ( ESNU48  + BHAT * KNNUC * ESNU20 )
            VGHAT3NU  = DCONST3NU * ( ESNU64  + BHAT * KNNUC * ESNU28 )

C i-mode
            DCHAT0AT  = DCONST1AT * ( ESAT04  + BHAT * KNATK * ESAT16 )
            DCHAT2AT  = DCONST1AT * ( ESATM12 + BHAT * KNATK * ESATM16 )
            DCHAT3AT  = DCONST1AT * ( ESATM20 + BHAT * KNATK * ESATM32 )
            VGHAT0AT  = DCONST3AT * ( ESAT16  + BHAT * KNATK * ESAT04 )
            VGHAT2AT  = DCONST3AT * ( ESAT48  + BHAT * KNATK * ESAT20 )
            VGHAT3AT  = DCONST3AT * ( ESAT64  + BHAT * KNATK * ESAT28 )

C j-mode
            DCHAT0AC  = DCONST1AC * ( ESAC04  + BHAT * KNACC * ESAC16 )
            DCHAT2AC  = DCONST1AC * ( ESACM12 + BHAT * KNACC * ESACM16 )
            DCHAT3AC  = DCONST1AC * ( ESACM20 + BHAT * KNACC * ESACM32 )
            VGHAT0AC  = DCONST3AC * ( ESAC16  + BHAT * KNACC * ESAC04 )
            VGHAT2AC  = DCONST3AC * ( ESAC48  + BHAT * KNACC * ESAC20 )
            VGHAT3AC  = DCONST3AC * ( ESAC64  + BHAT * KNACC * ESAC28 )

C coarse mode
            DCHAT0CO  = DCONST1CO * ( ESCO04  + BHAT * KNCOR * ESCO16 )
            DCHAT2CO  = DCONST1CO * ( ESCOM12 + BHAT * KNCOR * ESCOM16 )
            DCHAT3CO  = DCONST1CO * ( ESCOM20 + BHAT * KNCOR * ESCOM32 )
            VGHAT0CO  = DCONST3CO * ( ESCO16  + BHAT * KNCOR * ESCO04 )
            VGHAT2CO  = DCONST3CO * ( ESCO48  + BHAT * KNCOR * ESCO20 )
            VGHAT3CO  = DCONST3CO * ( ESCO64  + BHAT * KNCOR * ESCO28 )

C now calculate the deposition velocities

            NU = AMU( C,R ) / Met_Data%DENS1( C,R )
            USTFAC = Met_Data%USTAR( C,R ) ** 2 / ( GRAV * NU )
            STOKENU = DCONST3NU * USTFAC
            STOKEAT = DCONST3AT * USTFAC
            STOKEAC = DCONST3AC * USTFAC
            STOKECO = DCONST3CO * USTFAC
            UTSCALE = Met_Data%USTAR( C,R )
     &              + 0.24 * Met_Data%WSTAR( C,R )**2
     &              /        Met_Data%USTAR( C,R )

C first do 0th moment for the deposition of number

C  Nucleation mode
            SC0NU = NU / DCHAT0NU
            EIM = STOKENU ** 2 / 400.0 * ESNU64
            EIM = MIN( EIM, 1.0 )
            RD0NU = 1.0 / ( UTSCALE * ( SC0NU ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDNNUC ) = VGHAT0NU
     &             / ( 1.0 - EXP( -VGHAT0NU * ( Met_Data%RA( C,R ) + RD0NU ) ) )
 
C  Aitken mode
            SC0AT = NU / DCHAT0AT
            EIM = STOKEAT ** 2 / 400.0 * ESAT64
            EIM = MIN( EIM, 1.0 )
            RD0AT = 1.0 / ( UTSCALE * ( SC0AT ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDNATK ) = VGHAT0AT
     &             / ( 1.0 - EXP( -VGHAT0AT * ( Met_Data%RA( C,R ) + RD0AT ) ) )

C accumulation mode
            SC0AC = NU / DCHAT0AC
            EIM = STOKEAC ** 2 / 400.0 * ESAC64
            EIM = MIN( EIM, 1.0 )
            RD0AC = 1.0 / ( UTSCALE * ( SC0AC ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDNACC ) = VGHAT0AC
     &             / ( 1.0 - EXP( -VGHAT0AC * ( Met_Data%RA( C,R ) + RD0AC ) ) )

C coarse mode
            SC0CO = NU / DCHAT0CO
            EIM = STOKECO ** 2 / 400.0 * ESCO64
            EIM = MIN( EIM, 1.0 )
            RD0CO = 1.0 / ( UTSCALE * ( SC0CO ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDNCOR ) = VGHAT0CO
     &             / ( 1.0 - EXP( -VGHAT0CO * ( Met_Data%RA( C,R ) + RD0CO ) ) )

C now do 2nd moment for the deposition of surface area
 
C  Nucleation mode
            SC2NU = NU / DCHAT2NU
            EIM = STOKENU ** 2 / 400.0 * ESNU128
            EIM = MIN( EIM, 1.0 )
            RD2NU = 1.0 / ( UTSCALE * ( SC2NU ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDSNUC ) = VGHAT2NU
     &             / ( 1.0 - EXP( -VGHAT2NU * ( Met_Data%RA( C,R ) + RD2NU ) ) )
 
C  Aitken mode
            SC2AT = NU / DCHAT2AT
            EIM = STOKEAT ** 2 / 400.0 * ESAT128
            EIM = MIN( EIM, 1.0 )
            RD2AT = 1.0 / ( UTSCALE * ( SC2AT ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDSATK ) = VGHAT2AT
     &             / ( 1.0 - EXP( -VGHAT2AT * ( Met_Data%RA( C,R ) + RD2AT ) ) )

C accumulation mode
            SC2AC = NU / DCHAT2AC
            EIM = STOKEAC ** 2 / 400.0 * ESAC128
            EIM = MIN( EIM, 1.0 )
            RD2AC = 1.0 / ( UTSCALE * ( SC2AC ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDSACC ) = VGHAT2AC
     &             / ( 1.0 - EXP( -VGHAT2AC * ( Met_Data%RA( C,R ) + RD2AC ) ) )

C coarse mode
            SC2CO = NU / DCHAT2CO
            EIM = STOKECO ** 2 / 400.0 * ESCO128
            EIM = MIN( EIM, 1.0 )
            RD2CO = 1.0 / ( UTSCALE * ( SC2CO ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDSCOR ) = VGHAT2CO
     &             / ( 1.0 - EXP( -VGHAT2CO * ( Met_Data%RA( C,R ) + RD2CO ) ) )

C now do 3rd moment for the deposition of mass

C  Nucleation mode
            SC3NU = NU / DCHAT3NU
            EIM = STOKENU ** 2 / 400.0 * ESNU160
            EIM = MIN( EIM, 1.0 )
            RD3NU = 1.0 / ( UTSCALE * ( SC3NU ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDMNUC ) = VGHAT3NU
     &             / ( 1.0 - EXP( -VGHAT3NU * ( Met_Data%RA( C,R ) + RD3NU ) ) )
 
C  Aitken mode
            SC3AT = NU / DCHAT3AT
            EIM = STOKEAT ** 2 / 400.0 * ESAT160
            EIM = MIN( EIM, 1.0 )
            RD3AT = 1.0 / ( UTSCALE * ( SC3AT ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDMATK ) = VGHAT3AT
     &             / ( 1.0 - EXP( -VGHAT3AT * ( Met_Data%RA( C,R ) + RD3AT ) ) )

C accumulation mode
            SC3AC = NU / DCHAT3AC
            EIM = STOKEAC ** 2 / 400.0 * ESAC160
            EIM = MIN( EIM, 1.0 )
            RD3AC = 1.0 / ( UTSCALE * ( SC3AC ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDMACC ) = VGHAT3AC
     &             / ( 1.0 - EXP( -VGHAT3AC * ( Met_Data%RA( C,R ) + RD3AC ) ) )

C coarse mode
            SC3CO = NU / DCHAT3CO
            EIM = STOKECO ** 2 / 400.0 * ESCO160
            EIM = MIN( EIM, 1.0 )
            RD3CO = 1.0 / ( UTSCALE * ( SC3CO ** ( -TWO3 ) + EIM ) )

            VDEP( C,R,VDMCOR ) = VGHAT3CO
     &             / ( 1.0 - EXP( -VGHAT3CO * ( Met_Data%RA( C,R ) + RD3CO ) ) )


C Do mosaic calculations - essentially a repeat of the above, using the mosaic vars
            IF ( MOSAIC ) THEN

               DO J = 1, N_LUFRAC
                  IF ( Mosaic_Data%RA( C,R,J ) .GT. 0.0 .And. Grid_Data%LUFRAC( C,R,J ) .GT. 0.0 ) THEN

C now calculate the deposition velocities

                     NU = AMU( C,R ) / Met_Data%DENS1( C,R )
                     USTFAC = Mosaic_Data%USTAR( C,R,J )**2 / ( GRAV * NU )
                     STOKENU = DCONST3NU * USTFAC
                     STOKEAT = DCONST3AT * USTFAC
                     STOKEAC = DCONST3AC * USTFAC
                     STOKECO = DCONST3CO * USTFAC
                     UTSCALE = Mosaic_Data%USTAR( C,R,J )
     &                       + 0.24 * Met_Data%WSTAR( C,R )**2
     &                       /        Mosaic_Data%USTAR( C,R,J )

C first do 0th moment for the deposition of number

C  Nucleation mode

                     SC0NU = NU / DCHAT0NU
                     EIM = STOKENU ** 2 / 400.0 * ESNU64
                     EIM = MIN( EIM, 1.0 )
                     RD0NU = 1.0 / ( UTSCALE * ( SC0NU ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDNNUC ) = VGHAT0NU
     &                       / ( 1.0 - EXP(-VGHAT0NU * ( Mosaic_Data%RA( C,R,J ) + RD0NU ) ) )
 
C  Aitken mode

                     SC0AT = NU / DCHAT0AT
                     EIM = STOKEAT ** 2 / 400.0 * ESAT64
                     EIM = MIN( EIM, 1.0 )
                     RD0AT = 1.0 / ( UTSCALE * ( SC0AT ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDNATK ) = VGHAT0AT
     &                       / ( 1.0 - EXP(-VGHAT0AT * ( Mosaic_Data%RA( C,R,J ) + RD0AT ) ) )

C accumulation mode

                     SC0AC = NU / DCHAT0AC
                     EIM = STOKEAC ** 2 / 400.0 * ESAC64
                     EIM = MIN( EIM, 1.0 )
                     RD0AC = 1.0 / ( UTSCALE * ( SC0AC ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDNACC ) = VGHAT0A
     &                       / ( 1.0 - EXP(-VGHAT0A * ( Mosaic_Data%RA( C,R,J ) + RD0A ) ) )

C coarse mode

                     SC0CO = NU / DCHAT0CO
                     EIM = STOKECO ** 2 / 400.0 * ESCO64
                     EIM = MIN( EIM, 1.0 )
                     RD0CO = 1.0 / ( UTSCALE * ( SC0CO ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDNCOR ) = VGHAT0CO
     &                       / ( 1.0 - EXP(-VGHAT0CO * ( Mosaic_Data%RA( C,R,J ) + RD0CO ) ) )

C now do 2nd moment for the deposition of surface area

C  Nucleation mode

                     SC2NU = NU / DCHAT2NU
                     EIM = STOKENU ** 2 / 400.0 * ESNU128
                     EIM = MIN( EIM, 1.0 )
                     RD2NU = 1.0 / ( UTSCALE * ( SC2NU ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDSNUC ) = VGHAT2NU
     &                     / ( 1.0 - EXP(-VGHAT2NU * ( Mosaic_Data%RA( C,R,J ) + RD2NU ) ) )

C  Aitken mode

                     SC2AT = NU / DCHAT2AT
                     EIM = STOKEAT ** 2 / 400.0 * ESAT128
                     EIM = MIN( EIM, 1.0 )
                     RD2AT = 1.0 / ( UTSCALE * ( SC2AT ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDSATK ) = VGHAT2AT
     &                     / ( 1.0 - EXP(-VGHAT2AT * ( Mosaic_Data%RA( C,R,J ) + RD2AT ) ) )

C accumulation mode

                     SC2AC = NU / DCHAT2AC
                     EIM = STOKEAC ** 2 / 400.0 * ESAC128
                     EIM = MIN( EIM, 1.0 )
                     RD2AC = 1.0 / ( UTSCALE * ( SC2AC ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDSACC ) = VGHAT2AC
     &                       / ( 1.0 - EXP(-VGHAT2AC * ( Mosaic_Data%RA( C,R,J ) + RD2AC ) ) )

C coarse mode

                     SC2CO = NU / DCHAT2CO
                     EIM = STOKECO ** 2 / 400.0 * ESCO128
                     EIM = MIN( EIM, 1.0 )
                     RD2CO = 1.0 / ( UTSCALE * ( SC2CO ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDSCOR ) = VGHAT2CO
     &                       / ( 1.0 - EXP(-VGHAT2CO * ( Mosaic_Data%RA( C,R,J ) + RD2CO ) ) )

C now do 3rd moment for the deposition of mass

C  Nucleation mode

                     SC3NU = NU / DCHAT3NU
                     EIM = STOKENU ** 2 / 400.0 * ESNU160
                     EIM = MIN( EIM, 1.0 )
                     RD3NU = 1.0 / ( UTSCALE * ( SC3NU ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDMATK ) = VGHAT3NU
     &                       / ( 1.0 - EXP(-VGHAT3NU * ( Mosaic_Data%RA( C,R,J ) + RD3NU ) ) )
 
C  Aitken mode

                     SC3AT = NU / DCHAT3AT
                     EIM = STOKEAT ** 2 / 400.0 * ESAT160
                     EIM = MIN( EIM, 1.0 )
                     RD3AT = 1.0 / ( UTSCALE * ( SC3AT ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDMATK ) = VGHAT3AT
     &                       / ( 1.0 - EXP(-VGHAT3AT * ( Mosaic_Data%RA( C,R,J ) + RD3AT ) ) )

C accumulation mode

                     SC3AC = NU / DCHAT3AC
                     EIM = STOKEAC ** 2 / 400.0 * ESAC160
                     EIM = MIN( EIM, 1.0 )
                     RD3AC = 1.0 / ( UTSCALE * ( SC3AC ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDMACC ) = VGHAT3AC
     &                       / ( 1.0 - EXP(-VGHAT3AC * ( Mosaic_Data%RA( C,R,J ) + RD3AC ) ) )

C coarse mode

                     SC3CO = NU / DCHAT3CO
                     EIM = STOKECO ** 2 / 400.0 * ESCO160
                     EIM = MIN( EIM, 1.0 )
                     RD3CO = 1.0 / ( UTSCALE * ( SC3CO ** ( -TWO3 ) + EIM ) )

                     VDEPJ( J,C,R,VDMCOR ) = VGHAT3CO
     &                     / ( 1.0 - EXP(-VGHAT3CO * ( Mosaic_Data%RA( C,R,J ) + RD3CO ) ) )

                  END IF ! RA > 0

               END DO ! n_lufrac

            END IF ! mosaic
         END DO ! end loop on C
      END DO ! end loop on R

      RETURN
      END SUBROUTINE GETDEP_V
