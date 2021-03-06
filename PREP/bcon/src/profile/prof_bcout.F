
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
C $Header: /project/yoj/arc/BCON/src/profile/prof_bcout.F,v 1.2 2011/10/21 16:52:36 yoj Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%


      SUBROUTINE BCOUT( LOGUNIT, SDATE, STIME, NHRS, NSPCS_IN, NSPCS_OUT,
     &                  SPMAP, PROF_FL_NAME, PROF_SP_NAME, SPNAME_OUT )

C***********************************************************************
 
C  Function: Reads the input profile concentration file and opens and
C            writes the output BC file
              
C  Preconditions: None
  
C  Key Subroutines/Functions Called: OPN_BC_FILE
C                                    PROF_VINTERP 
 
C  Revision History:
C    Prototype created by Jerry Gipson, January, 1998
C    Modified April 1998 by JG to remove ppb option for input profile file
C    02/25/02 Steve Howard (Jeff Young) - dynamic allocation
C    01/05/05 J.Young: vert dyn alloc - Use VGRD_DEFN
C    13 Jul 11 J.Young: Replaced I/O API include files with M3UTILIO and
C                       Namelist for species definitions
C    23 May 12 J.Young: Replaced BC_PARMS include file with an F90 module
 
C***********************************************************************

      USE HGRD_DEFN   ! Module to store and load the horizontal grid variables
      USE VGRD_DEFN   ! vertical layer specifications
      USE M3UTILIO    ! IOAPI module
      USE CGRID_SPCS
      USE BC_PARMS    ! BCON parameters

      IMPLICIT NONE     

C Include Files: None

C Arguments: 
      INTEGER LOGUNIT          ! Unit number for output log
      INTEGER SDATE            ! Date for IC Output
      INTEGER STIME            ! Time for IC output
      INTEGER NHRS             ! No. of hours for BC output
      INTEGER NSPCS_IN         ! Total No. of species in input conc file(s)
      INTEGER NSPCS_OUT        ! Number of IC species on output file(s)
      INTEGER :: SPMAP( * )    ! Map to input file species from IC species
      CHARACTER( * ) :: PROF_FL_NAME       ! Input profile file name
      CHARACTER( * ) :: PROF_SP_NAME( * )  ! Profile species names
      CHARACTER( * ) :: SPNAME_OUT( * )    ! Ouput file species names

C Parameters: None

C External Functions: None
 
C Local Variables:
      CHARACTER( 80 ) :: MSG      ! Log message
      CHARACTER( 16 ) :: PNAME = 'BCOUT'    ! Procedure Name
      CHARACTER( 16 ) :: PROF_SP  ! Species name on profile file
      CHARACTER( 16 ) :: VNAME    ! Species name on CTM conc file

      CHARACTER( 16 ) :: BC_FNAME( MXCTMS )  ! Logical names of IC Output file(s) 

      INTEGER C             ! Column loop indices
      INTEGER EDG           ! Boundary edge loop index
      INTEGER ES, EE        ! East boundary loop indices (East start, East end)
      INTEGER FLN           ! IC output file number
      INTEGER IND           ! Array indices for species
      INTEGER L             ! Layer loop index
      INTEGER N, SPC        ! Loop indices for species
      INTEGER NLEVS_IN      ! No. of layers in input conc file
      INTEGER NS, NE        ! North boundary loop indices (North start, North end)
      INTEGER PFILE         ! Unit number of profile file   
      INTEGER SS, SE        ! South boundary loop indices (South start, South end)
      INTEGER WS, WE        ! West boundary loop indices (West start, West end)

      LOGICAL LNEG          ! Flag for negative concentrations

      REAL    :: VGLVS_IN( MXLAYS3 + 1 )            ! Input vertical levels
      REAL    :: COUT( NBNDY,NLAYS )                ! Output BC conc
      REAL    :: INPROF( MXLAYS3,4,MX_INFL_SP )     ! Input conc profiles  
      REAL    :: VIPROF( NLAYS,4,MX_INFL_SP )       ! Vertically interp profiles
                            
C***********************************************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Compute loop paramters for boundary conc array
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( NTHIK .EQ. 1 ) THEN

         SS = 1
         SE = SS + NCOLS     ! = ncols + 1

         ES = SE + 1         ! = ncols + 2
         EE = ES + NROWS     ! = ncols + 2 + nrows

         NS = EE + 1         ! = ncols + 2 + nrows + 1
         NE = NS + NCOLS     ! = 2*ncols + nrows + 3

         WS = NE + 1         ! = 2*ncols + nrows + 4
         WE = WS + NROWS     ! = 2*(ncols + nrows) + 4  = nbndy

      ELSE

         MSG = 'Only NTHIK = 1 currently allowed'
         CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT2 )

      END IF 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Read the input profile file data  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PFILE = JUNIT()

      OPEN( UNIT = PFILE, FILE = PROF_FL_NAME, ERR = 999 )

C Consume the 3-line header
      DO N = 1, 3
         READ( PFILE, 94000 ) 
      END DO

C Get file info (no. of species already retrieved)
      READ( PFILE, * ) NLEVS_IN, SPC, ( VGLVS_IN( L ), L = 1, NLEVS_IN + 1 )
!     READ( PFILE, * ) SPC, NLEVS_IN, ( VGLVS_IN( L ), L = 1, NLEVS_IN + 1 )

C Consume a date and time line
      READ( PFILE, * ) 

C Get file data
      DO EDG = 1, 4
         READ( PFILE, * )
         DO SPC = 1, NSPCS_IN
            READ( PFILE, * ) PROF_SP, ( INPROF( L,EDG,SPC ), L = 1, NLEVS_IN )
         END DO
      END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Possible vertical interpolation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL PROF_VINTERP( LOGUNIT, NSPCS_IN, NLEVS_IN, VGLVS_IN, INPROF,
     &                   VIPROF )      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Open the Models3 BC output file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL OPN_BC_FILE( LOGUNIT, SDATE, STIME, NHRS, NSPCS_OUT, 
     &                  SPNAME_OUT, BC_FNAME, 1 )  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Write the output BC concentrations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      LNEG = .FALSE.

      DO N = 1, NSPCS_OUT

         FLN = ( N - 1 ) / MXVARS3 + 1
         VNAME = SPNAME_OUT( N )
         IND = SPMAP( N )

         DO L = 1, NLAYS

            DO C = 1, NBNDY
               COUT( C,L ) = 0.0
            END DO                   

c..North
            DO C = NS, NE
               COUT( C,L ) = COUT( C,L ) + VIPROF( L,1, IND )
               IF ( COUT( C,L ) .LT. 0.0 ) LNEG = .TRUE.
            END DO

c..East
            DO C = ES, EE
               COUT( C,L ) = COUT( C,L ) + VIPROF( L,2,IND )
               IF ( COUT( C,L ) .LT. 0.0 ) LNEG = .TRUE.
            END DO

c..South
            DO C = SS, SE
               COUT( C,L ) = COUT( C,L ) + VIPROF( L,3,IND )
               IF ( COUT( C,L ) .LT. 0.0 ) LNEG = .TRUE.
            END DO

c..West
            DO C = WS, WE
               COUT( C,L ) = COUT( C,L ) + VIPROF( L,4,IND )
               IF ( COUT( C,L ) .LT. 0.0 ) LNEG = .TRUE.
            END DO

         END DO

         IF ( .NOT. WRITE3( BC_FNAME( FLN ), VNAME, SDATE, STIME,
     &                      COUT ) ) THEN
            MSG =  'Could not WRITE species ' //  VNAME // 
     &             'to file ' // BC_FNAME( FLN ) 
            CALL M3EXIT( PNAME, SDATE, STIME, MSG, XSTAT2 )
         END IF

      END DO

      IF ( LNEG ) THEN
         MSG = 'Negative ICs output'
         CALL M3EXIT( PNAME, 0, 0, ' ', XSTAT2 ) 
      END IF

      RETURN

999   CONTINUE 

      MSG = 'Could not open file ' // PROF_FL_NAME
      CALL M3EXIT( PNAME, 0, 0, ' ', XSTAT2 ) 

C************************* FORMAT STATEMENTS ***************************

94000 FORMAT( A )

      END
