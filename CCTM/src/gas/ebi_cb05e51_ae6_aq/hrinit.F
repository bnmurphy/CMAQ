
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
C $Header$

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE HRINIT
C***********************************************************************
C
C  FUNCTION: To initialize species tolerances, arrays, and indices
C
C  PRECONDITIONS: For the CB05E51 family of mechanisms
C
C  RETURN VALUES: None
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Created by EBI solver program, Aug 26, 2017
C   21 Jun 10 J.Young:   convert for Namelist redesign
C   18 Jul 14 B.Hutzell: revised to use the RXNS_DATA and RXNS_FUNCTION
C                        modules instead of include files
C***********************************************************************

      USE UTILIO_DEFN
      USE HRDATA
      USE RXNS_DATA
      USE RXNS_FUNCTION

      IMPLICIT NONE

C.....INCLUDES:
      Include SUBST_CONST          ! common constants

C.....ARGUMENTS: NONE

C.....PARAMETERS: NONE


C.....SAVED VARIABLES:
      CHARACTER( 16 ), SAVE  ::  PNAME = 'HRINIT'   ! Program name


C.....LOCAL VARIABLES:
      CHARACTER( 132 ) ::  MSG       ! Log message

      INTEGER N               ! Loop index
      INTEGER IOS             ! Allocate error status

C***********************************************************************

      N_SPEC = NUMB_MECH_SPC
      N_RXNS = NRXNS
      N_SPCL_RKS =   0
      ALLOCATE( RKI( NRXNS ),
     &          RKI_SAV( MAX_CELLS_EBI, NRXNS ),
     &          FORWARD_CONV( N_SPEC),
     &          REVERSE_CONV( N_SPEC),
     &          RXRAT( NRXNS ),
     &          RTOL( N_SPEC),
     &          YC(   N_SPEC),
     &          YC0(  N_SPEC),
     &          YCP(  N_SPEC),
     &          PROD( N_SPEC),
     &          LOSS( N_SPEC),
     &          STAT = IOS )

       IF ( IOS .NE. 0 ) THEN
          MSG = 'Error allocating RKI, RXRAT, RTOL, YC, YC0, YCP, FORWARD_CONV, '
     &       // 'REVERSE_CONV, PROD, LOSS, or PNEG'
          CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
       END IF

       ALLOCATE( ATMPRES( MAX_CELLS_EBI ),
     &           H2O    ( MAX_CELLS_EBI ),
     &           TEMP   ( MAX_CELLS_EBI ),
     &           HET    ( MAX_CELLS_EBI, NHETERO ),
     &           RJIN   ( MAX_CELLS_EBI, NPHOTAB ),
     &           LAND   ( MAX_CELLS_EBI ) ,
     &           STAT = IOS )

       IF ( IOS .NE. 0 ) THEN
          MSG = 'Error allocating ATMPRES, H2O, TEMP, HET, RJIN, LAND '
          CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
       END IF

       IF( NSPECIAL .GT. 0 )THEN
          ALLOCATE( SYC( MAX_CELLS_EBI, N_SPEC ), STAT = IOS )
          IF ( IOS .NE. 0 ) THEN
               MSG = 'Error allocating SYC'
               CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
          END IF
          ALLOCATE( IOLD2NEW( N_SPEC, 1 ), STAT = IOS )
          IF ( IOS .NE. 0 ) THEN
               MSG = 'Error allocating IOLD2NEW'
               CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
          END IF
          DO N = 1, NUMB_MECH_SPC
             IOLD2NEW( N, 1 ) = N
          END DO
       END IF

       IF( .NOT. MAP_CHEMISTRY_SPECIES() )THEN
          MSG = 'Detected above error(s) when mapping Chemistry species from CGRID species'
          CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT2 )
       END IF

       DO N = 1, NUMB_MECH_SPC
          FORWARD_CONV( N ) = 1.0E-3 * MWAIR / SPECIES_MOLWT( N )
          REVERSE_CONV( N ) = 1.0E+3 / MWAIR * SPECIES_MOLWT( N )
       END DO

! initialize all reaction rates to zero so unused array elements are defined

       RXRAT = 0.0D0

! initialize all production and loss rates to zero so unused array elements are defined

       PROD  = 0.0D0
       LOSS  = 0.0D0

       NING1 = 4
       NING2 = 4



c..Allocate and define ebi species

      N_EBISP  = 145
      ALLOCATE( EBISP( N_EBISP ), STAT = IOS ) 
      IF ( IOS .NE. 0 ) THEN
         MSG = 'Error allocating EBISP' 
         CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
      END IF

      EBISP(   1 ) = HNO3
      EBISP(   2 ) = H2O2
      EBISP(   3 ) = XO2
      EBISP(   4 ) = XO2N
      EBISP(   5 ) = NTROH
      EBISP(   6 ) = NTRALK
      EBISP(   7 ) = ROOH
      EBISP(   8 ) = ALD2
      EBISP(   9 ) = ALDX
      EBISP(  10 ) = ISOPX
      EBISP(  11 ) = IEPOX
      EBISP(  12 ) = ISOPO2
      EBISP(  13 ) = IOLE
      EBISP(  14 ) = IEPXO2
      EBISP(  15 ) = MGLY
      EBISP(  16 ) = FORM
      EBISP(  17 ) = FACD
      EBISP(  18 ) = CO
      EBISP(  19 ) = PAR
      EBISP(  20 ) = MEO2
      EBISP(  21 ) = AACD
      EBISP(  22 ) = MEPX
      EBISP(  23 ) = MEOH
      EBISP(  24 ) = HCO3
      EBISP(  25 ) = PACD
      EBISP(  26 ) = CXO3
      EBISP(  27 ) = PANX
      EBISP(  28 ) = ROR
      EBISP(  29 ) = OLE
      EBISP(  30 ) = ETH
      EBISP(  31 ) = TOL
      EBISP(  32 ) = CRES
      EBISP(  33 ) = TO2
      EBISP(  34 ) = TOLRO2
      EBISP(  35 ) = OPEN
      EBISP(  36 ) = CRO
      EBISP(  37 ) = CAT1
      EBISP(  38 ) = CRON
      EBISP(  39 ) = CRNO
      EBISP(  40 ) = CRN2
      EBISP(  41 ) = CRPX
      EBISP(  42 ) = OPO3
      EBISP(  43 ) = CAO2
      EBISP(  44 ) = OPAN
      EBISP(  45 ) = XYLMN
      EBISP(  46 ) = XYLRO2
      EBISP(  47 ) = NAPH
      EBISP(  48 ) = PAHRO2
      EBISP(  49 ) = ISOP
      EBISP(  50 ) = ISPD
      EBISP(  51 ) = ISOPRXN
      EBISP(  52 ) = IHOM
      EBISP(  53 ) = NTRM
      EBISP(  54 ) = MACO3
      EBISP(  55 ) = NTRI
      EBISP(  56 ) = TERP
      EBISP(  57 ) = TRPRXN
      EBISP(  58 ) = XO2T
      EBISP(  59 ) = THOM
      EBISP(  60 ) = SO2
      EBISP(  61 ) = SULF
      EBISP(  62 ) = SULRXN
      EBISP(  63 ) = ETOH
      EBISP(  64 ) = ETHA
      EBISP(  65 ) = CL2
      EBISP(  66 ) = CL
      EBISP(  67 ) = HOCL
      EBISP(  68 ) = CLO
      EBISP(  69 ) = FMCL
      EBISP(  70 ) = HCL
      EBISP(  71 ) = CLNO2
      EBISP(  72 ) = TOLNRXN
      EBISP(  73 ) = TOLHRXN
      EBISP(  74 ) = XYLNRXN
      EBISP(  75 ) = XYLHRXN
      EBISP(  76 ) = BENZENE
      EBISP(  77 ) = BENZRO2
      EBISP(  78 ) = BNZNRXN
      EBISP(  79 ) = BNZHRXN
      EBISP(  80 ) = SESQ
      EBISP(  81 ) = SESQRXN
      EBISP(  82 ) = PAHNRXN
      EBISP(  83 ) = PAHHRXN
      EBISP(  84 ) = SOAALK
      EBISP(  85 ) = ALKRXN
      EBISP(  86 ) = MAPAN
      EBISP(  87 ) = NALKO2
      EBISP(  88 ) = NTRCN
      EBISP(  89 ) = NTRCNOH
      EBISP(  90 ) = NTRPX
      EBISP(  91 ) = NOHO2
      EBISP(  92 ) = NCNO2
      EBISP(  93 ) = NCNOHO2
      EBISP(  94 ) = NTRMO2
      EBISP(  95 ) = NTRIO2
      EBISP(  96 ) = H2NO3PIJ
      EBISP(  97 ) = H2NO3PK
      EBISP(  98 ) = ACLI
      EBISP(  99 ) = ACLJ
      EBISP( 100 ) = ACLK
      EBISP( 101 ) = AISO3J
      EBISP( 102 ) = AXYL1J
      EBISP( 103 ) = AOLGAJ
      EBISP( 104 ) = AXYL2J
      EBISP( 105 ) = ATOL1J
      EBISP( 106 ) = ATOL2J
      EBISP( 107 ) = ABNZ1J
      EBISP( 108 ) = ABNZ2J
      EBISP( 109 ) = ATRP1J
      EBISP( 110 ) = AOLGBJ
      EBISP( 111 ) = ATRP2J
      EBISP( 112 ) = AISO1J
      EBISP( 113 ) = AISO2J
      EBISP( 114 ) = ASQTJ
      EBISP( 115 ) = APAH1J
      EBISP( 116 ) = APAH2J
      EBISP( 117 ) = AALK1J
      EBISP( 118 ) = AALK2J
      EBISP( 119 ) = DMA
      EBISP( 120 ) = PCVOC
      EBISP( 121 ) = PCSOARXN
      EBISP( 122 ) = VLVPO1
      EBISP( 123 ) = VSVPO1
      EBISP( 124 ) = VSVPO2
      EBISP( 125 ) = VSVPO3
      EBISP( 126 ) = VIVPO1
      EBISP( 127 ) = VLVOO1
      EBISP( 128 ) = VLVOO2
      EBISP( 129 ) = VSVOO2
      EBISP( 130 ) = VSVOO3
      EBISP( 131 ) = VSVOO1
      EBISP( 132 ) = FORM_PRIMARY
      EBISP( 133 ) = ALD2_PRIMARY
      EBISP( 134 ) = BUTADIENE13
      EBISP( 135 ) = ACROLEIN
      EBISP( 136 ) = ACRO_PRIMARY
      EBISP( 137 ) = TOLU
      EBISP( 138 ) = MXYL
      EBISP( 139 ) = OXYL
      EBISP( 140 ) = PXYL
      EBISP( 141 ) = APIN
      EBISP( 142 ) = BPIN
      EBISP( 143 ) = HG
      EBISP( 144 ) = HGIIAER
      EBISP( 145 ) = HGIIGAS


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set species tolerances
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      RTOL( NO2          ) = 5.0D-04
      RTOL( NO           ) = 5.0D-04
      RTOL( O            ) = 1.0D+00
      RTOL( O3           ) = 5.0D-04
      RTOL( NO3          ) = 5.0D-04
      RTOL( O1D          ) = 1.0D+00
      RTOL( OH           ) = 5.0D-04
      RTOL( HO2          ) = 5.0D-04
      RTOL( N2O5         ) = 5.0D-04
      RTOL( HNO3         ) = 5.0D-04
      RTOL( HONO         ) = 5.0D-04
      RTOL( PNA          ) = 5.0D-04
      RTOL( H2O2         ) = 5.0D-04
      RTOL( XO2          ) = 5.0D-04
      RTOL( XO2N         ) = 5.0D-04
      RTOL( NTROH        ) = 5.0D-04
      RTOL( NTRALK       ) = 5.0D-04
      RTOL( ROOH         ) = 5.0D-04
      RTOL( ALD2         ) = 5.0D-04
      RTOL( ALDX         ) = 5.0D-04
      RTOL( ISOPX        ) = 5.0D-04
      RTOL( IEPOX        ) = 5.0D-04
      RTOL( ISOPO2       ) = 5.0D-04
      RTOL( IOLE         ) = 5.0D-04
      RTOL( IEPXO2       ) = 5.0D-04
      RTOL( MGLY         ) = 5.0D-04
      RTOL( FORM         ) = 5.0D-04
      RTOL( FACD         ) = 5.0D-04
      RTOL( CO           ) = 5.0D-04
      RTOL( PAR          ) = 5.0D-04
      RTOL( C2O3         ) = 5.0D-04
      RTOL( MEO2         ) = 5.0D-04
      RTOL( AACD         ) = 5.0D-04
      RTOL( MEPX         ) = 5.0D-04
      RTOL( MEOH         ) = 5.0D-04
      RTOL( HCO3         ) = 1.0D+00
      RTOL( PAN          ) = 5.0D-04
      RTOL( PACD         ) = 5.0D-04
      RTOL( CXO3         ) = 5.0D-04
      RTOL( PANX         ) = 5.0D-04
      RTOL( ROR          ) = 5.0D-04
      RTOL( OLE          ) = 5.0D-04
      RTOL( ETH          ) = 5.0D-04
      RTOL( TOL          ) = 5.0D-04
      RTOL( CRES         ) = 5.0D-04
      RTOL( TO2          ) = 5.0D-04
      RTOL( TOLRO2       ) = 5.0D-04
      RTOL( OPEN         ) = 5.0D-04
      RTOL( CRO          ) = 5.0D-04
      RTOL( CAT1         ) = 5.0D-04
      RTOL( CRON         ) = 5.0D-04
      RTOL( CRNO         ) = 5.0D-04
      RTOL( CRN2         ) = 5.0D-04
      RTOL( CRPX         ) = 5.0D-04
      RTOL( OPO3         ) = 5.0D-04
      RTOL( CAO2         ) = 5.0D-04
      RTOL( OPAN         ) = 5.0D-04
      RTOL( XYLMN        ) = 5.0D-04
      RTOL( XYLRO2       ) = 5.0D-04
      RTOL( NAPH         ) = 5.0D-04
      RTOL( PAHRO2       ) = 1.0D+00
      RTOL( ISOP         ) = 5.0D-04
      RTOL( ISPD         ) = 5.0D-04
      RTOL( ISOPRXN      ) = 1.0D+00
      RTOL( IHOM         ) = 5.0D-04
      RTOL( NTRM         ) = 5.0D-04
      RTOL( MACO3        ) = 5.0D-04
      RTOL( NTRI         ) = 5.0D-04
      RTOL( TERP         ) = 5.0D-04
      RTOL( TRPRXN       ) = 1.0D+00
      RTOL( XO2T         ) = 5.0D-04
      RTOL( THOM         ) = 5.0D-04
      RTOL( SO2          ) = 5.0D-04
      RTOL( SULF         ) = 1.0D+00
      RTOL( SULRXN       ) = 1.0D+00
      RTOL( ETOH         ) = 5.0D-04
      RTOL( ETHA         ) = 5.0D-04
      RTOL( CL2          ) = 5.0D-04
      RTOL( CL           ) = 5.0D-04
      RTOL( HOCL         ) = 5.0D-04
      RTOL( CLO          ) = 5.0D-04
      RTOL( FMCL         ) = 5.0D-04
      RTOL( HCL          ) = 5.0D-04
      RTOL( CLNO2        ) = 5.0D-04
      RTOL( TOLNRXN      ) = 1.0D+00
      RTOL( TOLHRXN      ) = 1.0D+00
      RTOL( XYLNRXN      ) = 1.0D+00
      RTOL( XYLHRXN      ) = 1.0D+00
      RTOL( BENZENE      ) = 5.0D-04
      RTOL( BENZRO2      ) = 5.0D-04
      RTOL( BNZNRXN      ) = 1.0D+00
      RTOL( BNZHRXN      ) = 1.0D+00
      RTOL( SESQ         ) = 5.0D-04
      RTOL( SESQRXN      ) = 1.0D+00
      RTOL( PAHNRXN      ) = 1.0D+00
      RTOL( PAHHRXN      ) = 1.0D+00
      RTOL( SOAALK       ) = 5.0D-04
      RTOL( ALKRXN       ) = 1.0D+00
      RTOL( MAPAN        ) = 5.0D-04
      RTOL( NALKO2       ) = 5.0D-04
      RTOL( NTRCN        ) = 5.0D-04
      RTOL( NTRCNOH      ) = 5.0D-04
      RTOL( NTRPX        ) = 5.0D-04
      RTOL( NOHO2        ) = 5.0D-04
      RTOL( NCNO2        ) = 5.0D-04
      RTOL( NCNOHO2      ) = 5.0D-04
      RTOL( NTRMO2       ) = 5.0D-04
      RTOL( NTRIO2       ) = 5.0D-04
      RTOL( H2NO3PIJ     ) = 5.0D-04
      RTOL( H2NO3PK      ) = 5.0D-04
      RTOL( ACLI         ) = 5.0D-04
      RTOL( ACLJ         ) = 5.0D-04
      RTOL( ACLK         ) = 5.0D-04
      RTOL( AISO3J       ) = 5.0D-04
      RTOL( AXYL1J       ) = 5.0D-04
      RTOL( AOLGAJ       ) = 5.0D-04
      RTOL( AXYL2J       ) = 5.0D-04
      RTOL( ATOL1J       ) = 5.0D-04
      RTOL( ATOL2J       ) = 5.0D-04
      RTOL( ABNZ1J       ) = 5.0D-04
      RTOL( ABNZ2J       ) = 5.0D-04
      RTOL( ATRP1J       ) = 5.0D-04
      RTOL( AOLGBJ       ) = 5.0D-04
      RTOL( ATRP2J       ) = 5.0D-04
      RTOL( AISO1J       ) = 5.0D-04
      RTOL( AISO2J       ) = 5.0D-04
      RTOL( ASQTJ        ) = 5.0D-04
      RTOL( APAH1J       ) = 5.0D-04
      RTOL( APAH2J       ) = 5.0D-04
      RTOL( AALK1J       ) = 5.0D-04
      RTOL( AALK2J       ) = 5.0D-04
      RTOL( DMA          ) = 5.0D-04
      RTOL( PCVOC        ) = 5.0D-04
      RTOL( PCSOARXN     ) = 5.0D-04
      RTOL( VLVPO1       ) = 5.0D-04
      RTOL( VSVPO1       ) = 5.0D-04
      RTOL( VSVPO2       ) = 5.0D-04
      RTOL( VSVPO3       ) = 5.0D-04
      RTOL( VIVPO1       ) = 5.0D-04
      RTOL( VLVOO1       ) = 5.0D-04
      RTOL( VLVOO2       ) = 5.0D-04
      RTOL( VSVOO2       ) = 5.0D-04
      RTOL( VSVOO3       ) = 5.0D-04
      RTOL( VSVOO1       ) = 5.0D-04
      RTOL( FORM_PRIMARY ) = 5.0D-04
      RTOL( ALD2_PRIMARY ) = 5.0D-04
      RTOL( BUTADIENE13  ) = 5.0D-04
      RTOL( ACROLEIN     ) = 5.0D-04
      RTOL( ACRO_PRIMARY ) = 5.0D-04
      RTOL( TOLU         ) = 5.0D-04
      RTOL( MXYL         ) = 5.0D-04
      RTOL( OXYL         ) = 5.0D-04
      RTOL( PXYL         ) = 5.0D-04
      RTOL( APIN         ) = 5.0D-04
      RTOL( BPIN         ) = 5.0D-04
      RTOL( HG           ) = 5.0D-04
      RTOL( HGIIAER      ) = 5.0D-04
      RTOL( HGIIGAS      ) = 5.0D-04


      RETURN

      END
