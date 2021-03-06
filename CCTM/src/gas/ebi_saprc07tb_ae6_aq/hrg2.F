
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

      SUBROUTINE HRG2( DTC )
C**********************************************************************
C
C  FUNCTION: To solve for the concentration of HO, HO2, HONO, and
C            HNO4 alebraically.
C
C  PRECONDITIONS: For the SAPRC07TB family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Created by EBI solver program, May 26, 2017
C
C   18 Jul 14 B.Hutzell: revised to use real(8) variables
C**********************************************************************
      USE HRDATA

      IMPLICIT NONE


C..INCLUDES: None


C..ARGUMENTS:
      REAL( 8 ), INTENT( IN ) :: DTC                      ! Time step

C..PARAMETERS: None


C..EXTERNAL FUNCTIONS: NONE


C..SAVED LOCAL VARIABLES:
!     CHARACTER( 16 ), SAVE :: PNAME = 'HRG2'             ! Program name


C..SCRATCH LOCAL VARIABLES:
      REAL( 8 ) ::    O1D_S                 ! sum of O1D loss frequencies
      REAL( 8 ) ::    OH_S                  ! stoich coeff for OH from O1D+H2O or H2
      REAL( 8 ) ::    HO2_S                 ! stoich coeff for HO2 from O1D+H2 rxn
      REAL( 8 ) ::    XOH_S                 ! stoich coeff for OH & HONO from NO2EX loss rxn
      REAL( 8 ) ::    R4_19                 ! production of OH from HONO
      REAL( 8 ) ::    R19_4                 ! production of HONO from OH
      REAL( 8 ) ::    R4_5                  ! production of OH from HO2
      REAL( 8 ) ::    R5_4                  ! production of HO2 from OH
      REAL( 8 ) ::    R5_21                 ! production of HO2 from HNO4
      REAL( 8 ) ::    R21_5                 ! production of HNO4 from HO2
      REAL( 8 ) ::    P4, P5, P19, P21      ! species production form other paths
      REAL( 8 ) ::    L4, L5, L19, L21      ! species loss terms

      REAL( 8 ) ::    A, B, C               ! coeffs of quadratic eq. for HO2
      REAL( 8 ) ::    Q                     ! intermediate term

      REAL( 8 ) ::    T1, T2, T3            ! intermediate terms

      REAL( 8 ) ::    L21_INV               ! reciprocal of HNO4 loss term

C**********************************************************************


c..stoichiometric coefficient for production of HOX from O3 via O1D
      O1D_S =                 RKI(  20 )                        ! O1D=0.2000D+01*OH
     &      +                 RKI(  21 )                        ! O1D=O3P

      OH_S =    2.0000D+00 * RKI(  20 )                        ! O1D=0.2000D+01*OH

      OH_S  = OH_S / O1D_S


c..Production of HO from HO2 (r4,5 terms )
      R4_5 =                 RKI(  31 ) * YC( NO           )   ! HO2+NO=OH+NO2
     &     +                 RKI(  36 ) * YC( O3           )   ! HO2+O3=OH
     &     +    8.0000D-01 * RKI(  39 ) * YC( NO3          )   ! HO2+NO3=0.8000D+...
     &     +    4.4000D-01 * RKI(  67 ) * YC( MECO3        )   ! HO2+MECO3=...
     &     +    4.4000D-01 * RKI(  77 ) * YC( RCO3         )   ! HO2+RCO3=0.4400D+...
     &     +    4.4000D-01 * RKI(  88 ) * YC( BZCO3        )   ! HO2+BZCO3=...
     &     +    4.4000D-01 * RKI( 100 ) * YC( MACO3        )   ! HO2+MACO3=...
     &     +    4.4000D-01 * RKI( 212 ) * YC( HCOCO3       )   ! HO2+HCOCO3=...
     &     +    1.7500D-01 * RKI( 287 ) * YC( IEPOXOO      )   ! HO2+IEPOXOO=...
     &     +                 RKI( 338 ) * YC( CL           )   ! HO2+CL=OH+CLO

      R4_5  = R4_5 * DTC

c..Production of HO from HONO (r4,19 terms )
      R4_19 =                 RKI(  23 )                        ! HONO=OH+NO

      R4_19 = R4_19 * DTC

c..Remaining HO production
      P4 =   OH_S  *       RXRAT(  18 )      ! O1D=0.2000D+01*OH
     &   +                 RXRAT(  28 )      ! HNO3=OH+NO2
     &   +    3.9000D-01 * RXRAT(  34 )      ! HNO4=0.3900D+00*OH+0.6100D+...
     &   +    2.0000D+00 * RXRAT(  41 )      ! HO2H=0.2000D+01*OH
     &   +                 RXRAT( 116 )      ! xOH=OH
     &   +                 RXRAT( 151 )      ! COOH=OH+HO2+HCHO
     &   +                 RXRAT( 153 )      ! ROOH=OH+HO2+RCHO
     &   +                 RXRAT( 155 )      ! R6OOH=OH+0.1420D+00*HO2+...
     &   +                 RXRAT( 157 )      ! RAOOH=OH+HO2+0.5000D+...
     &   +    8.2600D-01 * RXRAT( 175 )      ! AFG1+O3=0.8260D+00*OH+...
     &   +    8.2600D-01 * RXRAT( 178 )      ! AFG2+O3=0.8260D+00*OH+...
     &   +    4.7100D-01 * RXRAT( 181 )      ! AFG3+O3=0.4710D+00*OH+...
     &   +    2.0800D-01 * RXRAT( 183 )      ! MACR+O3=0.2080D+00*OH+...
     &   +    3.3000D-01 * RXRAT( 186 )      ! MACR=0.3300D+00*OH+0.6700D+...
     &   +    1.6400D-01 * RXRAT( 188 )      ! MVK+O3=0.1640D+00*OH+...
     &   +    2.8500D-01 * RXRAT( 192 )      ! IPRD+O3=0.2850D+00*OH+...
     &   +    3.3000D-01 * RXRAT( 203 )      ! ACROLEIN+O3=0.3300D+00*OH+...
     &   +    1.7800D-01 * RXRAT( 206 )      ! ACROLEIN=0.1780D+00*OH+...
     &   +                 RXRAT( 208 )      ! CCOOOH=OH+CO2+MEO2
     &   +    1.6000D-01 * RXRAT( 268 )      ! ETHENE+O3=0.1600D+00*OH+...
     &   +    3.5000D-01 * RXRAT( 272 )      ! PROPENE+O3=0.3500D+00*OH+...
     &   +    8.0000D-02 * RXRAT( 276 )      ! BUTADIENE13+O3=0.8000D-...
     &   +    2.6600D-01 * RXRAT( 280 )      ! ISOPRENE+O3=0.2660D+00*OH+...
     &   +                 RXRAT( 285 )      ! ISOPOOH=OH+0.9100D+00*HO2+...
     &   +    9.5000D-01 * RXRAT( 287 )      ! HO2+IEPOXOO=0.1125D+01*OH+...
     &   +    1.2500D-01 * RXRAT( 288 )      ! IEPOXOO+NO=0.1250D+00*OH+...
     &   +    6.3000D-02 * RXRAT( 289 )      ! IEPOXOO+MEO2=0.6300D-01*OH+...
     &   +    6.3000D-02 * RXRAT( 290 )      ! IEPOXOO+RO2C=0.6300D-01*OH+...
     &   +    1.2500D-01 * RXRAT( 291 )      ! IEPOXOO+MECO3=0.1250D+...
     &   +    7.2800D-01 * RXRAT( 293 )      ! APIN+O3=0.7280D+00*OH+...
     &   +    5.0000D-01 * RXRAT( 297 )      ! ACETYLENE+O3=0.5000D+00*OH+...
     &   +    1.2800D-01 * RXRAT( 312 )      ! OLE1+O3=0.1280D+00*OH+...
     &   +    4.4300D-01 * RXRAT( 316 )      ! OLE2+O3=0.4430D+00*OH+...
     &   +    4.9900D-01 * RXRAT( 323 )      ! TERP+O3=0.4990D+00*OH+...
     &   +    4.9900D-01 * RXRAT( 327 )      ! SESQ+O3=0.4990D+00*OH+...
     &   +                 RXRAT( 348 )      ! HOCL=OH+CL
     &   +    4.1400D-01 * RXRAT( 364 )      ! ROOH+CL=0.4140D+00*OH+HCL+...
     &   +    1.4500D-01 * RXRAT( 365 )      ! R6OOH+CL=0.1450D+00*OH+HCL+...
     &   +    1.3900D-01 * RXRAT( 366 )      ! RAOOH+CL=0.1390D+00*OH+...

      P4 = YC0( OH ) + P4 * DTC

c..Production of HO2 from OH ( r5,4 terms )
      R5_4 =                 RKI(  26 ) * YC( NO3          )   ! OH+NO3=HO2+NO2
     &     +                 RKI(  29 ) * YC( CO           )   ! OH+CO=HO2+CO2
     &     +                 RKI(  30 ) * YC( O3           )   ! OH+O3=HO2
     &     +                 RKI(  42 ) * YC( HO2H         )   ! OH+HO2H=HO2
     &     +                 RKI(  44 ) * YC( SO2          )   ! OH+SO2=HO2+SULF+...
     &     +                 RKI(  45 )                        ! OH=HO2
     &     +                 RKI( 134 ) * YC( HCHO         )   ! OH+HCHO=HO2+CO
     &     +                 RKI( 146 ) * YC( MEOH         )   ! OH+MEOH=HO2+HCHO
     &     +                 RKI( 147 ) * YC( HCOOH        )   ! OH+HCOOH=HO2+CO2
     &     +    1.4800D-01 * RKI( 156 ) * YC( RAOOH        )   ! OH+RAOOH=0.1480D+...
     &     +    7.0000D-01 * RKI( 160 ) * YC( GLY          )   ! OH+GLY=0.7000D+...
     &     +    4.7200D-01 * RKI( 195 ) * YC( PRD2         )   ! OH+PRD2=0.4720D+...
     &     +    1.8900D-01 * RKI( 197 ) * YC( RNO3         )   ! OH+RNO3=0.1890D+...
     &     +    3.0000D-01 * RKI( 296 ) * YC( ACETYLENE    )   ! OH+ACETYLENE=...
     &     +    5.7000D-01 * RKI( 298 ) * YC( BENZENE      )   ! OH+BENZENE=...
     &     +    1.8100D-01 * RKI( 299 ) * YC( TOLUENE      )   ! OH+TOLUENE=...
     &     +    1.5900D-01 * RKI( 300 ) * YC( MXYL         )   ! OH+MXYL=0.1590D+...
     &     +    1.6100D-01 * RKI( 301 ) * YC( OXYL         )   ! OH+OXYL=0.1610D+...
     &     +    1.5900D-01 * RKI( 302 ) * YC( PXYL         )   ! OH+PXYL=0.1590D+...
     &     +    2.2000D-02 * RKI( 303 ) * YC( TMBENZ124    )   ! OH+TMBENZ124=...
     &     +    9.5000D-01 * RKI( 304 ) * YC( ETOH         )   ! OH+ETOH=0.9500D+...
     &     +    1.2300D-01 * RKI( 319 ) * YC( ARO1         )   ! OH+ARO1=0.1230D+...
     &     +    7.7000D-02 * RKI( 320 ) * YC( ARO2MN       )   ! OH+ARO2MN=...
     &     +    7.7000D-02 * RKI( 321 ) * YC( NAPHTHAL     )   ! OH+NAPHTHAL=...

      R5_4  = R5_4 * DTC

c..Production of HO2 from HNO4 (r5,21 term )
      R5_21 =                 RKI(  33 )                        ! HNO4=HO2+NO2
     &      +    6.1000D-01 * RKI(  34 )                        ! HNO4=0.6100D+...

      R5_21 = R5_21 * DTC

c..Remaining HO2 production terms
      P5 =                 RXRAT(  46 )      ! MEO2+NO=HO2+HCHO+NO2
     &   +                 RXRAT(  49 )      ! MEO2+NO3=HO2+HCHO+NO2
     &   +    2.0000D+00 * RXRAT(  51 )      ! MEO2+MEO2=0.2000D+01*HO2+...
     &   +    5.0000D-01 * RXRAT(  55 )      ! RO2C+MEO2=0.5000D+00*HO2+...
     &   +    5.0000D-01 * RXRAT(  60 )      ! RO2XC+MEO2=0.5000D+00*HO2+...
     &   +    9.0000D-01 * RXRAT(  69 )      ! MECO3+MEO2=0.9000D+00*HO2+...
     &   +                 RXRAT(  79 )      ! RCO3+MEO2=HO2+HCHO+RO2C+...
     &   +                 RXRAT(  90 )      ! BZCO3+MEO2=HO2+HCHO+RO2C+...
     &   +                 RXRAT( 102 )      ! MACO3+MEO2=HO2+0.2000D+...
     &   +                 RXRAT( 114 )      ! xHO2=HO2
     &   +    2.0000D+00 * RXRAT( 132 )      ! HCHO=0.2000D+01*HO2+CO
     &   +                 RXRAT( 135 )      ! HCHO+NO3=HO2+HNO3+CO
     &   +                 RXRAT( 137 )      ! CCHO=HO2+CO+MEO2
     &   +                 RXRAT( 140 )      ! RCHO=HO2+xHO2+yROOH+xCCHO+...
     &   +                 RXRAT( 151 )      ! COOH=HO2+OH+HCHO
     &   +                 RXRAT( 153 )      ! ROOH=HO2+OH+RCHO
     &   +    1.4200D-01 * RXRAT( 155 )      ! R6OOH=0.1420D+00*HO2+OH+...
     &   +                 RXRAT( 157 )      ! RAOOH=HO2+OH+0.5000D+...
     &   +    2.0000D+00 * RXRAT( 158 )      ! GLY=0.2000D+01*HO2+0.2000D+...
     &   +    7.0000D-01 * RXRAT( 161 )      ! GLY+NO3=0.7000D+00*HO2+...
     &   +                 RXRAT( 162 )      ! MGLY=HO2+CO+MECO3
     &   +    5.2200D-01 * RXRAT( 175 )      ! AFG1+O3=0.5220D+00*HO2+...
     &   +    1.0230D+00 * RXRAT( 176 )      ! AFG1=0.1023D+01*HO2+...
     &   +    5.2200D-01 * RXRAT( 178 )      ! AFG2+O3=0.5220D+00*HO2+...
     &   +    5.5400D-01 * RXRAT( 181 )      ! AFG3+O3=0.5540D+00*HO2+...
     &   +    1.0800D-01 * RXRAT( 183 )      ! MACR+O3=0.1080D+00*HO2+...
     &   +    6.7000D-01 * RXRAT( 186 )      ! MACR=0.6700D+00*HO2+...
     &   +    6.4000D-02 * RXRAT( 188 )      ! MVK+O3=0.6400D-01*HO2+...
     &   +    4.0000D-01 * RXRAT( 192 )      ! IPRD+O3=0.4000D+00*HO2+...
     &   +    1.2330D+00 * RXRAT( 194 )      ! IPRD=0.1233D+01*HO2+...
     &   +    3.4400D-01 * RXRAT( 198 )      ! RNO3=0.3440D+00*HO2+...
     &   +    2.0000D+00 * RXRAT( 200 )      ! HOCCHO=0.2000D+01*HO2+CO+HCHO
     &   +    8.3000D-01 * RXRAT( 203 )      ! ACROLEIN+O3=0.8300D+00*HO2+...
     &   +    1.0660D+00 * RXRAT( 206 )      ! ACROLEIN=0.1066D+01*HO2+...
     &   +                 RXRAT( 210 )      ! HCOCO3+NO=HO2+CO+CO2+NO2
     &   +                 RXRAT( 211 )      ! HCOCO3+NO2=HO2+CO+CO2+NO3
     &   +                 RXRAT( 248 )      ! zRNO3=HO2+PRD2
     &   +    1.6000D-01 * RXRAT( 268 )      ! ETHENE+O3=0.1600D+00*HO2+...
     &   +    8.0000D-01 * RXRAT( 270 )      ! ETHENE+O3P=0.8000D+00*HO2+...
     &   +    1.6500D-01 * RXRAT( 272 )      ! PROPENE+O3=0.1650D+00*HO2+...
     &   +    8.0000D-02 * RXRAT( 276 )      ! BUTADIENE13+O3=0.8000D-...
     &   +    2.5000D-01 * RXRAT( 278 )      ! BUTADIENE13+O3P=0.2500D+...
     &   +    6.6000D-02 * RXRAT( 280 )      ! ISOPRENE+O3=0.6600D-01*HO2+...
     &   +    9.1000D-01 * RXRAT( 285 )      ! ISOPOOH=0.9100D+00*HO2+OH+...
     &   +    8.2500D-01 * RXRAT( 288 )      ! IEPOXOO+NO=0.8250D+00*HO2+...
     &   +    9.1300D-01 * RXRAT( 289 )      ! IEPOXOO+MEO2=0.9130D+...
     &   +    4.1300D-01 * RXRAT( 290 )      ! IEPOXOO+RO2C=0.4130D+...
     &   +    8.2500D-01 * RXRAT( 291 )      ! IEPOXOO+MECO3=0.8250D+...
     &   +    9.0000D-03 * RXRAT( 293 )      ! APIN+O3=0.9000D-02*HO2+...
     &   +    1.5000D+00 * RXRAT( 297 )      ! ACETYLENE+O3=0.1500D+...
     &   +    9.5000D-02 * RXRAT( 312 )      ! OLE1+O3=0.9500D-01*HO2+...
     &   +    9.4000D-02 * RXRAT( 316 )      ! OLE2+O3=0.9400D-01*HO2+...
     &   +    7.8000D-02 * RXRAT( 323 )      ! TERP+O3=0.7800D-01*HO2+...
     &   +    7.8000D-02 * RXRAT( 327 )      ! SESQ+O3=0.7800D-01*HO2+...
     &   +                 RXRAT( 351 )      ! CL=HO2+HCL
     &   +                 RXRAT( 352 )      ! HCHO+CL=HO2+HCL+CO
     &   +                 RXRAT( 354 )      ! MEOH+CL=HO2+HCHO+HCL
     &   +    5.5000D-02 * RXRAT( 358 )      ! RNO3+CL=0.5500D-01*HO2+...
     &   +    3.1400D-01 * RXRAT( 359 )      ! PRD2+CL=0.3140D+00*HO2+HCL+...
     &   +    6.3000D-01 * RXRAT( 360 )      ! GLY+CL=0.6300D+00*HO2+HCL+...
     &   +    1.4800D-01 * RXRAT( 366 )      ! RAOOH+CL=0.1480D+00*HO2+...
     &   +    8.4000D-02 * RXRAT( 370 )      ! IPRD+CL=0.8400D-01*HO2+...
     &   +                 RXRAT( 371 )      ! CLCCHO=HO2+CO+RO2C+xCL+...
     &   +                 RXRAT( 387 )      ! ACETYLENE+CL=HO2+CO
     &   +    6.8800D-01 * RXRAT( 393 )      ! ETOH+CL=0.6880D+00*HO2+HCL+...

      P5 = YC0( HO2 ) + P5 * DTC

c..Production of HONO from OH (r19,4 terms )

      R19_4 =   RKI(  22 ) * YC( NO           ) * DTC      ! OH+NO=HONO

c..Remaining HONO production terms
      P19 =                 RXRAT( 169 )      ! NPHE=HONO
     &    +    5.0000D-01 * RXRAT( 429 )      ! NO2=0.5000D+00*HONO+...

      P19 = YC0( HONO ) + P19 * DTC

c..Production of HNO4 from HO2 (r21,5 term )

      R21_5 =   RKI(  32 ) * YC( NO2          ) * DTC      ! HO2+NO2=HNO4

c..Remaining HNO4 production terms
      P21   =   YC0( HNO4 )

c..OH loss terms not in R5_4 & R19_4
      L4 =                 RKI(  24 ) * YC( HONO         )   ! OH+HONO=NO2
     &   +                 RKI(  25 ) * YC( NO2          )   ! OH+NO2=HNO3
     &   +                 RKI(  27 ) * YC( HNO3         )   ! OH+HNO3=NO3
     &   +                 RKI(  35 ) * YC( HNO4         )   ! OH+HNO4=NO2
     &   +                 RKI(  43 ) * YC( HO2          )   ! OH+HO2=
     &   +                 RKI( 136 ) * YC( CCHO         )   ! OH+CCHO=MECO3
     &   +                 RKI( 139 ) * YC( RCHO         )   ! OH+RCHO=0.9650D+...
     &   +                 RKI( 142 ) * YC( ACETONE      )   ! OH+ACETONE=RO2C+...
     &   +                 RKI( 144 ) * YC( MEK          )   ! OH+MEK=0.9670D+...
     &   +                 RKI( 148 ) * YC( CCOOH        )   ! OH+CCOOH=0.5090D+...
     &   +                 RKI( 149 ) * YC( RCOOH        )   ! OH+RCOOH=RO2C+...
     &   +    7.0000D-01 * RKI( 150 ) * YC( COOH         )   ! OH+COOH=0.3000D+...
     &   +    2.5600D-01 * RKI( 152 ) * YC( ROOH         )   ! OH+ROOH=0.7440D+...
     &   +    1.6000D-01 * RKI( 154 ) * YC( R6OOH        )   ! OH+R6OOH=0.8400D+...
     &   +    7.1300D-01 * RKI( 156 ) * YC( RAOOH        )   ! OH+RAOOH=0.1480D+...
     &   +    3.0000D-01 * RKI( 160 ) * YC( GLY          )   ! OH+GLY=0.7000D+...
     &   +                 RKI( 163 ) * YC( MGLY         )   ! OH+MGLY=CO+MECO3
     &   +                 RKI( 166 ) * YC( CRES         )   ! OH+CRES=0.2000D+...
     &   +                 RKI( 168 ) * YC( NPHE         )   ! OH+NPHE=BZO
     &   +                 RKI( 171 ) * YC( BALD         )   ! OH+BALD=BZCO3
     &   +                 RKI( 174 ) * YC( AFG1         )   ! OH+AFG1=0.2170D+...
     &   +                 RKI( 177 ) * YC( AFG2         )   ! OH+AFG2=0.2170D+...
     &   +                 RKI( 180 ) * YC( AFG3         )   ! OH+AFG3=0.2060D+...
     &   +                 RKI( 182 ) * YC( MACR         )   ! OH+MACR=0.5000D+...
     &   +                 RKI( 187 ) * YC( MVK          )   ! OH+MVK=0.9750D+...
     &   +                 RKI( 191 ) * YC( IPRD         )   ! OH+IPRD=0.2890D+...
     &   +    5.2800D-01 * RKI( 195 ) * YC( PRD2         )   ! OH+PRD2=0.4720D+...
     &   +    8.1100D-01 * RKI( 197 ) * YC( RNO3         )   ! OH+RNO3=0.1890D+...
     &   +                 RKI( 199 ) * YC( HOCCHO       )   ! OH+HOCCHO=MECO3
     &   +                 RKI( 202 ) * YC( ACROLEIN     )   ! OH+ACROLEIN=...
     &   +                 RKI( 207 ) * YC( CCOOOH       )   ! OH+CCOOOH=...
     &   +                 RKI( 209 ) * YC( RCOOOH       )   ! OH+RCOOOH=...
     &   +                 RKI( 266 )                        ! OH=MEO2
     &   +                 RKI( 267 ) * YC( ETHENE       )   ! OH+ETHENE=xHO2+...
     &   +                 RKI( 271 ) * YC( PROPENE      )   ! OH+PROPENE=...
     &   +                 RKI( 275 ) * YC( BUTADIENE13  )   ! OH+BUTADIENE13=...
     &   +                 RKI( 279 ) * YC( ISOPRENE     )   ! OH+ISOPRENE=...
     &   +    3.9000D-01 * RKI( 284 ) * YC( ISOPOOH      )   ! OH+ISOPOOH=...
     &   +                 RKI( 286 ) * YC( IEPOX        )   ! OH+IEPOX=IEPOXOO
     &   +                 RKI( 292 ) * YC( APIN         )   ! OH+APIN=0.7990D+...
     &   +    3.1400D-01 * RKI( 298 ) * YC( BENZENE      )   ! OH+BENZENE=...
     &   +    5.0700D-01 * RKI( 299 ) * YC( TOLUENE      )   ! OH+TOLUENE=...
     &   +    6.0200D-01 * RKI( 300 ) * YC( MXYL         )   ! OH+MXYL=0.1590D+...
     &   +    6.4100D-01 * RKI( 301 ) * YC( OXYL         )   ! OH+OXYL=0.1610D+...
     &   +    5.6300D-01 * RKI( 302 ) * YC( PXYL         )   ! OH+PXYL=0.1590D+...
     &   +    7.4800D-01 * RKI( 303 ) * YC( TMBENZ124    )   ! OH+TMBENZ124=...
     &   +    5.0000D-02 * RKI( 304 ) * YC( ETOH         )   ! OH+ETOH=0.9500D+...
     &   +                 RKI( 305 ) * YC( ALK1         )   ! OH+ALK1=xHO2+...
     &   +                 RKI( 306 ) * YC( ALK2         )   ! OH+ALK2=0.9650D+...
     &   +                 RKI( 307 ) * YC( ALK3         )   ! OH+ALK3=0.6950D+...
     &   +                 RKI( 308 ) * YC( ALK4         )   ! OH+ALK4=0.8300D+...
     &   +                 RKI( 309 ) * YC( ALK5         )   ! OH+ALK5=0.6470D+...
     &   +                 RKI( 311 ) * YC( OLE1         )   ! OH+OLE1=0.8710D+...
     &   +                 RKI( 315 ) * YC( OLE2         )   ! OH+OLE2=0.9120D+...
     &   +    6.7500D-01 * RKI( 319 ) * YC( ARO1         )   ! OH+ARO1=0.1230D+...
     &   +    7.4500D-01 * RKI( 320 ) * YC( ARO2MN       )   ! OH+ARO2MN=...
     &   +    7.4500D-01 * RKI( 321 ) * YC( NAPHTHAL     )   ! OH+NAPHTHAL=...
     &   +                 RKI( 322 ) * YC( TERP         )   ! OH+TERP=0.7340D+...
     &   +                 RKI( 326 ) * YC( SESQ         )   ! OH+SESQ=0.7340D+...
     &   +                 RKI( 350 ) * YC( HCL          )   ! OH+HCL=CL
     &   +                 RKI( 372 ) * YC( CLCCHO       )   ! OH+CLCCHO=RCO3

      L4    = 1.0D0 + L4 * DTC + R5_4 + R19_4

c..HO2 loss terms not included in R4_5 & R21_5 (except for HO2+HO2 )
      L5 =    2.0000D-01 * RKI(  39 ) * YC( NO3          )   ! HO2+NO3=0.8000D+...
     &   +                 RKI(  43 ) * YC( OH           )   ! HO2+OH=
     &   +                 RKI(  47 ) * YC( MEO2         )   ! HO2+MEO2=COOH
     &   +                 RKI(  48 ) * YC( MEO2         )   ! HO2+MEO2=HCHO
     &   +                 RKI(  53 ) * YC( RO2C         )   ! HO2+RO2C=
     &   +                 RKI(  58 ) * YC( RO2XC        )   ! HO2+RO2XC=
     &   +    5.6000D-01 * RKI(  67 ) * YC( MECO3        )   ! HO2+MECO3=...
     &   +    5.6000D-01 * RKI(  77 ) * YC( RCO3         )   ! HO2+RCO3=0.4400D+...
     &   +    5.6000D-01 * RKI(  88 ) * YC( BZCO3        )   ! HO2+BZCO3=...
     &   +    5.6000D-01 * RKI( 100 ) * YC( MACO3        )   ! HO2+MACO3=...
     &   +                 RKI( 112 ) * YC( BZO          )   ! HO2+BZO=CRES
     &   +    1.2000D-01 * RKI( 212 ) * YC( HCOCO3       )   ! HO2+HCOCO3=...
     &   +                 RKI( 337 ) * YC( CL           )   ! HO2+CL=HCL
     &   +                 RKI( 347 ) * YC( CLO          )   ! HO2+CLO=HOCL

      L5    = 1.0D0 + L5 * DTC + R4_5 + R21_5

c..HONO loss terms not included in R4_19
      L19 =                 RKI(  24 ) * YC( OH           )   ! HONO+OH=NO2

      L19   = 1.0D0 + L19 * DTC + R4_19

c..HNO4 loss terms not inluded in R5_21
      L21 =    3.9000D-01 * RKI(  34 )                        ! HNO4=0.6100D+...
     &    +                 RKI(  35 ) * YC( OH           )   ! HNO4+OH=NO2

      L21   = 1.0D0 + L21 * DTC + R5_21

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Solution section
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c..compute terms used to calculate a,b & c
      L21_INV = 1.0D0 / L21
      T1 = 1.0D0 / ( L4 * L19 - R4_19 * R19_4 )
      T2 = R5_4 * T1
      T3 = R5_21 * L21_INV

c..solve quadratic equation for HO2
      A = 2.0D0 * ( RKI(  37 ) + RKI(  38 ) ) * DTC

      B = L5 - T3 * R21_5 - T2 * R4_5 * L19

      C = P5 + T3 * P21 + T2 * ( P4 * L19 + P19 * R4_19 )

      Q = -0.5D0 * ( B + SIGN( 1.0D0, B ) * SQRT( B * B + 4.0D0 * A * C ) )

      YCP( HO2 ) = MAX( Q / A , -C / Q  )

c..compute remaining species concentrations
      YCP( OH ) = ( ( P4 + R4_5 * YCP( HO2 ) ) * L19 + R4_19 * P19 ) * T1

      YCP( HNO4 ) = ( P21 + R21_5 * YCP( HO2 ) ) * L21_INV

      YCP( HONO ) = ( P19 + R19_4 * YCP( OH ) ) / L19

      RETURN

      END


