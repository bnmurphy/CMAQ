
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
C $Header: /project/yoj/arc/CCTM/src/aero/aero6/NPF.F,v 1.8 2012/01/19 13:21:28 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C routines for aerosol new particle formation 


C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE NPF_MULTI ( CS, ASURF, RH, PRESS, TEMP, XSA, SO4RATE,
     &                       XNH3, XDMA, XMBTCA, XELVOC,
     &                       JRATE, M2RATE, DMDT_SO4_TOT, DMDT_NH3_TOT,
     &                       DMDT_DMA_TOT, DMDT_MBTCA_TOT, DMDT_ELVOC_TOT)
 
      USE AERO_DATA
      USE AEROMET_DATA   ! Includes CONST.EXT
      USE PRECURSOR_DATA, ONLY: PRECURSOR_MW, SULF_IDX, NH3_IDX 
      USE UTILIO_DEFN
      
      IMPLICIT NONE

!  Comments:
!  This subroutine calculates the New Particle Formation Rate,
!    NPF, using published empirical rate coefficients. Two paths
!    are calculated; Sulfuric Acid (SA) with Dimathylamine, DMA, using
!    the method of Glasoe et al. (2014), and SA with and organics
!    species
!    using the method of Paasonen et al. (2010).
!    FSB code August 24, 2015 by Dr. Francis S. Binkowski
!       University of North Carolina at Chapel Hill, NC
!    
!     References:
!     Glasoe, W.A. et al. 2014, Sulfuric acid nucleation: An
!      experimental study of the effect of seven bases, Journal
!      of Geophysical Research: Atmospheres, doe: 10.1002/2014JD0
!
!     Paasonen, P et al., 2010, On the role of sulphuric acid and
!      low-volatility organic vapors in the initial steps of
!      atmospheric new particle production, Atmospheric Chemistry
!      and Physics, doi: 105194/acp-10-11223-2010.
!
!     Riccobono,F. et al. (2014) Oxidation of Biogenic Emissions
!     Contribute to
!     nucleation of atmospheric particles, Science,doi:
!     10.1126/science.1243527
!
!     Yu, F., Ion-mediated nucleation in the atmosphere: Key controlling parameters,
!     implications, and look-up table, J. Geophy. Res., 115, D03206, 
!     doi:10.1029/2009JD012630, 2010.
!
!     Yu, F., From molecular clusters to nanoparticles: Second-generation ion-mediated
!     nucleation model, Atmos. Chem. Phys.,6, 5193-5211, 2006.
!
!     Yu, F. and R. P. Turco, From molecular clusters to nanoparticles: The role of 
!     ambient ionization in tropospheric aerosol formation, J. Geophy. Res., 106, 
!     4797-4814, 2001.
!
!
!  REVISION HISTORY:
!  2015/09/15 Installed by Ben Murphy. In this first version,
!             only H2SO4 is passed in to the routine. The other
!             precursor species are assumed constant.
!

      REAL, INTENT(IN)  :: PRESS     ! ambient pressure [PA]
      REAL, INTENT(IN)  :: TEMP      ! ambient temperature [K]
      REAL, INTENT(IN)  :: CS        ! condensation sink [s-1]
      REAL, INTENT(IN)  :: ASURF     ! aerosol surface area [m2 m-3]
      REAL, INTENT(IN)  :: RH        ! ambient Relative Humidity
      REAL, INTENT(IN)  :: XSA       ! H2SO4 volume concentration [ug/m3]
      REAL, INTENT(IN)  :: XNH3      ! NH3 volume concentration [ug/m3]
      REAL, INTENT(IN)  :: XDMA      ! DMA volume concentration [ug/m3]
      REAL, INTENT(IN)  :: XMBTCA    ! MBTCA volume concentration [ug/m3]
      REAL, INTENT(IN)  :: XELVOC    ! ELVOC volume concentration [ug/m3]
      REAL, INTENT(IN)  :: SO4RATE   ! H2SO4 production rate [ug/m3/s]
      REAL              :: DMA_ppmv  ! DMA mixing ratio [ppmv]
      REAL              :: XORG_ppmv ! organic mixing ratio [ppmv]
      REAL              :: YORG_ppmv ! organic mixing ratio [ppmv]
      REAL(8), INTENT(OUT) :: JRATE     ! NPF RATE [N / m**3 / sec ]
      REAL(8), INTENT(OUT) :: M2RATE    ! RATE OF 2nd Moment Chanfe [ m**2 / m**3 / sec ]
      REAL(8), INTENT(OUT) :: DMDT_SO4_TOT  ! NPF RATE [ ug / m**3 / sec ]
      REAL(8), INTENT(OUT) :: DMDT_NH3_TOT  ! NPF RATE [ ug / m**3 / sec ]
      REAL(8), INTENT(OUT) :: DMDT_DMA_TOT  ! NPF RATE [ ug / m**3 / sec ]
      REAL(8), INTENT(OUT) :: DMDT_MBTCA_TOT  ! NPF RATE [ ug / m**3 / sec ]
      REAL(8), INTENT(OUT) :: DMDT_ELVOC_TOT  ! NPF RATE [ ug / m**3 / sec ]

! coefficient from Glasoe et al. (2014)
      REAL(8), PARAMETER :: K_GLASOE = 1.0d-20 ! (SA) * (SA) * (DMA)
      REAL(8), PARAMETER :: Dg_GLASOE = 3.6d-9 ! New Particle Diam [m]
      REAL(8), PARAMETER :: Sg_GLASOE = 1.28 ! New Particle Standard Dev
      REAL(8), PARAMETER :: Dg_GLASOEsq = Dg_GLASOE * Dg_GLASOE  !Diameter Squared [m^2/Particle]
      REAL(8), PARAMETER :: Dg_GLASOEcu = Dg_GLASOE * Dg_GLASOEsq !New Particle Volume Cubed [m3]

! coefficients from Passonen et al (2010) 
      REAL(8), PARAMETER   ::  KPAAS_2 = 2.6d-14      ! (SA) * (XORG)
      REAL(8), PARAMETER   ::  Dg_PAAS2 = 3.6d-9 ! New Particle Diam (m)
      REAL(8), PARAMETER   ::  Sg_PAAS2 = 1.28 ! New Particle Standard Dev
      REAL(8), PARAMETER   ::  Dg_PAAS2sq = Dg_PAAS2 * Dg_PAAS2  !Diameter Squared [m^2/Particle]
      REAL(8), PARAMETER   ::  Dg_PAAS2cu = Dg_PAAS2 * Dg_PAAS2sq !New Particle Volume Cubed [m3]
      
      REAL(8), PARAMETER   ::  KPAAS_3 = 0.037d-14    ! (YORG) * (YORG)
      REAL(8), PARAMETER   ::  Dg_PAAS3 = 3.6d-9 ! New Particle Diam (m)
      REAL(8), PARAMETER   ::  Sg_PAAS3 = 1.28 ! New Particle Standard Dev
      REAL(8), PARAMETER   ::  Dg_PAAS3sq = Dg_PAAS3 * Dg_PAAS3  !Diameter Squared [m^2/Particle]
      REAL(8), PARAMETER   ::  Dg_PAAS3cu = Dg_PAAS3 * Dg_PAAS3sq !New Particle Volume Cubed [m3]

! Coefficients from Kirkby et al (2016)
      REAL(8), PARAMETER   ::  Dg_ORG   = 1.70d-9 ! New Particle Diam (m)
      REAL(8), PARAMETER   ::  Sg_ORG   = 1.28 ! New Particle Standard Dev
      REAL(8), PARAMETER   ::  Dg_ORGsq = Dg_ORG * Dg_ORG  !Diameter Squared [m^2/Particle]
      REAL(8), PARAMETER   ::  Dg_ORGcu = Dg_ORG * Dg_ORGsq !New Particle Volume Cubed [m3]
      REAL(8), PARAMETER   ::  K_KIRKBY  = 0.04001  ! [ cm-3 s-1 ]
      REAL(8), PARAMETER   ::  E1_KIRKBY = 1.848    ! Exponential Term
      REAL(8), PARAMETER   ::  E2_KIRKBY = 0.1863   ! Exponential Term
      
! Coefficient from Riccobono et al. (2014)
      REAL(8), PARAMETER   ::  K_RICCO = 3.27d-21     ! (SA)*(SA) * (XORG) 
      REAL(8), PARAMETER   ::  Dg_RICCO = 3.6d-9 ! New Particle Diam (m)
      REAL(8), PARAMETER   ::  Sg_RICCO = 1.28 ! New Particle Standard Dev
      REAL(8), PARAMETER   ::  Dg_RICCOsq = Dg_RICCO * Dg_RICCO  !Diameter Squared [m^2/Particle]
      REAL(8), PARAMETER   ::  Dg_RICCOcu = Dg_RICCO * Dg_RICCOsq !Particle Volume Cubed [m3]

! Variables for Ion-Mediated Calculations
      REAL(8) :: ION_RATE  ! Ionization Rate (ion-pairs/cm3)
           
! Define internal variables. See usage.       
      REAL(8)            :: Dg_BINARYsq, Dg_TERNARYsq, Dg_AMINEsq
      REAL(8)            :: Dg_IMNSO4sq, Dg_IMNORGsq
      REAL(8)            :: DMDT_SO4_BINARY, DMDT_SO4_TERNARY, DMDT_SO4_AMINE
      REAL(8)            :: DMDT_SO4_GLASOE, DMDT_SO4_PAAS2, DMDT_SO4_IMN
      REAL(8)            :: DMDT_DMA, DMDT_XORG, DMDT_YORG
      REAL(8)            :: DMDT_NH3, DMDT_MBTCA, DMDT_ELVOC
      REAL(8)            :: SA_RO, DMA_RO, XORG_RO, YORG_RO
      REAL(8)            :: SA_MW, DMA_MW, XORG_MW, YORG_MW
      REAL(8)            :: NH3_MW, MBTCA_MW, ELVOC_MW
      REAL(8)            :: M_GLASOE, M_PAAS2, M_PAAS3, M_RICCO, M_ORG
      REAL(8)            :: M3_GLASOE, M3_PAAS2, M3_PAAS3, M3_RICCO
      REAL(8)            :: M3_ORG
      REAL(8)            :: RO_GLASOE, RO_PAAS2, RO_PAAS3, RO_RICCO
      REAL(8)            :: RO_ORG
      REAL(8)            :: SA_M_BINARY, SA_M_TERNARY, SA_M_AMINE, SA_M_GLASOE, SA_M_PAAS2
      REAL(8)            :: SA_M_IMN
      REAL(8)            :: DMA_M_AMINE, DMA_M_GLASOE, XORG_M_PAAS2, YORG_M_PAAS3
      REAL(8)            :: MBTCA_M_PAAS2, ELVOC_M_PAAS2, ELVOC_M_ORG, ELVOC_M_IMN
      REAL(8)            :: NH3_M_TERNARY
      REAL(8)            :: M2_BINARY, M2_TERNARY, M2_AMINE, M2_GLASOE, M2_PAAS3
      REAL(8)            :: M2_PAAS2, M2_ORG, M2_RICCO, M2_IMN_SO4, M2_IMN_ORG
      REAL(8)            :: MWtot_GLASOE, MWtot_PAAS2, MWtot_PAAS3, MWtot_RICCO, MWtot_ORG
      REAL(8)            :: FAC0, FAC, FAC1, A1TOT
      REAL(8)            :: DMA, XORG, YORG, NH3, MBTCA, ELVOC
      REAL(8)            :: SASQRD, SAORGX, SA2ORGX, ORGYORGY
      REAL(8)            :: J0, J1, J2, J3, J4, J5, J6, J_PAAS3, J_RICCO, J_CGS, J_MKS
      REAL(8)            :: DIFF, A1DIFF, NORM_SO4
      LOGICAL            :: L_NAPARI, L_ACDC_TERNARY
      LOGICAL            :: L_GLASOE, L_ACDC_DMA

! Begin Code
      !Choose Formulation for H2SO4-NH4-H2O Ternary NPF
      !One of these should be false
      L_NAPARI       = .FALSE.       
      L_ACDC_TERNARY = .TRUE.  

      !Choose Formulation for H2SO4-Amine NPF
      L_GLASOE   = .FALSE.
      L_ACDC_DMA = .TRUE.      

! Initialize Diagnostic Data if Aerosol Diagnostics are Turned On
      JDIAG        = 0.0
      JM2DIAG      = 0.0
      JDGASDIAG    = 0.0

! Initialize NPF Rates
      J0 = 0.0
      J1 = 0.0
      J2 = 0.0
      J3 = 0.0
      J4 = 0.0
      J5 = 0.0
      J6 = 0.0
      DMDT_MBTCA_TOT = 0.0

! Set key variable FAC0 to convert from ug/m3 to number/cm**3
      FAC0     = AVO * 1.0E-9 * 1.0E-6 
! Set key variable Fac to convert from ppmv to [number / m**3] 
      FAC      = (1.0e-6) * AVO * PRESS / ( RGASUNIV * TEMP ) ! [Number / m**3]
! Adjust FAC to be [Number / cm**3] to be consistent with 
!     Paasonen et al (2010) and Glasoe et al. (2014)       
      FAC1     = (1.0E-6) * FAC ! [Number  / cm**3]
       
! Retrieve Intrinsic Properties of All Precursors
      SA_RO   = AEROSPC( ASO4_IDX )%DENSITY  !kg/m3
      DMA_RO  = 670.0  !kg/m3
      XORG_RO = 1250.0 !kg/m3  Reflecting MBTCA
      YORG_RO = 1200.0 !kg/m3  Paasonen et al. (2010); section 2.3.2 

      SA_MW = PRECURSOR_MW( SULF_IDX )/1000. !kg/mol
      NH3_MW = 0.0458  !kg/mol 
      DMA_MW = 0.0458  !kg/mol 
      XORG_MW = 0.204  !kg/mol Reflecting MBTCA
      YORG_MW = 0.168  !kg/mol Paasonen et al. (2010); section 2.3.2 (168 g/mol)
      MBTCA_MW= 0.204  !kg/mol
      ELVOC_MW= 0.248  !kg/mol Jokinen et al. (2015); ELVOCs (248 g/mol)

! Calculate the Mass-Averaged Density of each New Particle assuming
! ideal mixing of the component species. (the mixtures are probably NOT
! ideal)  [= kg / m3 ]
      MWtot_GLASOE = 2.0*SA_MW + DMA_MW
      RO_GLASOE    = ( 2.0*SA_MW*SA_RO + 1.0*DMA_MW*DMA_RO ) / MWtot_GLASOE
      !MWtot_PAAS2  = SA_MW + MBTCA_MW
      !RO_PAAS2    = ( 1.0*SA_MW*SA_RO + 1.0*MBTCA_MW*XORG_RO ) / MWtot_PAAS2
      MWtot_PAAS2  = SA_MW + ELVOC_MW
      RO_PAAS2     = ( 1.0*SA_MW*SA_RO + 1.0*ELVOC_MW*XORG_RO ) / MWtot_PAAS2
      !MWtot_PAAS3  = 2.0 * YORG_MW
      MWtot_PAAS3  = 2.0 * ELVOC_MW
      MWtot_ORG    = ELVOC_MW
      RO_ORG       = YORG_RO
      MWtot_RICCO  = 2.0*SA_MW + MBTCA_MW
      RO_RICCO     = ( 2.0*SA_MW*SA_RO + 1.0*MBTCA_MW*XORG_RO ) / MWtot_RICCO

! Calculate the mass of each new particle and the mass of each species
! per particle
      M3_GLASOE = Dg_GLASOEcu * exp( 9./2. * ( log( Sg_GLASOE )**2 ) )
      M3_PAAS2  = Dg_PAAS2cu  * exp( 9./2. * ( log( Sg_PAAS2  )**2 ) )
      M3_PAAS3  = Dg_PAAS3cu  * exp( 9./2. * ( log( Sg_PAAS3  )**2 ) )
      M3_ORG    = Dg_ORGcu    * exp( 9./2. * ( log( Sg_ORG    )**2 ) )
      M3_RICCO  = Dg_RICCOcu  * exp( 9./2. * ( log( Sg_RICCO  )**2 ) )

      M_GLASOE = M3_GLASOE * PI/6. * RO_GLASOE * 1.0E9  ![ug/particle]
      M_PAAS2  = M3_PAAS2  * PI/6. * RO_PAAS2  * 1.0E9  ![ug/particle]
      M_PAAS3  = M3_PAAS3  * PI/6. * RO_ORG    * 1.0E9  ![ug/particle]
      M_ORG    = M3_ORG    * PI/6. * RO_ORG    * 1.0E9  ![ug/particle]
      M_RICCO  = M3_RICCO  * PI/6. * RO_RICCO  * 1.0E9  ![ug/particle]

      SA_M_GLASOE  = SA_MW *  2.0/MWtot_GLASOE * M_GLASOE 
      SA_M_PAAS2   = SA_MW *  1.0/MWtot_PAAS2  * M_PAAS2   ![ug/particle]
      DMA_M_GLASOE = DMA_MW*  1.0/MWtot_GLASOE * M_GLASOE  ![ug/particle]
      !MBTCA_M_PAAS2 = MBTCA_MW* 1.0/MWtot_PAAS2  * M_PAAS2    ![ug/particle]
      ELVOC_M_PAAS2= ELVOC_MW* 1.0/MWtot_PAAS2 * M_PAAS2   ![ug/particle]
      !ELVOC_M_PAAS3= ELVOC_MW* 2.0/MWtot_PAAS3 * M_PAAS3   ![ug/particle]
      ELVOC_M_ORG  = ELVOC_MW* 1.0/MWtot_ORG   * M_ORG     ![ug/particle]

! Calculate number concentrations in [molecules / cm**3]      
      A1TOT    = XSA      * FAC0 / SA_MW    ! total number of SA 
      NH3      = XNH3     * FAC0 / NH3_MW   ! total number of DMA 
      DMA      = XDMA     * FAC0 / DMA_MW   ! total number of DMA 
      MBTCA    = XMBTCA   * FAC0 / MBTCA_MW ! total number of MBTCA
      ELVOC    = XELVOC   * FAC0 / ELVOC_MW ! total number of ELVOC
      DIFF     = A1TOT - DMA
      A1DIFF   = max( 0.0, DIFF )

! Calculate Binary H2SO4-H2O NPF from Vehkamaki et al. 2002
      CALL NPF_BINARY ( RH, TEMP, XSA, 
     &                   J0, SA_M_BINARY, Dg_BINARYsq)
 
! Calculate Ternary H2SO4-H2O-NH3 NPF 
      IF (L_NAPARI) THEN
        !Use Napari et al. 2004 Parameterization derived from theory
        CALL NPF_NAPARI ( PRESS, RH, TEMP, XSA, XNH3, 
     &                    J1, SA_M_TERNARY, NH3_M_TERNARY,
     &                    Dg_TERNARYsq)
      ELSE IF (L_ACDC_TERNARY) THEN
        !Use Look-Up table from Atmospheric Cluster Dynamics Code
	!which solves the birth-death equations
        CALL ACDC_TERNARY ( CS, RH, TEMP, XSA, XNH3, 
     &                     J1, SA_M_TERNARY, NH3_M_TERNARY,
     &                     Dg_TERNARYsq)
      END IF   

! Calculate H2SO4-Amine NPF
      IF (L_GLASOE) THEN
        !Use parameterization from Glasoe et al. (2014) from
	!flow-tube studies at University of Minnesota
        J2       = (1.0E6) * K_GLASOE * A1TOT**2 * DMA   ! For Glasoe et al. (2014)   J1_Glasoe
        SA_M_AMINE  = SA_M_GLASOE
        DMA_M_AMINE = DMA_M_GLASOE
        Dg_AMINEsq  = Dg_GLASOEsq
      ELSE IF (L_ACDC_DMA) THEN
        !Use Look-Up table from Atmospheric Cluster Dynamics Code
        !which solves the birth-death equations. Particles are 
        !brought into CMAQ at Dp = 1.2 nm
        CALL ACDC_AMINE ( CS, TEMP, XSA, XDMA, 
     &                    J2, SA_M_AMINE, DMA_M_AMINE,
     &                    Dg_AMINEsq)
       END IF
        
! Set up for calculation of NPF rate from other pathways
      !SAORGX   = A1DIFF   * MBTCA  !For Paasonen et al. (2010) J2
      SAORGX   = A1DIFF   * ELVOC  !For Paasonen et al. (2010) J2
      ORGYORGY = ELVOC    * ELVOC  !For Paasonen et al. (2010) J3
      SA2ORGX  = A1TOT**2 * MBTCA  !For Riccobono et al. (2014) J1_Ricco

! Calculate NPF Number Concentration rates in [ N / m**3 / sec ]
      J3       = (1.0E6) * KPAAS_2  * SAORGX         ! For Paasonen et al. (2010) J3
      J4       = (1.0E6) * K_KIRKBY * ELVOC ** ( E1_KIRKBY + E2_KIRKBY / ELVOC ) ! From Kirkby et al. (2016)
      J_RICCO  = (1.0E6) * K_RICCO  * SA2ORGX        ! For Riccobono et al. (2014) J1_Ricco
      J_PAAS3  = (1.0E6) * KPAAS_3  * ORGYORGY       ! For Paasonen et al. (2010) J4       

! Calculate Ion-Mediated NPF Rates
      CALL ION_MEDIATED_NPF( A1DIFF, ELVOC, PRESS, TEMP, RH, ASURF, 
     &                       J5, J6, SA_M_IMN, ELVOC_M_IMN, 
     &                       Dg_IMNSO4sq, Dg_IMNORGsq )

! Calculate Sulfate Mass Change to Each Chemical Species [ ug / m3 / s ]
      DMDT_SO4_BINARY  = J0 * SA_M_BINARY
      DMDT_SO4_TERNARY = J1 * SA_M_TERNARY
      DMDT_SO4_AMINE   = J2 * SA_M_AMINE 
      DMDT_SO4_PAAS2   = J3 * SA_M_PAAS2     ![ug/m3/s]
      DMDT_SO4_IMN     = J5 * SA_M_IMN       ![ug/m3/s]
      
! Compare Total SO4 Removal with SO4 Production Rate and Normalize NPF
! Rates if Necessary
      DMDT_SO4_TOT = DMDT_SO4_BINARY + DMDT_SO4_TERNARY + 
     &               DMDT_SO4_AMINE  + DMDT_SO4_PAAS2   +
     &               DMDT_SO4_IMN

      IF ( DMDT_SO4_TOT .GT. SO4RATE ) THEN
         NORM_SO4 = SO4RATE / DMDT_SO4_TOT
         DMDT_SO4_BINARY  = DMDT_SO4_BINARY  * NORM_SO4
         DMDT_SO4_TERNARY = DMDT_SO4_TERNARY * NORM_SO4
         DMDT_SO4_AMINE   = DMDT_SO4_AMINE   * NORM_SO4
         DMDT_SO4_PAAS2   = DMDT_SO4_PAAS2   * NORM_SO4
         DMDT_SO4_IMN     = DMDT_SO4_IMN     * NORM_SO4
         DMDT_SO4_TOT = SO4RATE

         J0 = DMDT_SO4_BINARY  / SA_M_BINARY
         J1 = DMDT_SO4_TERNARY / SA_M_TERNARY
         J2 = DMDT_SO4_AMINE   / SA_M_AMINE
         J3 = DMDT_SO4_PAAS2   / SA_M_PAAS2
         J5 = DMDT_SO4_IMN     / SA_M_IMN
      END IF
 
! Calculate the Mass Change for the other Species Involved
      DMDT_NH3_TOT   = J1 * NH3_M_TERNARY  ![ug/m3/s]
      DMDT_DMA_TOT   = J2 * DMA_M_AMINE    ![ug/m3/s]
      !DMDT_MBTCA_TOT = J3 * MBTCA_M_PAAS2 ![ug/m3/s]
      DMDT_ELVOC_TOT = J3 * ELVOC_M_PAAS2 + J4 * ELVOC_M_ORG 
     &               + J6 * ELVOC_M_IMN    ![ug/m3/s]

! Calculate combined rate
      JRATE    = J0 + J1 + J2  + J3 + J4 + J5 + J6 ![ N/m**3/s ]

! Calculate NPF Second Moment Rates
      M2_BINARY  = J0 * Dg_BINARYsq     ![m2/m3/s]
      M2_TERNARY = J1 * Dg_TERNARYsq    ![m2/m3/s]
      M2_AMINE   = J2 * Dg_AMINEsq      ![m2/m3/s]
      M2_PAAS2   = J3 * Dg_PAAS2sq      ![m2/m3/s]
      M2_ORG     = J4 * Dg_ORGsq        ![m2/m3/s]
      M2_IMN_SO4 = J4 * Dg_IMNSO4sq     ![m2/m3/s]
      M2_IMN_ORG = J4 * Dg_IMNORGsq     ![m2/m3/s]
      M2_RICCO   = J_RICCO * Dg_RICCOsq ![m2/m3/s]
      M2_PAAS3   = J_PAAS3 * Dg_PAAS3sq ![m2/m3/s]

      M2RATE   = M2_BINARY + M2_TERNARY 
     &         + M2_AMINE  + M2_PAAS2 + M2_ORG
     &         + M2_IMN_SO4 + M2_IMN_ORG

! Load Diagnostic Values into Multi-Element Vectors
      JDIAG(1) = J0 * 1.0E-6  !N cm-3 s-1
      JDIAG(2) = J1 * 1.0E-6  !N cm-3 s-1
      JDIAG(3) = J2 * 1.0E-6  !N cm-3 s-1
      JDIAG(4) = J3 * 1.0E-6  !N cm-3 s-1
      JDIAG(5) = J4 * 1.0E-6  !N cm-3 s-1
      JDIAG(6) = J4 * 1.0E-6  !N cm-3 s-1
      JDIAG(7) = J4 * 1.0E-6  !N cm-3 s-1
      JDIAG(8) = J_RICCO * 1.0E-6  !N cm-3 s-1
      JDIAG(9) = J_PAAS3 * 1.0E-6  !N cm-3 s-1
      JM2DIAG(1) = M2_BINARY
      JM2DIAG(2) = M2_TERNARY
      JM2DIAG(3) = M2_AMINE
      JM2DIAG(4) = M2_PAAS2
      JM2DIAG(5) = M2_ORG
      JM2DIAG(6) = M2_IMN_SO4
      JM2DIAG(7) = M2_IMN_ORG
      JM2DIAG(8) = M2_RICCO
      JM2DIAG(9) = M2_PAAS3
      JDGASDIAG(1) = DMDT_SO4_BINARY
      JDGASDIAG(2) = DMDT_SO4_TERNARY
      JDGASDIAG(3) = DMDT_SO4_AMINE
      JDGASDIAG(4) = DMDT_SO4_PAAS2
      JDGASDIAG(5) = DMDT_NH3_TOT
      JDGASDIAG(6) = DMDT_DMA_TOT
      JDGASDIAG(7) = DMDT_MBTCA_TOT
      JDGASDIAG(8) = DMDT_ELVOC_TOT

       
      RETURN
      END SUBROUTINE NPF_MULTI

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE NPF_BINARY ( RH, T, XH2SO4, DNDT, sulfmass, m2_20 )

      USE AERO_DATA
      USE AEROMET_DATA   ! Includes CONST.EXT
      USE PRECURSOR_DATA, ONLY: PRECURSOR_MW, SULF_IDX, NH3_IDX 
      USE UTILIO_DEFN
      
      IMPLICIT NONE

C   REVISION HISTORY:
C     Replacement of Kulmala et al., 1998 nucleation scheme with scheme 
C     of Vehkamaki et al. (2002) by G. Sarwar and K. Fahey - 03/2014
C

C.. References:
C     Vehkamaki, H., Kulmala, M, Napari, I., Lehtinen, K.E.J., Timmreck, C.,
C     Noppel, M., and A. Laaksonen. (2002) An improved parameterization for 
C     sulfuric acid-water nucleation rates for tropospheric and stratospheric
C     conditions.  JGR, v107(D22).    
   
C..Inputs: 
      REAL,    INTENT(IN ) :: RH              ! fractional relative humidity      
      REAL,    INTENT(IN)  :: T               ! Ambient temperature [ K ]
      REAL,    INTENT(IN)  :: XH2SO4          ! sulfuric acid concentration [ ug/m**3 ]

C.. Outputs:
      REAL(8), INTENT( OUT ) :: DNDT       ! particle number production rate [ m^-3/s ]
      REAL(8), INTENT( OUT ) :: sulfmass   ! mass of a new particle [ug]
      REAL(8), INTENT( OUT ) :: m2_20      !new-particle diameter squared [m**2]

C.. Parameters
      CHARACTER( 16 ), PARAMETER :: PNAME = 'NPF_BIN'                 

C.. Particle size parameters:
      REAL(8), PARAMETER :: d20 = 3.6E-07           ! diameter of a new particle [cm]
      REAL(8), PARAMETER :: d20sq = d20 * d20       ! new-particle diameter squared [cm**2]
      REAL(8), PARAMETER :: d20cu = d20 * d20sq     ! diameter of a new particle cubed [cm**3]
      REAL(8), PARAMETER :: sg  = 1.28              ! Standard Dev of new particle distirbution
      
      REAL(8)       :: sulfmass1                    ! inverse of sulfmass [ug**-1]

C.. Set constants and local variables for Vehkamaki et al. (2002) scheme
      REAL(8), PARAMETER :: C1 = 0.740997
      REAL(8), PARAMETER :: C2 = -0.00266379
      REAL(8), PARAMETER :: C3 = -0.00349998
      REAL(8), PARAMETER :: C4 = 0.0000504022
      REAL(8), PARAMETER :: C5 = 0.00201048
      REAL(8), PARAMETER :: C6 = -0.000183289
      REAL(8), PARAMETER :: C7 = 0.00157407
      REAL(8), PARAMETER :: C8 = -0.0000179059
      REAL(8), PARAMETER :: C9 = 0.000184403
      REAL(8), PARAMETER :: C10 = -1.50345E-6

      REAL(8) :: XSTAR           ! mole fraction of sulfuric acid in the critical cluster
      REAL(8) :: XMASS           ! mass fraction of sulfuric acid in the critical cluster
      REAL(8) :: NA              ! total gas phase concentration of H2SO4  [ #/cm**3 ]
      REAL(8) :: TEMP            ! ambient temperature
      REAL(8) :: LNRH, LNNA      ! LN(RH), LN(Na)
      REAL(8) :: LNRH2, LNNA2    ! LN(RH)**2, LN(Na)**2
      REAL(8) :: LNRH3, LNNA3    ! LN(RH)**3, LN(Na)**3  
      REAL(8) :: TEMP2, TEMP3    ! TEMP**2, TEMP**3
      REAL(8) :: XFAC       ! exponential term for the nucleation rate  
      REAL(8) :: Jnuc       ! nucleation rate [ #/cm**3 s ]
      REAL(8) :: m3         ! 3rd moment of 1 particle [ cm3/particle ]
      REAL(8) :: R0         ! particle density [ kg/m3 ]
      
      REAL(8) :: A, B, C, D, E, F, G, H, I, J

      REAL(8) :: MW_H2SO4        ! MW of H2SO4 in [ g / mole ]        
      REAL(8) :: DENSITY_H2SO4   ! DENSITY of H2SO4 in [ kg / m**3 ]
      REAL(8) :: DENSITY_H2O     ! DENSITY of H2O   in [ kg / m**3 ]

      REAL(8), PARAMETER :: SCALEFAC = 1.0E-06           ! for [ 1 / m**3 ] ->  [ 1 / cm**3 ]
      REAL(8), PARAMETER :: MUG2G = 1.0E-6               ! [ ug  ] -> [ g ]                 

C.. Initialize variables
      DNDT   = 0.0
      m2_20  = 1.0E-4 * d20sq       ! new-particle diameter squared [m**2]
      TEMP = T

C.. Calculate molecular weight of H2SO4 [ g / mole ] 
      MW_H2SO4 = PRECURSOR_MW( SULF_IDX )  
         
C.. Calculate density of sulfuric acid [ kg / m**3 ] 
      DENSITY_H2SO4 = AEROSPC( ASO4_IDX )%DENSITY           
      DENSITY_H2O   = 1000.0           

C.. Calculate mass of a new particle [m3/particle]
      m3  = d20cu  * exp( 9./2. * ( log( sg  )**2 ) )

C.. Calculate inverse of sulfmass [ug**-1]      
      !sulfmass1 = 1.0 / sulfmass   
     
C.. Calculate sulfuric acid concentration in molecules/cm3 
      NA = XH2SO4 * MUG2G * AVO * SCALEFAC / MW_H2SO4     

C.. The parameterization is valid at sulfuric acid concentrations of 1.0E4 - 1.0E11 molecules cm-3
      NA  = MAX (NA, 1.0E4)
      NA  = MIN (1.0E11, NA)

C.. The parameterization is valid at temperatures of 190.00-305.15 K
      TEMP = MAX (TEMP, 190.00)
      TEMP = MIN (305.15, TEMP)

C.. The parameterization is valid at RH of 0.0001-1.0
C   aero_driver.f limits RH to 0.005-0.99; thus no additional constraint is needed

C.. Define convenient constants
      TEMP2 = TEMP * TEMP
      TEMP3 = TEMP * TEMP2
      
      LNRH = LOG( RH )  
      LNNA = LOG( NA )
      LNRH2 = LNRH * LNRH
      LNRH3 = LNRH * LNRH2
      LNNA2 = LNNA * LNNA
      LNNA3 = LNNA * LNNA2

C.. Calculate mole fraction of sulfuric acid in the critical cluster
      XSTAR = C1 + C2 * TEMP + C3 * LNNA + C4 * TEMP * LNNA + 
     &      C5 * LNRH + C6 * TEMP * LNRH +
     &      C7 * LNRH2 + C8 * TEMP * LNRH2 + C9 * LNRH3 + 
     &      C10 * TEMP * LNRH3

C.. Calculate coefficients needed for the nucleation rate [Eq-12 - Vehkamaki et al., 2002]           
      A = 0.14309 + 2.21956 * TEMP - 0.0273911 * TEMP2 + 
     &      0.0000722811 * TEMP3 + 5.91822 / XSTAR
     
      B = 0.117489 + 0.462532 * TEMP - 0.0118059 * TEMP2 + 
     &      0.0000404196 * TEMP3 + 15.7963 / XSTAR
     
      C = -0.215554 - 0.0810269 * TEMP + 0.00143581 * TEMP2 - 
     &      4.7758E-6 * TEMP3 - 2.91297 / XSTAR
     
      D = -3.58856 + 0.049508 * TEMP - 0.00021382 * TEMP2 + 
     &      3.10801E-7 * TEMP3 - 0.0293333 / XSTAR
     
      E = 1.14598 - 0.600796 * TEMP + 0.00864245 * TEMP2 - 
     &      0.0000228947 * TEMP3 - 8.44985 / XSTAR
     
      F = 2.15855 + 0.0808121 * TEMP - 0.000407382 * TEMP2 - 
     &      4.01957E-7 * TEMP3 + 0.721326 / XSTAR

      G = 1.6241 - 0.0160106 * TEMP + 0.0000377124 * TEMP2 + 
     &      3.21794E-8 * TEMP3 - 0.0113255 / XSTAR
      
      H = 9.71682 - 0.115048 * TEMP + 0.000157098 * TEMP2 + 
     &      4.00914E-7 * TEMP3 + 0.71186 / XSTAR
     
      I = -1.05611 + 0.00903378 * TEMP - 0.0000198417 * TEMP2 +
     &      2.46048E-8 * TEMP3 - 0.0579087 / XSTAR

      J = -0.148712 + 0.00283508 * TEMP - 9.24619E-6 * TEMP2 + 
     &      5.00427E-9 * TEMP3 - 0.0127081 / XSTAR

C.. Calculate the exponential term for the nucleation rate [Eq-12 - Vehkamaki et al., 2002] 
      XFAC = A + B * LNRH + C * LNRH2 + D * LNRH3 + 
     &      E * LNNA + F * LNRH * LNNA + G * LNRH2 * LNNA +     
     &      H * LNNA2 + I * LNRH * LNNA2 + J * LNNA3

C.. Calculate particle nucleation rate: unit [ 1 / cm**3 s] [Eq-12 - Vehkamaki et al., 2002] 
      Jnuc = EXP(XFAC)  
      
C.. The parameterization is valid for nucleation rates of 1.0E-7-1.0E10 [ 1 / cm**3 s]
      Jnuc = MAX (Jnuc, 1.0E-7)
      Jnuc = MIN (1.0E10, Jnuc)

C.. Convert the unit of particle nucleation rate into [ 1 / m**3 s] by multiplying it by 1.0E6  
      DNDT = Jnuc * 1.0E06  ! (1/(m**3 s))

C.. Calculate the Mass of Sulfuric Acid in the Particle
      XMASS = XSTAR * MW_H2SO4 / ( XSTAR*MW_H2SO4 + (1-XSTAR)*18.0 )
      R0    = XMASS * DENSITY_H2SO4 + (1-XMASS) * DENSITY_H2O !Particle Density [kg/m3]
      sulfmass = 1.0E+9 * R0 * PI/6. * m3 * XMASS  !Particle Mass[ug/particle]

      RETURN    
      END SUBROUTINE NPF_BINARY 
 
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE NPF_NAPARI ( PRESS, RHin, T, XH2SO4, XNH3, DNDT, 
     &                         sulfmass, nh4mass, m2_20 )

      USE AERO_DATA
      USE AEROMET_DATA   ! Includes CONST.EXT
      USE PRECURSOR_DATA, ONLY: PRECURSOR_MW, SULF_IDX, NH3_IDX 
      USE UTILIO_DEFN
      
      IMPLICIT NONE

C   REVISION HISTORY:
C     Added by Ben Murphy to capture New Particle Formation via
C     H2SO4-NH3-H2O interactions

C.. References:
C     Napari, I., M. Noppel, H. Vehkamaki, and M. Kulmala.
C     "Parametrization of
C     Ternary Nucleation Rates for H2so4-Nh3-H2o Vapors." Journal of
C     Geophysical Research-Atmospheres 107, no. D19 (2002).
   
C..Inputs: 
      REAL,    INTENT(IN ) :: PRESS           ! Air Pressure [ Pa ]
      REAL,    INTENT(IN ) :: RHin            ! fractional relative humidity      
      REAL,    INTENT(IN)  :: T               ! Ambient temperature [ K ]
      REAL,    INTENT(IN)  :: XH2SO4          ! sulfuric acid concentration [ ug/m**3 ]
      REAL,    INTENT(IN)  :: XNH3            ! ammonia concentration [ ug/m**3 ]

C.. Outputs:
      REAL(8), INTENT( OUT ) :: DNDT       ! particle number production rate [ m^-3/s ]
      REAL(8), INTENT( OUT ) :: sulfmass   ! mass of sulfuric acid in new particle [ug]
      REAL(8), INTENT( OUT ) :: nh4mass    ! mass of ammonium in new particle [ug]
      REAL(8), INTENT( OUT ) :: m2_20      !new-particle diameter squared [m**2]

C.. Parameters
      CHARACTER( 16 ), PARAMETER :: PNAME = 'NPF_TERN'                 

C.. Particle size parameters:
      REAL(8) :: RNUC           ! radius of a new particle [nm]
      REAL(8) :: FNL            ! Natural Logarithm of Formation Rate
      REAL(8) :: d20            ! diameter of a new particle [cm]
      REAL(8) :: sg  = 1.28     ! Standard Dev of new particle distirbution
      
      REAL(8) :: TEMP       ! ambient temperature
      REAL(8) :: XFAC       ! exponential term for the nucleation rate  
      REAL(8) :: Jnuc       ! nucleation rate [ #/cm**3 s ]
      REAL(8) :: m3         ! 3rd moment of 1 particle [ cm3/particle ]
      REAL(8) :: H2SO4      ! Molecule Concentration of Sulfuric Acid [molec cm-3]
      REAL(8) :: NH3        ! Mixing Ration of NH3 [PPTV]
      REAL(8) :: RH         ! Internal Relative Humidity
      REAL(8) :: lnH2SO4    ! Nat Log of H2SO4 Concentration
      REAL(8) :: lnNH3      ! Nat Log of NH3 Mixing Ratio
      REAL(8) :: lnRH       ! Nat Log of RH
      
      REAL(8) :: MW_H2SO4        ! MW of H2SO4 in [ g / mole ]        
      REAL(8) :: MW_NH3          ! MW of NH3 in [ g / mole ]        
      REAL(8) :: DENSITY_H2SO4   ! DENSITY of H2SO4 in [ kg / m**3 ]
      REAL(8) :: n_H2SO4         ! Number of H2SO4 molecules in Cluster
      REAL(8) :: n_NH4           ! Number of Ammonia molecules in Cluster

      REAL(8), PARAMETER :: SCALEFAC = 1.0E-06    ! for [ 1 / m**3 ] ->  [ 1 / cm**3 ]
      REAL(8), PARAMETER :: MUG2G = 1.0E-6        ! [ ug  ] -> [ g ]                 
      REAL(8), PARAMETER :: R = 8.314             ! m3 Pa K-1 mol-1                 

      INTEGER :: I

C.. Set Parameters for Napari et al. (2002) Scheme
      REAL(8) :: aa0(20), a1(20), a2(20), a3(20), fa(20) !Parameters

      data AA0 /-0.355297, 3.13735,    19.0359,   1.07605, 6.0916,
     &          0.31176,   -0.0200738, 0.165536,
     &          6.52645,   3.68024,    -0.066514, 0.65874,
     &          0.0599321, -0.732731,  0.728429,  41.3016,
     &          -0.160336, 8.57868,    0.0530167, -2.32736        /

      data A1 /-33.8449, -0.772861, -0.170957, 1.48932, -1.25378,
     &         1.64009, -0.752115, 3.26623, -0.258002, -0.204098,
     &         -7.82382, 0.190542, 5.96475, -0.0184179, 3.64736,
     &         -0.35752, 0.00889881, -0.112358, -1.98815, 0.0234646/

      data A2 /0.34536, 0.00561204, 0.000479808, -0.00796052,
     &         0.00939836, -0.00343852, 0.00525813, -0.0489703,
     &         0.00143456, 0.00106259, 0.0122938, -0.00165718,
     &         -0.0362432, 0.000147186, -0.027422, 0.000904383,
     &         -5.39514d-05, 0.000472626, 0.0157827, -0.000076519/

      data A3 /-0.000824007, -9.74576d-06, -4.14699d-07, 7.61229d-06,
     &         -1.74927d-05, -1.09753d-05, -8.98038d-06, 0.000146967,
     &         -2.02036d-06, -1.2656d-06, 6.18554d-05, 3.41744d-06,
     &         4.93337d-05, -2.37711d-07, 4.93478d-05, -5.73788d-07,
     &         8.39522d-08, -6.48365d-07, -2.93564d-05, 8.0459d-08   /


C.. Initialize variables
      DNDT   = 0.0
      RH = RHin
      TEMP = T

C.. Calculate molecular weight of precursors [ g / mole ] 
      MW_H2SO4 = PRECURSOR_MW( SULF_IDX )  
      MW_NH3   = PRECURSOR_MW( NH3_IDX )  
         
C.. Calculate sulfuric acid concentration in molecules/cm3 
      H2SO4 = XH2SO4 * MUG2G * AVO * SCALEFAC / MW_H2SO4

C.. Calculate Ammonia Concentration [ug/m3 --> PPTV]
      NH3   = XNH3   * MUG2G / MW_NH3 * R*TEMP/PRESS * 1.0e+12

C.. The parameterization is valid at sulfuric acid concentrations of 1.0E4 - 1.0E11 molecules cm-3
      H2SO4  = MAX (H2SO4, 1.0E4)
      H2SO4  = MIN (1.0E9, H2SO4)

C.. The parameterization is valid at ammonia mixing ratios of 0.1 - 100 pptv
      NH3  = MAX (NH3, 0.1)
      NH3  = MIN (100.0, NH3)

C.. The parameterization is valid at temperatures of 190.00-305.15 K
      TEMP = MAX (TEMP, 240.00)
      TEMP = MIN (300.00, TEMP)

C.. The parameterization is valid at RH of 0.05-0.95 K
      RH = MAX (RH, 0.05)
      RH = MIN (0.95, RH)

C.. Calculate New Particle Formation Rate
      DO I = 1,20
         FA(I) = AA0(I) + A1(I)*TEMP + A2(I)*TEMP**2. + A3(I)*TEMP**3.
      ENDDO
 
      lnH2SO4 = log(H2SO4)
      lnNH3   = log(NH3)
      lnRH    = log(RH)
     
      FNL = -84.7551                 + FA(1)/lnH2SO4 
     &      + FA(2)*lnH2SO4          + FA(3)*lnH2SO4**2
     &      + FA(4)*lnNH3            + FA(5)*lnNH3**2 
     &      + FA(6)*RH               + FA(7)*lnRH  
     &      + FA(8)*lnNH3/lnH2SO4    + FA(9)*lnNH3*lnH2SO4
     &      + FA(10)*RH*lnH2SO4      + FA(11)*RH/lnH2SO4 
     &      + FA(12)*RH*lnNH3        + FA(13)*lnRH/lnH2SO4 
     &      + FA(14)*lnRH*lnNH3      + FA(15)*lnNH3**2/lnH2SO4 
     &      + FA(16)*lnH2SO4*lnNH3**2+ FA(17)*lnH2SO4**2*lnNH3
     &      + FA(18)*RH*lnNH3**2     + FA(19)*RH*lnNH3/lnH2SO4 
     &      + FA(20)*lnH2SO4**2*lnNH3**2

C.. Calculate particle nucleation rate: unit [ 1 / cm**3 s] [Eq-12 - Vehkamaki et al., 2002] 
      Jnuc = EXP(FNL)
      Jnuc = Jnuc * 1.e-6  !Scaling Parameter from European Regional Model Runs
      
C.. The parameterization is valid for nucleation rates of 1.0E-7-1.0E10 [ 1 / cm**3 s]
      Jnuc = MAX (Jnuc, 1.0E-7)
      Jnuc = MIN (1.0E6, Jnuc)
      FNL = log(Jnuc)

C.. Calculate critical Radius in Nanometers
      rnuc = 0.141027-0.00122625*FNL-7.82211d-6*FNL**2
     &      -0.00156727*TEMP-0.00003076*temp*FNL
     &      +0.0000108375*TEMP**2
      rnuc = MAX (rnuc, 0.405)

C.. Calculate Composition of New Particles (molecules/particle)
      n_H2SO4 = 38.1645 + 0.774106*FNL + 0.00298879*FNL**2
     &          -0.357605*TEMP - 0.00366358*TEMP*FNL + 0.0008553*TEMP**2
      n_NH4   = 26.8982 + 0.682905*FNL + 0.00357521*FNL**2
     &          -0.265748*TEMP - 0.00341895*TEMP*FNL + 0.000673454*TEMP**2

C.. Convert the unit of particle nucleation rate into [ 1 / m**3 s] by multiplying it by 1.0E6  
      DNDT = Jnuc * 1.0E06  ! (1/(m**3 s))

C.. Calculate Size of New Particles
      d20 = rnuc * 2.           !New Particle Diamter [nm]
      m2_20  = 1.0E-18 * d20**2 ! new-particle diameter squared [m**2]

C.. Calculate Mass of Each Species In Particle (g/particle)
    !Use the molecular weight of the gas-phase species since
    !that is the actual mass that has condensed. Convert that
    !mass to the mass of the ions before adding to the aerosol
    !master variable in VOLINORG
      sulfmass = n_H2SO4 / AVO * MW_H2SO4 * 1.0e6
      nh4mass  = n_NH4   / AVO * MW_NH3   * 1.0e6

      RETURN    
      END SUBROUTINE NPF_NAPARI 
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ACDC_TERNARY ( CSin, RHin, TEMPin, XH2SO4, XNH3, DNDT, 
     &                         sulfmass, nh3mass, m2_20 )

      USE AERO_DATA
      USE AEROMET_DATA   ! Includes CONST.EXT
      USE PRECURSOR_DATA, ONLY: PRECURSOR_MW, SULF_IDX, NH3_IDX 
      
      IMPLICIT NONE

C   REVISION HISTORY:
C     Added by Ben Murphy to capture New Particle Formation via
C     H2SO4-NH3-H2O interactions. These rates are brought in from
C     a look-up table informed by the Atmospheric Cluster Dynamics
C     Code developed at the University of Helsinki.

C.. References:
C

C..Inputs: 
      REAL,    INTENT(IN ) :: CSin            ! Condesation Sink (s-1)      
      REAL,    INTENT(IN ) :: RHin            ! fractional relative humidity      
      REAL,    INTENT(IN)  :: TEMPin          ! Ambient temperature [ K ]
      REAL,    INTENT(IN)  :: XH2SO4          ! sulfuric acid concentration [ ug/m**3 ]
      REAL,    INTENT(IN)  :: XNH3            ! ammonia concentration [ ug/m**3 ]

C.. Outputs:
      REAL(8), INTENT( OUT ) :: DNDT       ! particle number production rate [ m-3 s-1 ]
      REAL(8), INTENT( OUT ) :: sulfmass   ! mass of sulfuric acid used to make new particle [ug]
      REAL(8), INTENT( OUT ) :: nh3mass    ! mass of ammonia used to make new particle [ug]
      REAL(8), INTENT( OUT ) :: m2_20      ! new-particle diameter squared [m**2]

C.. Parameters
      CHARACTER( 16 ), PARAMETER :: PNAME = 'NPF_ACDC_TERN'                 

C.. Particle size parameters:
      REAL :: CS              ! Condesation Sink (s-1)      
      REAL :: RH              ! fractional relative humidity      
      REAL :: TEMP            ! Ambient temperature [ K ]
      REAL :: H2SO4      ! Molecule Concentration of Sulfuric Acid [molec cm-3]
      REAL :: NH3        ! Mixing Ration of NH3 [PPTV]
      REAL :: JNUC       ! nucleation rate [particles cm-3 s-1]
      
      REAL :: MW_H2SO4        ! MW of H2SO4 in [ g / mole ]        
      REAL :: MW_NH3          ! MW of NH3 in [ g / mole ]       
      REAL, PARAMETER :: R = 8.314             ! m3 Pa K-1 mol-1                 
      REAL, PARAMETER :: SCALEFAC = 1.0E-06    ! for [ 1 / m**3 ] ->  [ 1 / cm**3 ]
      REAL, PARAMETER :: MUG2G = 1.0E-6        ! [ ug  ] -> [ g ]                 

C.. Assert Size of New Particles
      REAL, PARAMETER :: RNUC = 0.6      ! radius of a new particle [nm]
      REAL, PARAMETER :: d20 = RNUC * 2. ! New Particle Diamter [nm]

C.. Calculate Mass of Each Vapor Species Used For Particle (g/particle)
    !First compute the mass of the condensed species in the particle
    !for its size, using bulk density. Then convert to the mass of the
    !vapor species needed for that condensation. Assume Ammonium
    !Bisulfate (NH4-HSO4) for these particles.                                                              
      REAL, PARAMETER :: p_dens = 1.78E12   !NH4-HSO4 Density (ug m-3)                                      
      REAL, PARAMETER :: p_mass = 4./3.*3.1415 * (RNUC*1.0E-9)**3 * p_dens  !Particle Mass (ug)             
      REAL, PARAMETER :: sulf_mfrac = 0.8144  !Mass Fraction of Bisulfate                                   
      REAL, PARAMETER :: nh4_mfrac  = 0.1856  !Mass Fraction of Ammonium                                    
      REAL, PARAMETER :: MW_HSO4 = 97.06  !g mol-1                                                          
      REAL, PARAMETER :: MW_NH4  = 18.03  !g mol-1                                                          
      
C.. Temporary Variables for Interpolating LookUp Table to Model Conditions
      INTEGER :: IAC, INH, ITK, IRH, ICS
      REAL :: AC1, NH1, TK1, RH1, CS1
      REAL :: PAC, PNH, PTK, PRH, PCS
      REAL :: JAC, JNH, JTK, JRH, JCS
      REAL :: QAC, QNH, QTK, QRH, QCS

C.. Initialize variables                                                                                    
      DNDT   = 0.0                                                                                          
                                                                                                            
C.. Fetch molecular weight of precursors [ g / mole ]                                                       
      MW_H2SO4 = PRECURSOR_MW( SULF_IDX )                                                                   
      MW_NH3   = PRECURSOR_MW( NH3_IDX )                                                                    
                                                                                                            
C.. Calculate Mass of Precursors Used to Make One Particle                                                  
      sulfmass = sulf_mfrac * p_mass / MW_HSO4 * MW_H2SO4  !H2SO4 Mass (ug m-3 s-1)                         
      nh3mass  = nh4_mfrac  * p_mass / MW_NH4  * MW_NH3    !NH3 Mass (ug m-3 s-1                            
      m2_20  = 1.0D-18 * d20**2 ! new-particle diameter squared [m**2]
                                                                                                            
C.. Calculate sulfuric acid concentration in molecules cm-3                                                  
      H2SO4 = XH2SO4 * MUG2G * AVO * SCALEFAC / MW_H2SO4                                                    
                                                                                                            
C.. Calculate Ammonia Concentration [ug/m3 --> molec cm-3]                                                        
      NH3   = XNH3   * MUG2G *AVO * SCALEFAC / MW_NH3                                               
                                                                                                            
C.. The parameterization is valid at sulfuric acid conc. of 1.0E4 - 3.16E9 molec cm-3          
      H2SO4  = MAX (H2SO4, 1.0E4)                                                                           
      H2SO4  = MIN (3.162E9, H2SO4)                                                                           
                                                                                                            
C.. The parameterization is valid at ammonia conc. of 1.00e6 - 1.00E11 molec cm-3
      NH3  = MAX (NH3, 1.00E6)
      NH3  = MIN (1.00E11, NH3)

C.. The parameterization is valid at temperatures of 180.00-320.0 K
      TEMP = MAX (TEMPin, 180.00)
      TEMP = MIN (320.00, TEMP)

C.. The parameterization is valid at RH of 0.0-100.0
      RH = RHin * 100.0   !Convert to Percent RH
      RH = MAX (RH, 0.0)
      RH = MIN (100.0, RH)

C.. The parameterization is valid at Condensation Sink of 1.0E-5 to 1.0E-1 s-1
      CS = MAX (CSin, 1.0E-5)
      CS = MIN (1.0E-1, CS)

C.. Normalize Input Variables to LookUp Table Space
      AC1 = ( LOG10( H2SO4 ) - AC0 ) / DAC
      NH1 = ( LOG10( NH3   ) - NH0 ) / DNH
      TK1 = ( TEMP           - TK0 ) / DTK
      RH1 = ( RH             - RH0 ) / SRH
      CS1 = ( LOG10( CS    ) - CS0 ) / DCS_T

C.. Calculate Temporary Interpolating Variables
      IAC = 1 + INT( AC1 )
      INH = 1 + INT( NH1 )
      ITK = 1 + INT( TK1 )
      IRH = 1 + INT( RH1 )
      ICS = 1 + INT( CS1 )

      PAC = FLOAT( IAC ) - AC1
      PNH = FLOAT( INH ) - NH1
      PTK = FLOAT( ITK ) - TK1
      PRH = FLOAT( IRH ) - RH1
      PCS = FLOAT( ICS ) - CS1

      JAC = MIN( IAC + 1, SIZE( ACDC_Tern,1 ) )
      JNH = MIN( INH + 1, SIZE( ACDC_Tern,2 ) )
      JTK = MIN( ITK + 1, SIZE( ACDC_Tern,3 ) )
      JRH = MIN( IRH + 1, SIZE( ACDC_Tern,4 ) )
      JCS = MIN( ICS + 1, SIZE( ACDC_Tern,5 ) )

      QAC = 1.0 - PAC
      QNH = 1.0 - PNH
      QTK = 1.0 - PTK
      QRH = 1.0 - PRH
      QCS = 1.0 - PCS

C.. Calculate New Particle Formation Rate
      JNUC = PAC * PNH * PTK * PRH * PCS * ACDC_Tern( IAC, INH, ITK, IRH, ICS ) +   
     &       PAC * PNH * PTK * PRH * QCS * ACDC_Tern( IAC, INH, ITK, IRH, JCS ) +   
     &       PAC * PNH * PTK * QRH * PCS * ACDC_Tern( IAC, INH, ITK, JRH, ICS ) +   
     &       PAC * PNH * PTK * QRH * QCS * ACDC_Tern( IAC, INH, ITK, JRH, JCS ) +   
     &       PAC * PNH * QTK * PRH * PCS * ACDC_Tern( IAC, INH, JTK, IRH, ICS ) +   
     &       PAC * PNH * QTK * PRH * QCS * ACDC_Tern( IAC, INH, JTK, IRH, JCS ) +   
     &       PAC * PNH * QTK * QRH * PCS * ACDC_Tern( IAC, INH, JTK, JRH, ICS ) +   
     &       PAC * PNH * QTK * QRH * QCS * ACDC_Tern( IAC, INH, JTK, JRH, JCS ) +   
     &       PAC * QNH * PTK * PRH * PCS * ACDC_Tern( IAC, JNH, ITK, IRH, ICS ) +   
     &       PAC * QNH * PTK * PRH * QCS * ACDC_Tern( IAC, JNH, ITK, IRH, JCS ) +   
     &       PAC * QNH * PTK * QRH * PCS * ACDC_Tern( IAC, JNH, ITK, JRH, ICS ) +   
     &       PAC * QNH * PTK * QRH * QCS * ACDC_Tern( IAC, JNH, ITK, JRH, JCS ) +   
     &       PAC * QNH * QTK * PRH * PCS * ACDC_Tern( IAC, JNH, JTK, IRH, ICS ) +   
     &       PAC * QNH * QTK * PRH * QCS * ACDC_Tern( IAC, JNH, JTK, IRH, JCS ) +   
     &       PAC * QNH * QTK * QRH * PCS * ACDC_Tern( IAC, JNH, JTK, JRH, ICS ) +   
     &       PAC * QNH * QTK * QRH * QCS * ACDC_Tern( IAC, JNH, JTK, JRH, JCS ) +   
     &       QAC * PNH * PTK * PRH * PCS * ACDC_Tern( JAC, INH, ITK, IRH, ICS ) +   
     &       QAC * PNH * PTK * PRH * QCS * ACDC_Tern( JAC, INH, ITK, IRH, JCS ) +   
     &       QAC * PNH * PTK * QRH * PCS * ACDC_Tern( JAC, INH, ITK, JRH, ICS ) +   
     &       QAC * PNH * PTK * QRH * QCS * ACDC_Tern( JAC, INH, ITK, JRH, JCS ) +   
     &       QAC * PNH * QTK * PRH * PCS * ACDC_Tern( JAC, INH, JTK, IRH, ICS ) +   
     &       QAC * PNH * QTK * PRH * QCS * ACDC_Tern( JAC, INH, JTK, IRH, JCS ) +   
     &       QAC * PNH * QTK * QRH * PCS * ACDC_Tern( JAC, INH, JTK, JRH, ICS ) +   
     &       QAC * PNH * QTK * QRH * QCS * ACDC_Tern( JAC, INH, JTK, JRH, JCS ) +   
     &       QAC * QNH * PTK * PRH * PCS * ACDC_Tern( JAC, JNH, ITK, IRH, ICS ) +   
     &       QAC * QNH * PTK * PRH * QCS * ACDC_Tern( JAC, JNH, ITK, IRH, JCS ) +   
     &       QAC * QNH * PTK * QRH * PCS * ACDC_Tern( JAC, JNH, ITK, JRH, ICS ) +   
     &       QAC * QNH * PTK * QRH * QCS * ACDC_Tern( JAC, JNH, ITK, JRH, JCS ) +   
     &       QAC * QNH * QTK * PRH * PCS * ACDC_Tern( JAC, JNH, JTK, IRH, ICS ) +   
     &       QAC * QNH * QTK * PRH * QCS * ACDC_Tern( JAC, JNH, JTK, IRH, JCS ) +   
     &       QAC * QNH * QTK * QRH * PCS * ACDC_Tern( JAC, JNH, JTK, JRH, ICS ) +   
     &       QAC * QNH * QTK * QRH * QCS * ACDC_Tern( JAC, JNH, JTK, JRH, JCS )    


C.. Convert the unit of particle nucleation rate into [ 1 / m**3 s] by multiplying it by 1.0E6  
      DNDT = DBLE(JNUC) * 1.0D06  ! (1/(m**3 s))

      RETURN    
      END SUBROUTINE ACDC_TERNARY 


C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ACDC_AMINE  ( CSin, TEMPin, XH2SO4, XDMA, DNDT, 
     &                         sulfmass, dmamass, m2_20 )

      USE AERO_DATA
      USE AEROMET_DATA   ! Includes CONST.EXT
      USE PRECURSOR_DATA, ONLY: PRECURSOR_MW, SULF_IDX, DMA_IDX 
      
      IMPLICIT NONE

C   REVISION HISTORY:
C     Added by Ben Murphy to capture New Particle Formation via
C     H2SO4-NH3-H2O interactions. These rates are brought in from
C     a look-up table informed by the Atmospheric Cluster Dynamics
C     Code developed at the University of Helsinki.

C.. References:
C

C..Inputs: 
      REAL,    INTENT(IN ) :: CSin            ! Condesation Sink (s-1)      
      REAL,    INTENT(IN)  :: TEMPin          ! Ambient temperature [ K ]
      REAL,    INTENT(IN)  :: XH2SO4          ! sulfuric acid concentration [ ug/m**3 ]
      REAL,    INTENT(IN)  :: XDMA            ! ammonia concentration [ ug/m**3 ]

C.. Outputs:
      REAL(8), INTENT( OUT ) :: DNDT       ! particle number production rate [ m-3 s-1 ]
      REAL(8), INTENT( OUT ) :: sulfmass   ! mass of sulfuric acid used to make new particle [ug]
      REAL(8), INTENT( OUT ) :: dmamass    ! mass of ammonia used to make new particle [ug]
      REAL(8), INTENT( OUT ) :: m2_20      ! new-particle diameter squared [m**2]

C.. Parameters
      CHARACTER( 16 ), PARAMETER :: PNAME = 'NPF_ACDC_AMINE'                 

C.. Particle size parameters:
      REAL :: H2SO4      ! Molecule Concentration of Sulfuric Acid [molec cm-3]
      REAL :: DMA        ! Conc. of DMA [molec cm-3]
      REAL :: TEMP       ! Temperature
      REAL :: CS         ! Internal Condensation Sink
      REAL :: JNUC       ! nucleation rate [particles cm-3 s-1]
      
      REAL :: MW_H2SO4        ! MW of H2SO4 in [ g / mole ]        
      REAL :: MW_DMA          ! MW of DMA in [ g / mole ]       
      REAL, PARAMETER :: R = 8.314             ! m3 Pa K-1 mol-1                 
      REAL, PARAMETER :: SCALEFAC = 1.0E-06    ! for [ 1 / m**3 ] ->  [ 1 / cm**3 ]
      REAL, PARAMETER :: MUG2G = 1.0E-6        ! [ ug  ] -> [ g ]                 

C.. Assert Size of New Particles
      REAL, PARAMETER :: RNUC = 0.6      ! radius of a new particle [nm]
      REAL, PARAMETER :: d20 = RNUC * 2. ! New Particle Diamter [nm]

C.. Calculate Mass of Each Vapor Species Used For Particle (g/particle)
    !First compute the mass of the condensed species in the particle
    !for its size, using bulk density. Then convert to the mass of the
    !vapor species needed for that condensation. Assume Dimethylaminium
    !Bisulfate (CH3CH3NH2-HSO4) for these particles.                                                              
      REAL, PARAMETER :: p_dens = 1.47E12   !Guess at the DMAH-HSO4 Density (ug m-3) 
      REAL, PARAMETER :: p_mass = 4./3.*3.1415 * (RNUC*1.0E-9)**3 * p_dens  !Particle Mass (ug)             

C.. Temporary Variables for Interpolating LookUp Table to Model Conditions
      INTEGER :: IAC, INH, ITK, ICS
      REAL :: AC1, NH1, TK1, CS1
      REAL :: PAC, PNH, PTK, PCS
      REAL :: JAC, JNH, JTK, JCS
      REAL :: QAC, QNH, QTK, QCS

C.. Initialize variables                                                                                    
      DNDT   = 0.0                                                                                          
      
C.. Calculate Output Particle Properties
      sulfmass = 0.685 * p_mass  !Mass of Sulfuric Acid (ug ptcl-1)
      dmamass  = 0.315 * p_mass  !Mass of DMA (ug ptcl-1)
      m2_20  = 1.0D-18 * d20**2 ! new-particle diameter squared [m**2]

C.. Fetch molecular weight of precursors [ g / mole ]                                                       
      MW_H2SO4 = PRECURSOR_MW( SULF_IDX )                                                                   
      MW_DMA = PRECURSOR_MW( DMA_IDX ) 
                                                                                                            
C.. Calculate sulfuric acid concentration in molecules cm-3                                                  
      H2SO4 = XH2SO4 * MUG2G * AVO * SCALEFAC / MW_H2SO4                                                    
                                                                                                            
C.. Calculate Ammonia Concentration [ug/m3 --> molec cm-3]                                                        
      DMA   = XDMA   * MUG2G *AVO * SCALEFAC / MW_DMA                                               
                                                                                                            
C.. The parameterization is valid at sulfuric acid conc. of 1.0E4 - 3.16E9 molec cm-3          
      H2SO4  = MAX (H2SO4, 1.0E4)                                                                           
      H2SO4  = MIN (3.162E9, H2SO4)                                                                           
                                                                                                            
C.. The parameterization is valid at DMA conc. of 1.00e4 - 1.00E9 molec cm-3
      DMA  = MAX (DMA, 1.00E4)
      DMA  = MIN (1.00E9, DMA)

C.. The parameterization is valid at temperatures of 180.00-320.0 K
      TEMP = MAX (TEMPin, 180.00)
      TEMP = MIN (320.00, TEMP)

C.. The parameterization is valid at Condensation Sink of 1.0E-5 to 1.0E-1 s-1
      CS = MAX (CSin, 1.0E-5)
      CS = MIN (1.0E-1, CS)

C.. Normalize Input Variables to LookUp Table Space
      AC1 = ( LOG10( H2SO4 ) - AC0 ) / DAC
      NH1 = ( LOG10( DMA   ) - DM0 ) / DDM
      TK1 = ( TEMP           - TK0 ) / DTK
      CS1 = ( LOG10( CS    ) - CS0 ) / DCS_A

C.. Calculate Temporary Interpolating Variables
      IAC = 1 + INT( AC1 )
      INH = 1 + INT( NH1 )
      ITK = 1 + INT( TK1 )
      ICS = 1 + INT( CS1 )

      PAC = FLOAT( IAC ) - AC1
      PNH = FLOAT( INH ) - NH1
      PTK = FLOAT( ITK ) - TK1
      PCS = FLOAT( ICS ) - CS1

      JAC = MIN( IAC + 1, SIZE( ACDC_DMA,1 ) )
      JNH = MIN( INH + 1, SIZE( ACDC_DMA,2 ) )
      JTK = MIN( ITK + 1, SIZE( ACDC_DMA,3 ) )
      JCS = MIN( ICS + 1, SIZE( ACDC_DMA,4 ) )

      QAC = 1.0 - PAC
      QNH = 1.0 - PNH
      QTK = 1.0 - PTK
      QCS = 1.0 - PCS

C.. Calculate New Particle Formation Rate (particle cm-3 s-1)
      JNUC = PAC * PNH * PTK * PCS * ACDC_DMA( IAC, INH, ITK, ICS ) +   
     &       PAC * PNH * PTK * QCS * ACDC_DMA( IAC, INH, ITK, JCS ) +   
     &       PAC * PNH * QTK * PCS * ACDC_DMA( IAC, INH, JTK, ICS ) +   
     &       PAC * PNH * QTK * QCS * ACDC_DMA( IAC, INH, JTK, JCS ) +   
     &       PAC * QNH * PTK * PCS * ACDC_DMA( IAC, JNH, ITK, ICS ) +   
     &       PAC * QNH * PTK * QCS * ACDC_DMA( IAC, JNH, ITK, JCS ) +   
     &       PAC * QNH * QTK * PCS * ACDC_DMA( IAC, JNH, JTK, ICS ) +   
     &       PAC * QNH * QTK * QCS * ACDC_DMA( IAC, JNH, JTK, JCS ) +   
     &       QAC * PNH * PTK * PCS * ACDC_DMA( JAC, INH, ITK, ICS ) +   
     &       QAC * PNH * PTK * QCS * ACDC_DMA( JAC, INH, ITK, JCS ) +   
     &       QAC * PNH * QTK * PCS * ACDC_DMA( JAC, INH, JTK, ICS ) +   
     &       QAC * PNH * QTK * QCS * ACDC_DMA( JAC, INH, JTK, JCS ) +   
     &       QAC * QNH * PTK * PCS * ACDC_DMA( JAC, JNH, ITK, ICS ) +   
     &       QAC * QNH * PTK * QCS * ACDC_DMA( JAC, JNH, ITK, JCS ) +   
     &       QAC * QNH * QTK * PCS * ACDC_DMA( JAC, JNH, JTK, ICS ) +   
     &       QAC * QNH * QTK * QCS * ACDC_DMA( JAC, JNH, JTK, JCS ) 


C.. Convert the unit of particle nucleation rate into [ 1 / m**3 s] 
      DNDT = DBLE(JNUC) * 1.0D06  ! (1/(m**3 s))

      RETURN    
      END SUBROUTINE ACDC_AMINE 
 

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ION_MEDIATED_NPF (  H2SO4in, ELVOCin, PRESSin, TEMPin, RHin, ASURFin, 
     &                       J5, J6, sulfmass, elvocmass, J5_2m, J6_2m )

      USE AEROMET_DATA   ! Includes CONST.EXT
      USE PRECURSOR_DATA, ONLY : PRECURSOR_MW, SULF_IDX, DMA_IDX 
      USE APM_NUCL_MOD

      IMPLICIT NONE

C   REVISION HISTORY:
C     Added by Ben Murphy to capture New Particle Formation via
C     H2SO4 and ion interactions. These rates are brought in from
C     a look-up table informed by the work of Yu et al.

C.. References:
C

C..Inputs: 
      REAL,    INTENT(IN ) :: ASURFin         ! Condesation Sink (s-1)      
      REAL,    INTENT(IN)  :: TEMPin          ! Ambient temperature [ K ]
      REAL,    INTENT(IN)  :: PRESSin         ! Ambient Pressure [Pa]
      REAL,    INTENT(IN)  :: RHin            ! Fractional Relative Humidity
      REAL,    INTENT(IN)  :: H2SO4in         ! sulfuric acid concentration [ molec/cm**3 ]
      REAL,    INTENT(IN)  :: ELVOCin         ! organic concentration [ molec/cm**3 ]

C.. Outputs:
      REAL(8), INTENT( OUT ) :: J5         ! ION-H2SO4 particle number production rate [ m-3 s-1 ]
      REAL(8), INTENT( OUT ) :: J6         ! ION-ELVOC particle number production rate [ m-3 s-1 ]
      REAL(8), INTENT( OUT ) :: sulfmass   ! mass of sulfuric acid used to make new particle [ug]
      REAL(8), INTENT( OUT ) :: elvocmass  ! mass of elvoc used to make new particle [ug]
      REAL(8), INTENT( OUT ) :: J5_2m      ! new-particle diameter squared [m**2]
      REAL(8), INTENT( OUT ) :: J6_2m      ! new-particle diameter squared [m**2]

C.. Parameters
      CHARACTER( 16 ), PARAMETER :: PNAME = 'NPF_IMN'                 

C.. Sub-Module Input Parameters
      REAL(8) :: H2SO4      ! Molecule Concentration of Sulfuric Acid [molec cm-3]
      REAL(8) :: ELVOC      ! Conc. of ELVOC [molec cm-3]
      REAL(8) :: TEMP       ! Temperature
      REAL(8) :: ASURF      ! Internal Condensation Sink
      REAL(8) :: PRESS      ! Pressure [ mb ]
      REAL(8) :: PSURF      ! Surface Pressure [ mb ]
      REAL(8) :: RH         ! Relative Humidity in Percent
      
      REAL :: MW_H2SO4        ! MW of H2SO4 in [ g / mole ]        
      REAL :: MW_ELVOC        ! MW of DMA in [ g / mole ]       
      REAL, PARAMETER :: R = 8.314             ! m3 Pa K-1 mol-1                 
      REAL, PARAMETER :: SCALEFAC = 1.0E-06    ! for [ 1 / m**3 ] ->  [ 1 / cm**3 ]
      REAL, PARAMETER :: MUG2G = 1.0E-6        ! [ ug  ] -> [ g ]                 

C.. Ion-Pair Formation Rates
      REAL(8) :: ION_RATE  ! Ion-Pair formation Rate [ ions/cm3 ]
      REAL(8) :: ION_RATE_WATER  ! Ion-Pair formation Rate [ ions/cm3 ]
      REAL(8) :: ION_RATE_LAND  ! Ion-Pair formation Rate [ ions/cm3 ]

C.. Parameters for Kirkby Ion-Org Parameterization
      REAL(8), PARAMETER :: K_KIRKBY  = 0.001366D0
      REAL(8), PARAMETER :: E1_KIRKBY = 1.566D0
      REAL(8), PARAMETER :: E2_KIRKBY = 0.1863D0
      REAL(8), PARAMETER :: J6_R_KIRKBY = 0.85D0

C.. Size of New Particles
      REAL(8) :: J5_R      ! radius of a new Ion_H2SO4 particle [nm]
      REAL(8) :: J5_R_WATER! radius of a new Ion_H2SO4 particle [nm]
      REAL(8) :: J5_R_LAND ! radius of a new Ion_H2SO4 particle [nm]
      REAL(8) :: J5_WATER  ! radius of a new Ion_H2SO4 particle [nm]
      REAL(8) :: J5_LAND   ! radius of a new Ion_H2SO4 particle [nm]
      REAL(8) :: J5_d20    ! New Particle Diamter [nm]
      REAL(8) :: J6_d20    ! New Particle Diamter [nm]

C.. Calculate Mass of Each Vapor Species Used For Particle (g/particle)
      !First compute the mass of the condensed species in the particle
      !for its size, using bulk density. Then convert to the mass of the
      !vapor species needed for that condensation. Assume Dimethylaminium
      REAL, PARAMETER :: J5_dens = 2.00E12   !Guess at the HSO4 Density (ug m-3) 
      REAL, PARAMETER :: J6_dens = 1.40E12   !Guess at ELVOA Density (ug m-3)

C.. Initialize variables                                                                                    
      J5   = 0.0                                                                                          
      J6   = 0.0                                                                                          
      
C.. The parameterization is valid at sulfuric acid conc. of 1.0E4 - 3.16E9 molec cm-3          
      H2SO4  = REAL( MIN( MAX (H2SO4in, 5.0E5), 5.0E8), 8 )
                                                                                                            
C.. The parameterization is valid at DMA conc. of 1.00e4 - 1.00E9 molec cm-3
      ELVOC  = REAL( MIN( MAX (ELVOCin, 1.00E4), 1.00E9), 8 )

C.. The parameterization is valid at temperatures of 180.00-320.0 K
      TEMP = REAL( MIN( MAX (TEMPin, 190.0), 302.0), 8 )

C.. The parameterization is valid at pressures in mb
      PRESS = REAL( MIN( MAX( PRESSin, 5.0), 1015.0), 8 ) * 0.01D0
      PSURF = REAL( SURFPRESS, 8 ) * 0.01D0

C.. The parameterization is valid at Aerosol Surface Area from 1-1000 um2/m3
      ASURF = REAL( MIN( MAX (ASURFin * 1.0E12, 1.0), 1.0E3), 8 )

C.. The parameterization is valid at RH of 0% to 100%
      RH = REAL( MIN( MAX (RHin*100.0, 0.5), 99.5), 8 )

C.. Look up Ion Pair Formation Rate from Internal Table
      CALL IONRATE( 1, PSURF, DBLE(XLON), DBLE(XLAT), PRESS, ION_RATE_LAND)
      CALL IONRATE( 0, PSURF, DBLE(XLON), DBLE(XLAT), PRESS, ION_RATE_WATER)

C.. Calculate Ion-H2SO4 Formation Rate
      CALL YUJIMN( H2SO4, RH, TEMP, ION_RATE_LAND,  ASURF, J5_LAND,  J5_R_LAND )
      CALL YUJIMN( H2SO4, RH, TEMP, ION_RATE_WATER, ASURF, J5_WATER, J5_R_WATER )

      J5   = J5_LAND   * DBLE( 1.0 - LSURF ) + J5_WATER   * DBLE( LSURF )
      J5_R = J5_R_LAND * DBLE( 1.0 - LSURF ) + J5_R_WATER * DBLE( LSURF )

C.. Calculate Ion-ELVOC Formation Rate
      ION_RATE = ION_RATE_LAND * DBLE( 1.0 - LSURF ) + ION_RATE_WATER * DBLE( LSURF ) 
      J6 = 2.D0 * ION_RATE
     &          * K_KIRKBY * ELVOC ** ( E1_KIRKBY + E2_KIRKBY / ELVOC )

C.. Convert the unit of particle nucleation rate into [ 1 / m**3 s] 
      J5 = J5 * 1.0D06  ! (1/(m**3 s))
      J6 = J6 * 1.0D06  ! (1/(m**3 s))
      J5_2m = 1.0D-18 * J5_R ** 2  !New particle diameter squared
      J6_2m = 1.0D-18 * J6_R_KIRKBY ** 2  !New particle diameter squared
      sulfmass  = 4./3.*3.1415 * ( J5_R * 1.0D-9 )** 3 * J5_dens  !Particle Mass (ug)             
      elvocmass = 4./3.*3.1415 * ( J6_R_KIRKBY * 1.0D-9 )** 3 * J6_dens  !Particle Mass (ug)             

      RETURN    
      END SUBROUTINE ION_MEDIATED_NPF 
 
 
