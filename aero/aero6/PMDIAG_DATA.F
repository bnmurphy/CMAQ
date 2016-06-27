
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

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE PMDIAG_DATA_MODULE

c-----------------------------------------------------------------------
c Revision History:
c David Wong 11 May 16: created
c-----------------------------------------------------------------------

        IMPLICIT NONE

        TYPE PMDIAG_SPC_RECORD_TYPE
          CHARACTER (LEN = 16) :: SPC_NAME
          INTEGER              :: VAR_TYPE
          CHARACTER (LEN = 16) :: SPC_UNIT
          CHARACTER (LEN = 80) :: SPC_DESC
        END TYPE PMDIAG_SPC_RECORD_TYPE

        TYPE (PMDIAG_SPC_RECORD_TYPE), PARAMETER :: PMDIAG_SPC_RECORD(72) = (/
     &    PMDIAG_SPC_RECORD_TYPE( 'STDEVNUC        ', 5, 'na              ', 'NPF mode standard deviation (dry and wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'STDEVATKN       ', 5, 'na              ', 'Aitken mode standard deviation (dry and wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'STDEVACC        ', 5, 'na              ', 
     &                            'Accumulation mode standard deviation (dry and wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'STDEVCOR        ', 5, 'na              ', 'Coarse mode standard deviation (dry and wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DGNUC_DRY       ', 5, 'um              ', 'NPF mode mean diameter (dry)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DGATKN_DRY      ', 5, 'um              ', 'Aitken mode mean diameter (dry)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DGACC_DRY       ', 5, 'um              ', 'Accumulation mode mean diameter (dry)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DGCOR_DRY       ', 5, 'um              ', 'Coarse mode mean diameter (dry)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DGNUC_WET       ', 5, 'um              ', 'NPF mode mean diameter (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DGATKN_WET      ', 5, 'um              ', 'Aitken mode mean diameter (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DGACC_WET       ', 5, 'um              ', 'Accumulation mode mean diameter (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DGCOR_WET       ', 5, 'um              ', 'Coarse mode mean diameter (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M2NUC_WET       ', 5, 'm**2/m**3       ', 'NPF mode 2nd moment (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M2ATKN_WET      ', 5, 'm**2/m**3       ', 'Aitken mode 2nd moment (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M2ACC_WET       ', 5, 'm**2/m**3       ', 'Accumulation mode 2nd moment (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M2COR_WET       ', 5, 'm**2/m**3       ', 'Coarse mode 2nd moment (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M3NUC_DRY       ', 5, 'm**3/m**3       ', 'NPF mode 3rd moment (dry)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M3ATKN_DRY      ', 5, 'm**3/m**3       ', 'Aitken mode 3rd moment (dry)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M3ACC_DRY       ', 5, 'm**3/m**3       ', 'Accumulation mode 3rd moment (dry)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M3COR_DRY       ', 5, 'm**3/m**3       ', 'Coarse mode 3rd moment (dry)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M3NUC_WET       ', 5, 'm**3/m**3       ', 'NPF mode 3rd moment (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M3ATKN_WET      ', 5, 'm**3/m**3       ', 'Aitken mode 3rd moment (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M3ACC_WET       ', 5, 'm**3/m**3       ', 'Accumulation mode 3rd moment (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'M3COR_WET       ', 5, 'm**3/m**3       ', 'Coarse mode 3rd moment (wet)' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'RH              ', 5, 'na              ', 'Relative humidity' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM1NU           ', 5, 'fraction        ', 'PM1 fraction of NPF mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM1AT           ', 5, 'fraction        ', 'PM1 fraction of Aitken mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM1AC           ', 5, 'fraction        ', 'PM1 fraction of accumulation mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM1CO           ', 5, 'fraction        ', 'PM1 fraction of coarse mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM25NU          ', 5, 'fraction        ', 'PM2.5 fraction of NPF mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM25AT          ', 5, 'fraction        ', 'PM2.5 fraction of Aitken mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM25AC          ', 5, 'fraction        ', 'PM2.5 fraction of accumulation mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM25CO          ', 5, 'fraction        ', 'PM2.5 fraction of coarse mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM10NU          ', 5, 'fraction        ', 'PM10 fraction of NPF mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM10AT          ', 5, 'fraction        ', 'PM10 fraction of Aitken mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM10AC          ', 5, 'fraction        ', 'PM10 fraction of accumulation mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'PM10CO          ', 5, 'fraction        ', 'PM10 fraction of coarse mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'AMSNU           ', 5, 'fraction        ', 'AMS Transmission Fraction of NPF mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'AMSAT           ', 5, 'fraction        ', 'AMS Transmission Fraction of Aitken mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'AMSAC           ', 5, 'fraction        ', 'AMS Transmission Fraction of Accumulation mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'AMSCO           ', 5, 'fraction        ', 'AMS Transmission Fraction of Coarse mode' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'GAMMA_N2O5      ', 5, 'na              ', 
     &                            'Fine Mode N2O5 Heterogeneous reaction probability' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'GAMMA_N2O5K     ', 5, 'na              ', 
     &                            'Coarse Mode N2O5 Heterogeneous reaction probability' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'YIELD_CLNO2     ', 5, 'na              ', 'Fine Mode CLNO2 Heterogeneous reaction yield' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'YIELD_CLNO2K    ', 5, 'na              ', 'Coarse Mode CLNO2 Heterogeneous reaction yield' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'GAMMA_IEPOX     ', 5, 'na              ', 'IEPOX heterogeneous uptake coefficient' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DENS_NUC        ', 5, 'kg m-3          ', 'Bulk Density of NPF-Mode Particles' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DENS_ATKN       ', 5, 'kg m-3          ', 'Bulk Density of Aitken-Mode Particles' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DENS_ACC        ', 5, 'kg m-3          ', 'Bulk Density of Accumulation-Mode Particles' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'DENS_COR        ', 5, 'kg m-3          ', 'Bulk Density of Coarse-Mode Particles' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'NPF_BINARY      ', 5, 'N cm-3 s-1      ', 'NPF Rate from Binary Production' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'NPF_TERNARY     ', 5, 'N cm-3 s-1      ', 'NPF Rate from Ternary Production' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'NPF_SO4_DMA     ', 5, 'N cm-3 s-1      ', 'NPF Rate from DMA-H2SO4 Production' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'NPF_SO4_ORG     ', 5, 'N cm-3 s-1      ', 'NPF Rate from Org-H2SO4 Production' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'NPF_ORG_ORG     ', 5, 'N cm-3 s-1      ', 'NPF Rate from Org-Org Production' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'NPF_RICCO       ', 5, 'N cm-3 s-1      ', 'NPF Rate from H2SO4-Org (Riccobono) Production' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'M2_BINARY       ', 5, 'm2 m-3 s-1      ', '2nd Moment Change from Binary Production' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'M2_TERNARY      ', 5, 'm2 m-3 s-1      ', '2nd Moment Change from Ternary Production' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'M2_SO4_DMA      ', 5, 'm2 m-3 s-1      ', '2nd Moment Change from DMA-H2SO4 Production' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'M2_SO4_ORG      ', 5, 'm2 m-3 s-1      ', '2nd Moment Change from Org-H2SO4 Production' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'M2_ORG_ORG      ', 5, 'm2 m-3 s-1      ', '2nd Moment Change from Org-Org Production' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'M2_RICCO        ', 5, 'm2 m-3 s-1      ', '2nd Moment Change from H2SO4-Org (Riccobono) Prod' ), 
     &    PMDIAG_SPC_RECORD_TYPE( 'dH2SO4_Binary   ', 5, 'ug m-3          ', 'Change in H2SO4(g) from Binary NPF' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'dH2SO4_Ternary  ', 5, 'ug m-3          ', 'Change in H2SO4(g) from Ternary NPF' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'dH2SO4_SO4_DMA  ', 5, 'ug m-3          ', 'Change in H2SO4(g) from DMA NPF' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'dH2SO4_SO4_ORG  ', 5, 'ug m-3          ', 'Change in H2SO4(g) from H2SO4-ORG NPF' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'dNH3_Ternary    ', 5, 'ug m-3          ', 'Change in NH3(g) from Ternary NPF' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'dDMA_SO4_DMA    ', 5, 'ug m-3          ', 'Change in DMA(g) from SO4-DMA NPF' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'dORG_ORG_ORG    ', 5, 'ug m-3          ', 'Change in Org(g) from Org-Org NPF NPF' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'dORG_RICCO      ', 5, 'ug m-3          ', 'Change in Org(g) from Riccobono Params' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'k_prtcl_iepox   ', 5, '1/s             ', 'IEPOX pseudo particle phase reaction rate const' ),
     &    PMDIAG_SPC_RECORD_TYPE( 'GAMMA_IMAE      ', 5, 'na              ', 'IMAE heterogeneous uptake coefficient' ) /)

        INTEGER :: A_PMDIAG_SPC_MAP (72), NUM_PMDIAG_SPC

      END MODULE PMDIAG_DATA_MODULE
