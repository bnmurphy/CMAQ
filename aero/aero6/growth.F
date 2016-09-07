
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
C $Header: /project/yoj/arc/CCTM/src/aero/aero5/coags.f,v 1.4 2011/10/21 16:10:13 yoj Exp $

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE HCOND3( AM0, AM1, AM2, DV, ALPHA, CBAR, F )

C  calculates the size-dependent term in the condensational-growth rate
C  expression for the 2nd and 3rd moments of a lognormal aerosol mode
C  using the harmonic mean method.  This code follows Section A2 of
C  Binkowski & Shankar (1995).
 
C  Key Subroutines/Functions called:  none
 
C  Revision History:
C     coded November 7, 2003 by Dr. Francis S. Binkowski
C     Revised November 20, 2003 by F. Binkowski to have am1 and
C     am2 as inputs
 
C  Reference:
C   1. Binkowski, F.S. and U. Shankar, The regional particulate matter
C      model 1. Model description and preliminary results, J. Geophys.
C      Res., Vol 100, No D12, 26101-26209, 1995.

      IMPLICIT NONE

C *** Includes:

      INCLUDE SUBST_CONST     ! physical and mathematical constants

C *** Arguments:

      REAL( 8 ), INTENT( IN ) :: AM0   ! zeroth moment of mode  [ #/m**3 ]
      REAL( 8 ), INTENT( IN ) :: AM1   ! first moment of mode   [ m/m**3 ]
      REAL( 8 ), INTENT( IN ) :: AM2   ! second moment of mode  [ m**2/m**3 ]
      REAL,      INTENT( IN ) :: Dv    ! molecular diffusivity of the
                                       ! condensing vapor  [ m**2/s ]
      REAL,      INTENT( IN ) :: ALPHA ! accommodation coefficient
      REAL,      INTENT( IN ) :: CBAR  ! kinetic velocity of condensing vapor [ m/s ]

      REAL( 8 ), INTENT( OUT ) :: F( 2 ) ! size-dependent term in condensational-growth
                                         ! rate: F(1) = 2nd moment [ m**2/m**3 s ]
                                         !       F(2) = 3rd moment [ m**3/m**3 s ]

C *** Local Variables:

      REAL( 8 ) :: GNC2 ! integrals used to calculate F(1) [m^2 / m^3 s]
      REAL( 8 ) :: GFM2 !

      REAL( 8 ) :: GNC3 ! integrals used to calculate F(2) [m^3 / m^3 s]
      REAL( 8 ) :: GFM3 !

      REAL( 8 ), PARAMETER :: TWOPI = 2.0D0 * PI
      REAL( 8 ), PARAMETER :: PI4 = 0.25D0 * PI

C-----------------------------------------------------------------------

C *** Implement equation A15 of Binkowski & Shankar (1995) for the
C     2nd and 3rd moments of a lognormal mode of arbitrary size.

      GNC2 = TWOPI * DV * AM0          ! 2nd moment, near-continuum
      GNC3 = TWOPI * DV * AM1          ! 3rd moment, near-continuum
      GFM2 = PI4 * ALPHA * CBAR * AM1  ! 2nd moment, free-molecular
      GFM3 = PI4 * ALPHA * CBAR * AM2  ! 3rd moment, free-molecular

C *** Implement equation A13 of Binkowski & Shankar (1995) for a
C     lognormal mode of arbitrary size.  These are the size-dependent
C     terms in the condensational-growth rate expression, given in
C     equation 7a of B&S (1995).

      F( 1 ) = GNC2 * GFM2 / ( GNC2 + GFM2 )  ! 2nd moment
      F( 2 ) = GNC3 * GFM3 / ( GNC3 + GFM3 )  ! 3rd moment

      RETURN
      END SUBROUTINE HCOND3
 

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE HCOND4( AM0, AM1, AM2, DV, DP, CBAR_V, CBAR_P, 
     &                   DIA_V, DIA_P, GF )

!  SUBROUTINE HCOND4 calculates collision rate, GF,  between a vapor and
!  and a lognormal aerosol mode. Particle size is represented by the
!  Geometric Mean Diameter of the mode.
!
!  KEY SUBROUTINES/FUNCTIONS CALLED:  none
!
!  REVISION HISTORY:
!     coded March 22, 2016 by Dr. Francis S. Binkowski
!   The form is based upon HCOND3 coded November 20, 2003
!   by Francis S. Binkowski

!  REFERENCE:
!   Nieminen, T. et al., 2010, Sub-10 nm particle growth by vapor
!   condensation - effect of vapor molecule size and particle thermal speed,
!      Atmospheric Chemistry and Physics, doi: 10.5194/acp-10-9773-2010.
      

      IMPLICIT NONE

C *** Includes:

      INCLUDE SUBST_CONST     ! physical and mathematical constants

C *** Arguments:

      REAL( 8 ), INTENT( IN ) :: AM0   ! zeroth moment of mode  [ #/m**3 ]
      REAL( 8 ), INTENT( IN ) :: AM1   ! first moment of mode   [ m/m**3 ]
      REAL( 8 ), INTENT( IN ) :: AM2   ! second moment of mode  [ m**2/m**3 ]
      REAL( 4 ), INTENT( IN ) :: Dv    ! molecular diffusivity of the
                                       ! condensing vapor  [ m**2/s ]
      REAL( 4 ), INTENT( IN ) :: Dp    ! Diffusivity of particles [m2 s-1]
      REAL( 4 ), INTENT( IN ) :: CBAR_V! kinetic velocity of condensing vapor [ m/s ]
      REAL( 4 ), INTENT( IN ) :: CBAR_P! kinetic velocity of particles [ m/s ]
      REAL( 4 ), INTENT( IN ) :: DIA_V ! diameter of vapor molecules [ m ]
      REAL( 4 ), INTENT( IN ) :: DIA_P ! diameter of particles [ m ]

      REAL( 8 ), INTENT( OUT ) :: GF(2)   ! Collision rate between vapor and particle
                                       !       [ m3 m-3 s-1 ]

C *** Local Variables:

      REAL( 4 ) :: FAC_NC, FAC_FM
      REAL( 8 ) :: GNC2, GNC3 ! integrals used to calculate F(1) [m^2 / m^3 s]
      REAL( 8 ) :: GFM2, GFM3 !
      REAL( 8 ) :: DIA_V8 ! diameter of vapor molecules [ m ]

      REAL( 8 ), PARAMETER :: TWOPI = 2.0D0 * PI
      REAL( 8 ), PARAMETER :: PI4 = 0.25D0 * PI

C-----------------------------------------------------------------------

      DIA_V8 = REAL( DIA_V,8 )

C *** Implement equations 1 and 2 of Nieminen et al., 2010

      FAC_NC = TWOPI * ( Dv + Dp )
      FAC_FM = PI4 * SQRT( CBAR_V **2 + CBAR_P**2 )

      !GNC is for the near-continuum regime
      GNC2 = REAL( FAC_NC,8 ) * ( AM0 )
      GNC3 = REAL( FAC_NC,8 ) * ( AM1 + AM0*DIA_V8 )

      !GFM is for the free-molecular regime
      GFM2 = REAL( FAC_FM,8 ) * ( AM1 + 2.0 * AM0 * DIA_V8 )
      GFM3 = REAL( FAC_FM,8 ) * ( AM2 + 2.0 * AM1 * DIA_V8 + AM0 * DIA_V8 **2 )

C *** Use Harmonic Mean method to account for transition from
C     Free-molecular to continuum regimes
      GF(1) = gnc2 * gfm2 / ( gnc2 + gfm2 )
      GF(2) = gnc3 * gfm3 / ( gnc3 + gfm3 )

      RETURN
      END SUBROUTINE HCOND4

#ifdef AERO_NUC_GROWTH
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE GHCOND( AM0, DIA_P, DV, DP, CBAR_V, CBAR_P, 
     &                   DIA_V, GF )

!  SUBROUTINE GHCOND calculates the condensation sink of a particle mode
!     using the Gauss-Hermite inumerical quadrature approximation to the 
!     integral of the log-normal distribution. This code replaces the
!     growth method is Binkowski and Shanker (1995) FSB coded July 19, 
!     2016 by Dr. Francis S. Binkowski CEMPD, The Institute for the 
!     Environment, University of North Carolina at Chapel Hill.
!
!     The code is based upon HCOND7, which in turn is based upon
!     Lehtinen & Kulmala(2003)
!  REFERENCES:
!
!     Binkowski, F.S. and U. Shankar,1995,  The regional particulate matter
!      model 1. Model description and preliminary results, J. Geophys.
!      Res., Vol 100, No D12, 26101-26209, 1995
!
!     Lehtinen, K.E.J, and M. Kulmala, 2003, A model for particle
!      formation and growth in the atmosphere with molecular resolution 
!      in size, Atmospheric Chemistry and Physics, doi: 10.5194/acp-3-257-2003
!
!     Nieminen, T. et al., 2010, Sub-10 nm particle growth by vapor
!      condensation - effect of vapor molecule size and particle thermal 
!      speed, Atmospheric Chemistry and Physics, doi: 10.5194/acp-10-9773-2010.
!
!
!  KEY SUBROUTINES/FUNCTIONS CALLED:  none
!
!  REVISION HISTORY:
!     coded March 22, 2016 by Dr. Francis S. Binkowski
!   The form is based upon HCOND3 coded November 20, 2003
!   by Francis S. Binkowski

      IMPLICIT NONE

C *** Includes:

      INCLUDE SUBST_CONST     ! physical and mathematical constants

C *** Arguments:

      REAL( 8 ), INTENT( IN ) :: AM0   ! zeroth moment of mode  [ #/m**3 ]
      REAL( 8 ), INTENT( IN ) :: AM1   ! first moment of mode   [ m/m**3 ]
      REAL( 8 ), INTENT( IN ) :: AM2   ! second moment of mode  [ m**2/m**3 ]
      REAL( 4 ), INTENT( IN ) :: Dv    ! molecular diffusivity of the
                                       ! condensing vapor  [ m**2/s ]
      REAL( 4 ), INTENT( IN ) :: Dp    ! Diffusivity of particles [m2 s-1]
      REAL( 4 ), INTENT( IN ) :: CBAR_V! kinetic velocity of condensing vapor [ m/s ]
      REAL( 4 ), INTENT( IN ) :: CBAR_P! kinetic velocity of particles [ m/s ]
      REAL( 4 ), INTENT( IN ) :: DIA_V ! diameter of vapor molecules [ m ]
      REAL( 4 ), INTENT( IN ) :: DIA_P ! diameter of particles [ m ]

      REAL( 8 ), INTENT( OUT ) :: GF(2)   ! Collision rate between vapor and particle
                                       !       [ m3 m-3 s-1 ]

C *** Local Variables:

      REAL( 4 ) :: FAC_NC, FAC_FM
      REAL( 8 ) :: GNC2, GNC3 ! integrals used to calculate F(1) [m^2 / m^3 s]
      REAL( 8 ) :: GFM2, GFM3 !
      REAL( 8 ) :: DIA_V8 ! diameter of vapor molecules [ m ]

      REAL( 8 ), PARAMETER :: TWOPI = 2.0D0 * PI
      REAL( 8 ), PARAMETER :: PI4 = 0.25D0 * PI

C-----------------------------------------------------------------------
      real*8, intent(in)  :: M0     ! Modal number concentration           [ #/m**3 ]
      real*8, intent(in)  :: dia_p  ! diameter of particle                 [ m ]
      real*8, intent(in)  :: xlnsig  ! log(GSD)   
      real*8, intent(in)  :: Dv     ! diffusivity of the  vapor            [m**2 / s]
      real*8, intent(in)  :: Dp     ! diffusivity of the particles         [m**2 / s]          
      real*8, intent(in)  :: cbar_v ! kinetic velocity of condensing vapor [m / s]   
      real*8, intent(in)  :: cbar_p ! kinetic velocity of the particles    [m / s]
      real*8, intent(in)  :: dia_v  ! diameter of vapor molecule           [ m ]
     
! *** output: 
                
      real*8, intent(out) :: GF ! size-dependent term in 
                                           ! condensational-growth rate
                                           ! GF = 3rd moment [m^3 / m^3 s]      
 
! *** internal
      real*8, parameter :: one   = 1.0D0, two = 2.0D0, sqrt_2 = sqrt( two )
      real*8, parameter :: three = 3.0D0, four = 4.0D0
      real*8, parameter :: pi    = four * atan(one)
      real*8, Parameter :: sqrt_pi = sqrt( pi )
      real*8, parameter :: twopi = two * pi
      real*8, parameter :: four_thirds = four / three
      real*8, parameter :: A_FS = 1.77D0
      real*8            :: Kn0, Kn_m, Kn_p ! Knudsen number
      real*8            :: betafs0, betafs_m, betafs_p ! F-S factor
      real*8            :: gnc0, gnc_m, gnc_p
      real*8            :: fac1, fac2, fac3
      real*8            :: x0, x_m, x_p,xi, wxi,  xf, aa1
      real*8            :: sum
      
       
! FSB Pre-calculate Gauss-Hermite abscissas and weights for three point
!     numerical quadrature.
!     Coded June 25, 2016 by Dr. Francis S. Binkowski
! Reference: 
!      Hildebrand, F.B., 1956, Introduction to Numerical Analysis, 
!      McGraw-Hill Book Company =, New York, 511 pp. Citation is to 
!      Pages 327-329. Method is the same in Second edition of this pook
!      published in 1974 where the page numbers are 395-397.

       real*8, parameter :: one_half = one / two, two_thirds = two / three
       real*8, parameter :: six = 6.0d0, sqrt_six = sqrt( six ), one_sixth = one / six
       real*8, parameter :: ghx0_3 = 0.0d0  
       real*8, parameter :: ghw0_3 = two_thirds * sqrt_pi
       real*8, parameter :: ghxi_3 = one_half * sqrt_six
       real*8, parameter :: ghwi_3 = one_sixth * sqrt_pi
       
! FSB Start Code 

      fac1 = (Dp +Dv)
      fac2 = sqrt( cbar_v**2 + cbar_p**2 )        
       
       
 ! FSB now start the integration code
!    We are using an odd number of points points here

       aa1 = sqrt_2 * xlnsig   ! This 1.0 / Sqrt( A ) in derivation of the integral
                              ! where A = 1.0 / ( 2.0 * xlnsg**2 ) 
! Then  = dia_p * exp[ u / sqrt(A) ]
! For Gauss-Hermite Quadrature u = xi 
! Therefore, xf = exp( xi / sqrt(A) ),
!  or xf = exp( xi * aa1 )
 
! Do center point xi = 0.0 first

       xi       = ghx0_3
       wxi      = ghw0_3
!      xf       = exp( xi * aa1 )     
! FSB  xf = 1.0 because xi = 0.0
       x0      =  dia_p 
       
       Kn0     = six * fac1 / ( (fac2) * ( x0 + dia_v) )
       betafs0 = ( one + kn0 ) /  ( one + A_FS * kn0 + four_thirds * kn0 * kn0)
       gnc0    = twopi * fac1 * (x0 + dia_v) * betafs0
       
       sum = wxi * gnc0
       
! FSB now do remaining points 

       xi      = ghxi_3
       wxi     = ghwi_3
       xf      = exp( xi * aa1 )
       x_p     = dia_p *xf
       x_m     = dia_p / xf
       kn_p  = six * fac1 / ( (fac2) * ( x_p + dia_v) )
       kn_m  = six * fac1 / ( (fac2) * ( x_m + dia_v) )
       betafs_m = ( one + kn_m ) / ( one + A_FS * kn_m + four_thirds * kn_m * kn_m)
       betafs_p = ( one + kn_p ) / ( one + A_FS * kn_p + four_thirds * kn_p * kn_p)
       gnc_m    = twopi * fac1 * (x_m + dia_v) * betafs_m
       gnc_p    = twopi * fac1 * (x_p + dia_v) * betafs_p
       sum = wxi * ( gnc_p + gnc_m)
       GF = sum * M0     
              
       return
              
      end subroutine get_GF_GH
 


      DIA_V8 = REAL( DIA_V,8 )


      RETURN
      END SUBROUTINE GHCOND
#endif 




