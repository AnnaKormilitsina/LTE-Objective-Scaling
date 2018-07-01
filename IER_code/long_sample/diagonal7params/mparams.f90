MODULE mparams

IMPLICIT NONE

!INTEGER, PUBLIC, PARAMETER :: num_thetta = 10, num_MOM = 3,   num_shocks = 3, nlags = 4, size_data = 5000! 200!225 
INTEGER, PUBLIC, PARAMETER ::  num_thetta = 7,  num_MOM = 3,   num_shocks = 3, nlags = 4, size_data =  5000
INTEGER, PARAMETER :: size_beta = 5, size_gama = 1 , size_normal = 2, size_uniform = 0, size_igama = 3
!                                  x    ALFA_R      ALFA_PI ALFA_Y     KAPPA       RHO_TH   RHO_EPS  RHO_GAM  STD_TH     STD_EPS    STD_GAM %
!REAL(8), PARAMETER :: Lower_Bound(num_thetta) = (/ 0.0D+0, 0.0D+0, 0.0D+0, 0.0D+0, 0.0D+00,  0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00  /)
!REAL(8), PARAMETER :: Upper_Bound(num_thetta) = (/ 1.0D+0, 10000.0D+0, 100.0D+0, 100000.0D+0, 1.0D+0,  1.0D+00, 1.0D+00, 1000000.0D+0, 1000000.0D+0,  1000000.0D+0/)

REAL(8), PARAMETER :: Lower_Bound(num_thetta) = (/  0.0D+0, 0.0D+00,  0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00  /)
REAL(8), PARAMETER :: Upper_Bound(num_thetta) = (/  100000.0D+0, 1.0D+0,  1.0D+00, 1.0D+00, 1000000.0D+0, 1000000.0D+0,  1000000.0D+0/)

!                                            KAPPA   RHO_TH  RHO_EPS  RHO_GAM  STD_TH     STD_EPS    STD_GAM %
!REAL(8), PARAMETER :: Pmean(num_thetta) = (/ 0.7D+0, 0.5D+0, 0.15D+0, 0.7D+0, 0.8D+0, 0.8D+00, 0.8D+0,  1.0D+0, 1.00D+0,  1.0D+0 /)
!REAL(8), PARAMETER :: Pstd(num_thetta)  = (/ 0.2D+0, 0.2D+0, 0.05D+0, 0.2D+0,  0.1D+0, 0.1D+0, 0.1D+00,  0.5D+0,   0.5D+0,  0.5D+0 /)

REAL(8), PARAMETER :: Pmean(num_thetta) = (/  0.7D+0, 0.8D+0, 0.8D+00, 0.8D+0,  1.0D+0, 1.00D+0,  1.0D+0 /)
REAL(8), PARAMETER :: Pstd(num_thetta)  = (/  0.2D+0, 0.1D+0, 0.1D+0,  0.1D+00, 0.5D+0, 0.5D+0,   0.5D+0 /)

INTEGER, PUBLIC, PARAMETER :: v_beta(size_beta)       = 0!(/ 1, 5, 6 /)
INTEGER, PUBLIC, PARAMETER :: v_gama(size_gama)       = 0!(/2,3, 4, /)
INTEGER, PUBLIC, PARAMETER :: v_igama(size_igama)     = 0!(/  7, 8, 9 /)
INTEGER, PUBLIC, PARAMETER :: v_normal(size_normal)   = 0
INTEGER, PUBLIC, PARAMETER :: v_uniform(size_uniform) = 0

END MODULE mparams
