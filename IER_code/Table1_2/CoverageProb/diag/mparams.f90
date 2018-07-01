MODULE mparams

IMPLICIT NONE

INTEGER, PUBLIC, PARAMETER ::  num_thetta = 2,  num_MOM = 2, size_data = 200
INTEGER, PARAMETER :: size_beta = 0, size_gama = 1, size_normal = 1, size_uniform = 0, size_igama = 0

!                                                 Mean       StD
REAL(8), PARAMETER :: Lower_Bound(num_thetta) = (/-1000.0D+00,   0.0D+00  /)
REAL(8), PARAMETER :: Upper_Bound(num_thetta) = (/1000.0D+00, 1000.0D+00/)

REAL(8), PARAMETER :: Pmean(num_thetta) = (/ 0.5D+0,     0.25D+0/)
REAL(8), PARAMETER :: Pstd(num_thetta)  = (/ 0.2D+0,     0.1D+0/)

INTEGER, PUBLIC, PARAMETER :: v_beta(size_beta)       = 0
INTEGER, PUBLIC, PARAMETER :: v_gama(size_gama)       = 0!2
INTEGER, PUBLIC, PARAMETER :: v_igama(size_igama)     = 0
INTEGER, PUBLIC, PARAMETER :: v_normal(size_normal)   = 0!1
INTEGER, PUBLIC, PARAMETER :: v_uniform(size_uniform) = 0

END MODULE mparams
