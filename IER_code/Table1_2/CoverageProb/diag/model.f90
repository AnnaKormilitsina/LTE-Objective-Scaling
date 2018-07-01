MODULE MODEL  
USE MPARAMS
USE FUNCTIONS  
IMPLICIT NONE  
INTEGER, PUBLIC, PARAMETER :: noutput = 2,  nr = 1, npai = 3, sth = 1, seps = 2, sgam = 3 
INTEGER, PUBLIC, PARAMETER :: sizef = 6, sizex = 3,  sizey = 3 
INTEGER, PUBLIC, PARAMETER :: nSUBPER = 1, nz(nSUBPER) = (/ 0 /) , ntheta(nSUBPER) = (/ 0 /) , sizetheta = 0, sizez = 0
 
CONTAINS  
 
SUBROUTINE fxfyfxpfypf(Thetta, nf, nfx, nfxp, nfy, nfyp, signout)  
IMPLICIT NONE   
REAL(8), INTENT(IN)  :: Thetta(num_thetta)
REAL(8), INTENT(OUT) :: nf(sizef), nfx(sizef*(sizex-sizetheta)), nfxp(sizef*(sizex-sizetheta)), nfy(sizef*sizey), nfyp(sizef*sizey)  
REAL(8)              :: yout(10), BETTA, R, PAI, Y, RHO_th, RHO_eps, RHO_gam, ALFA_PAI, KAPA, ALFA_Y
INTEGER, INTENT(OUT) :: signout 

INTEGER :: i

CALL calibration(Thetta,  signout, yout) 

IF (signout == 1) THEN
   nf   = 0.0D+00
   nfx  = 0.0D+00
   nfy  = 0.0D+00
   nfxp = 0.0D+00
   nfyp = 0.0D+00
   return
END IF

BETTA    =  yout(1)
R        =  yout(2)
PAI      =  yout(3)
Y        =  yout(4)
RHO_th   =  yout(5)
RHO_eps  =  yout(6)
RHO_gam  =  yout(7)
ALFA_PAI =  yout(8)
KAPA     =  yout(9)
ALFA_Y   =  yout(10)

IF ( 1 == 2) THEN 
 
  write(*,*) 'BETA'     , BETTA
  write(*,*) 'R'        , R
  write(*,*) 'PAI'      , PAI
  write(*,*) 'Y'        , Y
  write(*,*) 'RHO_th'   , RHO_th
  write(*,*) 'RHO_eps'  , RHO_eps
  write(*,*) 'RHO_gam'  , RHO_gam
  write(*,*) 'ALFA_PAI' , ALFA_PAI
  write(*,*) 'KAPA'     , KAPA
  write(*,*) 'ALFA_Y'   , ALFA_Y

END IF

nfx    = 0.0D+0 
nfxp   = 0.0D+0 
nfy    = 0.0D+0 
nfyp   = 0.0D+0 


nf(1) = (BETTA*R)/(PAI*Y) - 1.0D+0/Y
DO i = 1, sizef 
 IF ( ABS(nf(i)) > 1.0D-5) THEN 
signout = 1 
 nfx = 0.0D00 
 nfy = 0.0D+00 
 nfxp = 0.0D+00 
 nfyp = 0.0D+00 
 write(*,*) 'nf = ', nf 
 RETURN 
 END IF 
 END DO 
nfx(3) = 1.0D+0 
nfx(4) = -RHO_th 
nfx(7) = 1.0D+0 
nfx(11) = -RHO_eps 
nfx(14) = 1.0D+0 
nfx(18) = -RHO_gam 
nfy(1) = (BETTA*R)/(PAI*Y) 
nfy(3) = -1.0D+0 
nfy(7) = 1.0D+0/Y 
nfy(8) = KAPA 
nfy(9) = ALFA_Y 
nfy(14) = -1.0D+0 
nfy(15) = ALFA_PAI 
nfy(15) = ALFA_PAI 
nfyp(7) = -(BETTA*R)/(PAI*Y) 
nfyp(13) = -(BETTA*R)/(PAI*Y) 
nfyp(14) = BETTA 
nfxp(18) = 1.0D+0
nfxp(11) = 1.0D+0
nfxp(4) = 1.0D+0

END SUBROUTINE  fxfyfxpfypf

SUBROUTINE calibration( Thetta, signout, params)
IMPLICIT NONE
REAL(8), INTENT(IN)  :: Thetta(num_thetta)
INTEGER, INTENT(OUT) :: signout
REAL(8), INTENT(OUT) :: params(10)
REAL(8)              :: BETTA, R, PAI, Y, RHO_th, RHO_eps, RHO_gam, ALFA_PAI, KAPA, ALFA_Y

REAL(8) :: ALFA, TH, GAM, EPS


signout = 0

ALFA_PAI = Thetta(1)
ALFA_Y   = 0.0D+0!Thetta(2)
BETTA    = Thetta(2)
ALFA     = Thetta(3)
KAPA     = (1.0D+0-ALFA)*(1.0D+0-ALFA*BETTA)/ALFA
RHO_th   = Thetta(4)
RHO_eps  = Thetta(5)
RHO_gam  = Thetta(6)
Y        = 1.0D+0
PAI      = 1.0D+0
R        = PAI/BETTA
TH       = 1.0D+0
GAM      = 1.0D+0
EPS      = 1.0D+0

params(1)  = BETTA    
params(2)  = R        
params(3)  = PAI      
params(4)  = Y        
params(5)  = RHO_th   
params(6)  = RHO_eps  
params(7)  = RHO_gam  
params(8)  = ALFA_PAI 
params(9)  = KAPA     
params(10) = ALFA_Y   

END SUBROUTINE calibration
END MODULE MODEL
