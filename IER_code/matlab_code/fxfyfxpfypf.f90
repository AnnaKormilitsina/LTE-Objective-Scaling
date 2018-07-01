SUBROUTINE fxfyfxpfypf(Thetta, nf, nfx, nfxp, nfy, nfyp, signout)  
sizef = 7 
sizex = 4 
sizey = 3 
sth  = 2 
seps = 3 
sgam = 4 
noutput = 2 
npai    = 3 
nr      = 1 
nf  = 0.0D+0  
 nfx  = 0.0D+00 
 nfxp = 0.0D+00 
nfy = 0.0D+00 
 nfyp = 0.0D+00 
nf(1) = (BETTA*R)/(PAI*Y) - 1/Y
DO i = 1, sizef 
 IF ( ABS(nf(i)) > 1.0D-5) THEN 
signout = 1 
 nfx = 0.0D00 
 nfy = 0.0D+00 
 nfxp = 0.0D+00 
 nfyp = 0.0D+00 
 write(*,*)  nf 
 RETURN 
 END IF 
 END DO 
nfx(3) = ALFA_R 
  
nfx(10) = 1 
  
nfx(11) = -RHO_th 
  
nfx(15) = 1 
  
nfx(19) = -RHO_eps 
  
nfx(23) = 1 
  
nfx(27) = -RHO_gam 
  
nfy(7) = R 
  
nfy(8) = 1/Y 
  
nfy(9) = KAPA 
  
nfy(10) = ALFA_Y 
  
nfy(16) = -1 
  
nfy(17) = ALFA_PAI 
  
nfxp(1) = (BETTA*R)/(PAI*Y) 
  
nfxp(3) = -1 
  
nfxp(7) = -R 
  
nfxp(11) = 1 
  
nfxp(19) = 1 
  
nfxp(27) = 1 
  
nfyp(8) = -(BETTA*R)/(PAI*Y) 
  
nfyp(15) = -(BETTA*R)/(PAI*Y) 
  
nfyp(16) = BETTA 
  
END SUBROUTINE fxfyfxpfypf