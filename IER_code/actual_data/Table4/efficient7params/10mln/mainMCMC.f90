PROGRAM mainMCMC

USE EMPIRICAL
USE MPARAMS
USE FUNCTIONS
USE MODEL_SOL
USE MODEL

IMPLICIT NONE

INTEGER, PARAMETER :: Ndraws = 10000000, Whentosave = Ndraws, SaveEach = 100
REAL(8) ::  sig_xi 
REAL(8), DIMENSION(num_thetta) :: xi, x_save, randn, Par1, Par2

REAL(8) :: x(Ndraws+1,num_thetta),  Ln_x(Ndraws+1), xstart(num_thetta)
INTEGER :: i, j, k, index_acc, equilibrium, signout,  acceptance, acceptance_all(Ndraws/10000 )
REAL(8) :: Ln_xi(1,1),  f_xi, f_x, rho(1,1)
REAL    :: rand, rand_temp(Ndraws)

REAL(8) :: MOMdata(num_MOM**2*(nlags+1) - num_MOM,1), invSIG_data(num_MOM**2*(nlags+1)-num_MOM,num_MOM**2*(nlags+1)-num_MOM), Vmom( num_MOM**2*(nlags+1) - num_MOM,num_MOM**2*(nlags+1) - num_MOM )
REAL(8) :: D(num_thetta,num_thetta), mdot(num_MOM**2*(nlags+1) - num_MOM,num_thetta), Vg(num_thetta,num_thetta) , mean_est(1,num_thetta)

REAL(8) :: gx(sizey,sizex), hx(sizex,sizex), grdh(num_thetta), Vtheta(num_thetta,num_thetta),  x_mean(Ndraws+1,num_thetta)

INTEGER :: INFO_HESS, INFO
REAL(8) :: multiplier, grad(num_MOM**2*(nlags+1)-num_MOM,num_thetta)

REAL(8) :: shock_place(sizex), shock(sizex), x0(sizex), data(sizef)
REAL(8) :: data_sample(size_data,num_MOM)

character(30) ::filename
CHARACTER(8)  :: arg, arg1, arg2, arg3
INTEGER       :: enter_i, enter_j, nsig_xi, enter_k

CALL get_command_argument(1,arg1)
read(arg1,FMT='(I5)') enter_i

CALL get_command_argument(2,arg1)
read(arg1,FMT='(I5)') enter_j
multiplier = dble(enter_j)
arg1 = adjustl(arg1)

CALL get_command_argument(3,arg3)
read(arg3,FMT='(I5)') nsig_xi

sig_xi = 0.0D+0
DO i = 1, nsig_xi
 sig_xi = sig_xi + 0.0001D+0
END DO

write(*,*) 'sig_xi',sig_xi
write(*,*) 'multiplier', multiplier

!INITIALIZATION
!=========================================================================================
IF (enter_i == 1) THEN
  x(1,:) = Pmean
xstart = Pmean
ELSE
  OPEN(unit = 1 ,file = 'xmean_est.txt', status = "old", action = "read")
  READ( 1, *, end = 7 ) xstart
  7 continue
  CLOSE(1)
  x(1,:) = xstart
END IF

write(*,*) 'x1 = ', x(1,:)
CALL gx_hx(gx,hx, x(1,:), signout, equilibrium)

write( arg ,'(I8)') size_data
arg = adjustl(arg)

IF ( 1 == 2) THEN

    CALL RANDOM_SEED()
     shock_place       = 0.0D+00

     x0 = 0.0D+00
     DO i = 1, size_data

     CALL gx_hx(gx,hx, x(1,:), signout, equilibrium)

     shock_place(sth)  = (x(1,num_thetta-2))*0.01D+0
     shock_place(seps) = (x(1,num_thetta-1))*0.01D+0
     shock_place(sgam) = (x(1,num_thetta))*0.01D+0

       DO j = 1, sizex
          shock(j)  = shock_place(j)*RANDOM_NORMAL()
       END DO
       x0                       = MATMUL(hx,x0) + shock! shock_place*RANDOM_NORMAL()
       data(1 : sizey )         = MATMUL(gx,x0)
       data( sizey + 1 : sizey+sizex) = x0
       data_sample(i,1:num_MOM) = dble(data( (/ nr, noutput, npai /)))
     END DO
     
     write(filename, '(A)') 'mdata_'//trim(arg)//'.txt'
     OPEN(unit = 22,file = filename,action = "write")
     DO i = 1, num_MOM
        DO j = 1, size_data
          WRITE(22,*)  data_sample(j,i)
        END DO
     END DO
     CLOSE(22)

ELSE

!  write(filename, '(A)') 'mdata_hpf_'//trim(arg)//'.txt'
   write(filename, '(A)') 'mdata_'//trim(arg)//'.txt'
   OPEN(unit = 2 ,file = filename, status = "old", action = "read")
   READ( 2, *, end = 19 ) data_sample
   19 continue
   CLOSE(2)

END IF

!=========================================================================================
CALL RANDOM_SEED() 
x(1,:) = xstart

CALL sample_mom(MOMdata,invSIG_data, Vmom, num_MOM**2*(nlags+1)-num_MOM , data_sample)
write(*,*) 'x 1', x(1,:)

!invSIG_data = Vmom
!DO i = 1, num_MOM**2*(nlags+1)-num_MOM
!     DO j = 1, num_MOM**2*(nlags+1)-num_MOM
!        IF (i == j) THEN
!        ELSE 
!          invSIG_data(i,j) = 0.0D+0
!      END IF
!     END DO
!END DO
!CALL MATINV(invSIG_data,  num_MOM**2*(nlags+1)-num_MOM, info)

CALL OBJECTIVE(Ln_xi, x(1,:), MOMdata, invSIG_data,  num_MOM**2*(nlags+1)-num_MOM, equilibrium, signout)
write(*,*) 'Ln_xi', Ln_xi

IF ( 1 == 1 ) THEN

  write(filename, '(A)') 'xinfoMOM_'//trim(arg1)//'.txt'
  OPEN(unit   = 33 ,file = filename, action = "write")

  Ln_x(1)  = Ln_xi(1,1)

  acceptance = 0
  acceptance_all = 0
  index_acc  = 1

  IF (signout == 1 .OR. equilibrium /= 1 ) THEN
      WRITE(33,*) '; signout ==' , signout, '; equilibrium = ', equilibrium
     x              = 0.0D+0
     Ln_x           = 0.0D+0
  ELSE
     grdh = 1.0D-04*x(1,:)
     CALL HESSIAN_SYM(D, x(1,:), grdh, MOMdata, invSIG_data, num_MOM**2*(nlags+1)-num_MOM, INFO_HESS)

     IF ( SUM(D) == 0.0D+00) THEN
        WRITE(*,*) 'Warning: D failed in  HESSIAN_SYM'
     END IF

     CALL RANDOM_SEED() 

     f_x = 1.0D+0
     DO j = 1,num_thetta
        IF (x(1,j) > Upper_Bound(j) .OR. x(1,j) < Lower_Bound(j) ) THEN
             f_xi = 0.0D+0
             EXIT
        END IF
     END DO
     
     equilibrium = 0
     signout     = 1 
     xi          = x(1,:)

     Par1 = 0.0D+00
     Par2 = 0.0D+00
     IF ( SUM( v_beta) /= 0 ) THEN
        Par1(v_beta) = (( Pmean(v_beta)/Pstd(v_beta) )**2 *( 1.0D+0/Pmean(v_beta) - 1.0D+0 )-1.0D+0)*Pmean(v_Beta)
        Par2(v_beta) = (1.0D+0/Pmean(v_beta)-1.0D+0)*Par1(v_beta)
     END IF
     IF ( SUM(v_gama) /=0 ) THEN
        Par2(v_gama)  = Pstd(v_gama)**2/(Pmean(v_gama) - Lower_Bound(v_gama))
        Par1(v_gama)  = (Pmean(v_gama) - Lower_Bound(v_gama))/Par2(v_gama)
     END IF
     IF ( SUM(v_normal) /= 0 ) THEN
        Par1(v_normal)= Pmean(v_normal)
        Par2(v_normal)= Pstd(v_normal)
     END IF
     IF (SUM(v_uniform) /= 0) THEN
        Par1(v_uniform) =  Lower_Bound(v_uniform)
        Par2(v_uniform) =  Upper_Bound(v_uniform)
     END IF
     IF ( SUM(v_igama) /=0 ) THEN
        Par1(v_igama) = 2.0D+0+(Pmean(v_igama)-Lower_Bound(v_igama)/(Pstd(v_igama)))**2
        Par2(v_igama) = (Pmean(v_igama) - Lower_Bound(v_igama))*(Par1(v_igama)-1.0D+0)
     END IF

     write(*,*) 'Par1', Par1
     write(*,*) 'Par2', Par2

     CALL RANDOM_NUMBER(rand_temp)
     DO i = 1, Ndraws
          IF (MOD(i,10000) == 0) THEN
                index_acc = index_acc + 1
                acceptance_all(index_acc) = acceptance
                write(*,*) 'i=', i ,'Dist', Ln_x(i),'a ac',real( acceptance)/real(i)*real(100),'m ac', real(acceptance_all(index_acc)-acceptance_all(index_acc-1))/real(100)
! write(*,*) 'Ln_xhat', Ln_xhat
write(*,*) 'x=',x(i,:)
           END IF

           equilibrium = 0
           DO WHILE (equilibrium /= 1 .OR. signout /=0)
              DO j = 1, num_thetta
                 randn(j) =  RANDOM_NORMAL()
              END DO
              xi   = x(i,:) + MATMUL(sig_xi*randn,D)
              f_xi = 1.0D+0
              DO j = 1,num_thetta
                 IF (xi(j) > Upper_Bound(j) .OR. xi(j) < Lower_Bound(j) ) THEN
                    f_xi = 0.0D+0
                    EXIT
                 END IF
              END DO
              DO WHILE (f_xi == 0.0D+00)
                 DO j = 1, num_thetta
                    randn(j) =  RANDOM_NORMAL()
                 END DO
                 xi   = x(i,:) + MATMUL(sig_xi*randn,D)
                 f_xi = 1.0D+0
                 DO j = 1,num_thetta
                    IF (xi(j) > Upper_Bound(j) .OR. xi(j) < Lower_Bound(j) ) THEN
                       f_xi = 0.0D+0
                       EXIT
                    END IF
                 END DO
              END DO
              CALL OBJECTIVE(Ln_xi, xi, MOMdata, invSIG_data, num_MOM**2*(nlags+1)-num_MOM, equilibrium, signout)
              IF (equilibrium /=1 .OR. signout /=0 ) THEN
                 rho = 0.0D+00
              ELSE
                 f_xi = fprior(x(i,:),xi, Par1, Par2 )
                 rho        = max( 0.0D+00, min( exp(multiplier*(Ln_xi(1,1) - Ln_x(i)))*f_xi, 1.0D+00 ) )
              END IF
              rand = rand_temp(i)
!write(*,*) 'L=',Ln_x(i),' Li=', Ln_xi, 'RHO = ', rho, '  rand =', rand, ' f_xi =',f_xi
              IF (rand <= rho(1,1)) THEN  
                 x(i+1,:)       = xi   
                 Ln_x(i+1)      = Ln_xi(1,1)              
                 acceptance     = acceptance + 1     
                 f_x            = f_xi 
              ELSEIF ( rand > rho(1,1) ) THEN
                 x(i+1,:)       = x(i,:)   
                 Ln_x(i+1)      = Ln_x(i)
              END IF
          END DO
      END DO

  write(filename, '(A)') 'xDist_'//trim(arg1)//'_'//'mom.txt'
  OPEN(unit = 22 ,file = filename, action = "write")
  DO j = 1,Ndraws/SaveEach+1
     WRITE(22,*) Ln_x((j-1)*SaveEach+1)
  END DO
  CLOSE(22)

        OPEN(unit = 24 ,file = "xUpper.txt", action = "write")
        OPEN(unit = 25 ,file = "xLower.txt", action = "write")
        OPEN(unit = 23 ,file = "xPar1.txt", action = "write")
        OPEN(unit = 22 ,file = "xPar2.txt", action = "write")
        OPEN(unit = 26 ,file = "xmean.txt", action = "write")
        OPEN(unit = 27 ,file = "xstd.txt", action = "write")
        DO j = 1, num_thetta
           WRITE(22,*) Par2(j)
           WRITE(23,*) Par1(j)
           WRITE(24,*) Upper_bound(j)
           WRITE(25,*) Lower_bound(j)
           WRITE(26,*) Pmean(j)
           WRITE(27,*) Pstd(j)
        END DO
        CLOSE(22)
        CLOSE(23)
        CLOSE(24)
        CLOSE(25)
        CLOSE(26)
        CLOSE(27)


       OPEN(unit = 22 ,file = "xv_beta.txt", action = "write")
        DO j = 1,size_beta
           WRITE(22,*) v_beta(j)
        END DO
        CLOSE(22)

        OPEN(unit = 22 ,file = "xv_gama.txt", action = "write")
        DO j = 1,size_gama
           WRITE(22,*) v_gama(j)
        END DO
        CLOSE(22)

        OPEN(unit = 22 ,file = "xv_igama.txt", action = "write")
        DO j = 1, size_igama
           WRITE(22,*) v_igama(j)
        END DO
        CLOSE(22)

        OPEN(unit = 22 ,file = "xv_uniform.txt", action = "write")
        DO j = 1, size_uniform
           WRITE(22,*) v_uniform(j)
        END DO
        CLOSE(22)

        OPEN(unit = 22 ,file = "xv_normal.txt", action = "write")
        DO j = 1, size_normal
           WRITE(22,*) v_normal(j)
        END DO
        CLOSE(22)

  IF (1 == 1) THEN
      DO i = 1, num_thetta
           write( arg ,'(I8)') i
           arg = adjustl(arg)
           write(filename, '(A)') 'x'//trim(arg)//'_'//trim(arg1)//'_mom.txt'
           OPEN(unit = i, file = filename, action = "write")
           DO j = 1, whentosave/Saveeach+1
              WRITE(i,'(ES11.4)') x( Ndraws - whentosave +  (j - 1)*Saveeach + 1, i  )
           END DO
           CLOSE(i)
      END DO
  END IF

        WRITE(33,*) 'sig_xi = ', sig_xi
        WRITE(33,*) 'size_data = ', size_data
        WRITE(33,*) 'Ndraws = ',Ndraws
        WRITE(33,*) 'INFO_HESS = ', INFO_HESS
        WRITE(33,*) 'multiplier', multiplier
        WRITE(33,*) 'acceptance', acceptance
        WRITE(33,*) 'acceptance_all', acceptance_all
        WRITE(33,*) 'x0 = ', x(1,:)
        WRITE(33,*) 'Ln_x0', Ln_x(1)
        WRITE(33,*) 'Ln_xlast=',Ln_x(Ndraws+1)
        CLOSE(33)

  DO j = 1,num_thetta
     mean_est(1,j) =  sum(x(Ndraws/2:Ndraws,j))/dble(Ndraws/2)
  END DO

 IF ( enter_i == 1 .OR. enter_i == 2) THEN 
    OPEN(unit = 22 ,file = "xmean_est.txt", action = "write")
    DO j = 1, num_thetta
      WRITE(22,*) mean_est(1,j)
    END DO
    CLOSE(22)
  END IF

!  CALL sample_mom( MOMdata, invSIG_data1, Vmom, num_MOM**2*(nlags+1)-num_MOM , mean_est, data_sample ) 
  CALL GRADIENT_MOM(mdot, mean_est,  num_MOM**2*(nlags+1)-num_MOM, grdh)

  DO i = 1, Ndraws+1
     x_mean(i,:) = x(i,:) - mean_est(1,:)
  END DO
  write(*,*) 'mean_est = ', mean_est

  Vg = dble(size_data)*MATMUL(TRANSPOSE(x_mean(Ndraws/2+2:Ndraws+1,:)), x_mean(Ndraws/2+2:Ndraws+1,:))/dble(Ndraws/2)
  Vtheta = dble(multiplier**2)*MATMUL(MATMUL(MATMUL(MATMUL(MATMUL( MATMUL( Vg, TRANSPOSE(mdot)), invSIG_data), Vmom ), invSIG_data ), mdot), Vg)
!  Vtheta = dble(multiplier**2)*MATMUL(MATMUL(MATMUL(MATMUL( Vg, TRANSPOSE(mdot)), invSIG_data),  mdot), Vg)
 
 write(filename, '(A)') 'xVg_'//trim(arg1)//'.txt'
  OPEN(unit = 22 ,file = filename, action = "write")
  DO j = 1, num_thetta
     DO i = 1, num_thetta
      WRITE(22,*) Vg(i,j)
     END DO
  END DO
  CLOSE(22)
  write(filename, '(A)') 'xVtheta_'//trim(arg1)//'.txt'
  OPEN(unit = 22 ,file = filename, action = "write")
  DO j = 1, num_thetta
     DO i = 1, num_thetta
      WRITE(22,*) Vtheta(i,j)
     END DO
  END DO
  CLOSE(22)

write(*,*) 'Vg'
DO i = 1,num_thetta
   write(*,*) Vg(i,i)
END DO

write(*,*) 'Vtheta'
DO i = 1,num_thetta
   write(*,*) Vtheta(i,i)
END DO

write(*,* ) 'mult', multiplier
write(*,*) 'enter_i', enter_i
write(*,*) 'sig', sig_xi
END IF
END IF
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) function fprior(Theta,Thetai, Par1, Par2 )
    IMPLICIT NONE
    REAL(8), PARAMETER   :: pi = 3.14159265358979
    REAL(8), DIMENSION(num_thetta), INTENT(IN) :: Par1, Par2, Theta, Thetai
    REAL(8)              :: f,  y_igama(size_igama), y_igamai(size_igama), y_beta(size_beta), y_gama(size_gama), y_normal(size_normal),  y_betai(size_beta),y_gamai(size_gama), y_normali(size_normal)
    INTEGER              :: i

    f = 0.0D+00
    DO i = 1,num_thetta
       IF (Thetai(i) > Upper_Bound(i) .OR. Thetai(i) < Lower_Bound(i) ) THEN
          RETURN
       END IF
    END DO
    f = 1.0D+00
    IF ( SUM( v_beta) /= 0 ) THEN
       y_beta  = ( Theta(v_beta) - Lower_Bound(v_beta) )/( Upper_Bound(v_beta) - Lower_Bound(v_beta) )
       y_betai = ( Thetai(v_beta) - Lower_Bound(v_beta) )/( Upper_Bound(v_beta) - Lower_Bound(v_beta) )
       f      = f * PRODUCT((y_betai/y_beta)**(Par1(v_beta)-1.0D+0)*((1.0D+0-y_betai)/(1.0D+0-y_beta))**(Par2(v_beta)-1.0D+0))
!       write(*,*) 'BETA' ,f
    END IF

    IF ( SUM( v_gama) /= 0 ) THEN
       y_gama = Theta(v_gama)- Lower_Bound(v_gama)
       y_gamai = Thetai(v_gama)- Lower_Bound(v_gama)
       f     = f * PRODUCT( (y_gamai/y_gama)**( Par1(v_gama)-1.0D+0)*exp(-(y_gamai-y_gama)/Par2(v_gama)))
!   write(*,*) 'GAMA' ,f
    END IF
!write(*,*) 'GAMA',f 
    IF ( SUM(v_normal) /= 0 ) THEN
       y_normal  = Theta(v_normal)
       y_normali = Thetai(v_normal)
       f      = f* PRODUCT( exp(-y_normali+y_normal )**2 / (2.0D+0*pi*Par2(v_normal)**2)  / Par2(v_normal) )
!       write(*,*) 'NORMAL', f
    END IF
    IF ( SUM( v_igama) /= 0 ) THEN
       y_igama  = Theta(v_igama)- Lower_Bound(v_igama)
       y_igamai = Thetai(v_igama)- Lower_Bound(v_igama)
       f     = f * PRODUCT((y_igamai/y_igama)**( -Par1(v_igama)-1.0D+0)*exp(-Par2(v_igama)*(1.0D+0/y_igamai-1.0D+0/y_igama)) )
!       write(*,*) 'iGAMA' ,f
    END IF
    fprior  = f
  END FUNCTION fprior
END PROGRAM mainMCMC
