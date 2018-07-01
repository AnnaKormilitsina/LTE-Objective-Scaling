PROGRAM mainMCMC

USE EMPIRICAL
USE MPARAMS
USE FUNCTIONS
USE MODEL_SOL
USE MODEL

IMPLICIT NONE

INTEGER, PARAMETER :: Ndraws = 500000, Whentosave = Ndraws, SaveEach = 100, num_rep = 100
REAL(8) ::  sig_xi 
REAL(8), DIMENSION(num_thetta) :: xi, randn, Par1, Par2

!!         B     KAPA    ALFA    ALFATIL  CHI      CHITIL      ALFA_R1    ALFA_PI  ALFA_Y   RHO  BETTA  THETA 
REAL(8) :: x(Ndraws+1,num_thetta),  Ln_x(Ndraws+1), xstart(num_thetta)
INTEGER :: i, j, k, index_acc, equilibrium, signout,  acceptance, acceptance_all(Ndraws/10000 )
REAL(8) :: Ln_xi(1,1),  f_xi, f_x, rho(1,1)
REAL    :: rand, rand_temp(Ndraws)

REAL(8) :: MOMdata(num_MOM**2*(nlags+1) - num_MOM,1), invSIG_data(num_MOM**2*(nlags+1)-num_MOM,num_MOM**2*(nlags+1)-num_MOM),   Vmom( num_MOM**2*(nlags+1) - num_MOM,num_MOM**2*(nlags+1) - num_MOM ), W0(num_thetta,num_MOM**2*(nlags+1) - num_MOM), W( num_MOM**2*(nlags+1) - num_MOM,num_MOM**2*(nlags+1) - num_MOM )
REAL(8) :: D(num_thetta,num_thetta), mdot(num_MOM**2*(nlags+1) - num_MOM,num_thetta), Vg(num_thetta,num_thetta) , mean_est(num_rep,num_thetta),  std_est(num_rep,num_thetta), mmean_est(1,num_thetta), mean_m(num_rep,num_thetta), std_mean(1,num_thetta)
REAL(8) ::  mdot1(num_MOM**2*(nlags+1) - num_MOM,num_thetta),  mdot2(num_MOM**2*(nlags+1) - num_MOM,num_thetta),  mdot3(num_MOM**2*(nlags+1) - num_MOM,num_thetta)
REAL(8) :: Vgmdot(num_thetta,num_MOM**2*(nlags+1) - num_MOM), VgmdotT(num_MOM**2*(nlags+1) - num_MOM,num_thetta)
REAL(8) :: gx(sizey,sizex), hx(sizex,sizex), grdh(num_thetta), Vtheta(num_thetta,num_thetta),  x_mean(Ndraws+1,num_thetta), Vm(num_thetta,num_thetta)

INTEGER :: INFO_HESS, INFO
REAL(8) :: multiplier

REAL(8) :: shock_place(sizex), shock(sizex), x0(sizex), data(sizef),  err(1,num_MOM)
REAL(8) :: data_sample(size_data,num_MOM)

character(30) ::filename
CHARACTER(8)  :: arg, arg1, arg2, arg3
INTEGER       :: enter_i, enter_j, nsig_xi, enter_k, count(num_thetta)

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

write( arg ,'(I8)') size_data
arg = adjustl(arg)

count = 0

DO k = 1, num_rep
!write(*,*) '',k
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

  CALL gx_hx(gx,hx, x(1,:), signout, equilibrium)

  IF ( 1 == 1) THEN

    CALL RANDOM_SEED()
    shock_place       = 0.0D+00

    x0 = 0.0D+00
    DO i = 1, size_data

      CALL gx_hx(gx,hx, x(1,:), signout, equilibrium)

      shock_place(sth)  = (x(1,5))*0.01D+0
      shock_place(seps) = (x(1,6))*0.01D+0
      shock_place(sgam) = (x(1,7))*0.01D+0

      DO j = 1, sizex
         shock(j)  = shock_place(j)*RANDOM_NORMAL()
      END DO
      x0                       = MATMUL(hx,x0) + shock! shock_place*RANDOM_NORMAL()
      data(1 : sizey )         = MATMUL(gx,x0)
      data( sizey + 1 : sizey+sizex) = x0
      data_sample(i,1:num_MOM) = dble(data( (/ nr, noutput, npai /)))
    END DO
     
!     write(filename, '(A)') 'mdata_'//trim(arg)//'.txt'
!     OPEN(unit = 22,file = filename,action = "write")
!     DO i = 1, num_MOM
!        DO j = 1, size_data
!          WRITE(22,*)  data_sample(j,i)
!        END DO
!     END DO
!     CLOSE(22)

  ELSE


     write(filename, '(A)') 'mdata_'//trim(arg)//'.txt'
     OPEN(unit = 2 ,file = filename, status = "old", action = "read")
     READ( 2, *, end = 19 ) data_sample
     19 continue
     CLOSE(2)

  END IF

!=========================================================================================
   CALL RANDOM_SEED() 
   x(1,:) = xstart

   CALL sample_mom(MOMdata,invSIG_data, Vmom, num_MOM**2*(nlags+1)-num_MOM , x(1,:), data_sample)

! invSIG_data = Vmom
!   DO i = 1, num_MOM**2*(nlags+1)-num_MOM
!     DO j = 1, num_MOM**2*(nlags+1)-num_MOM
!        IF (i == j) THEN
!        ELSE 
!           invSIG_data(i,j) = 0.0D+0
!      END IF
!     END DO
!  END DO
!CALL(invSIG_Data, num_MOM**2*(nlags+1)-num_MOM, INFO)

   CALL OBJECTIVE(Ln_xi, x(1,:), MOMdata, invSIG_data,  num_MOM**2*(nlags+1)-num_MOM, equilibrium, signout)

   write(*,*) 'Ln_xi', Ln_xi

  IF ( 1 == 1) THEN

    Ln_x(1)  = Ln_xi(1,1)

    acceptance = 0
    acceptance_all = 0
    index_acc  = 1

    IF (signout == 1 .OR. equilibrium /= 1 ) THEN
!      WRITE(33,*) '; signout ==' , signout, '; equilibrium = ', equilibrium
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
 
     CALL RANDOM_NUMBER(rand_temp)
     DO i = 1, Ndraws
          IF (MOD(i,10000) == 0 .AND. i<50000) THEN
                index_acc = index_acc + 1
                acceptance_all(index_acc) = acceptance
                write(*,*) 'k=', k, 'i=', i ,'Dist', Ln_x(i),'a ac',real( acceptance)/real(i)*real(100),'m ac', real(acceptance_all(index_acc)-acceptance_all(index_acc-1))/real(100)
!write(*,*) 'x=',x(i,:)
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

  DO j = 1,num_thetta
     mean_est(k,j) =  sum(x(Ndraws/2+2:Ndraws+1,j))/dble(Ndraws/2)
  END DO

  !IF ( enter_i == 1) THEN 
  !  OPEN(unit = 22 ,file = "xmean_est.txt", action = "write")
  !  DO j = 1, num_thetta
  !    WRITE(22,*) mean_est(1,j)
  !  END DO
  !  CLOSE(22)
  !END IF

  CALL GRADIENT_MOM(mdot,  mean_est(k,:),  num_MOM**2*(nlags+1)-num_MOM, grdh)

  DO i = 1, Ndraws+1
     x_mean(i,:) = x(i,:) - mean_est(k,:)
  END DO

  Vg = dble(size_data)*MATMUL(TRANSPOSE(x_mean(Ndraws/2+2:Ndraws+1,:) ), x_mean(Ndraws/2+2:Ndraws+1,:) )/dble(Ndraws/2)
  Vtheta = dble(multiplier)**2*MATMUL(MATMUL(MATMUL(MATMUL(MATMUL( MATMUL( Vg, TRANSPOSE(mdot)), invSIG_data), Vmom ), invSIG_data ), mdot), Vg)

  DO i = 1, num_thetta
      std_est(k,i) = sqrt(Vtheta(i,i))
  END DO
 
  write(*,*) 'std', std_est(k,:)

  DO i = 1, num_thetta
    IF ( abs(mean_est(k,i) - Pmean(i)) < 1.96D+0/sqrt(dble(size_data))*std_est(k,i) ) THEN
       count(i) = count(i)+1
    END IF
  END DO

END IF
END IF

END DO ! k = 1,100

write(*,*) 'mean of std'
DO j = 1, num_thetta
   write(*,*)   sum(std_est(:,j))/dble(num_rep)
END DO

DO j = 1,num_thetta
     mmean_est(1,j) =  sum(mean_est(:,j))/dble(num_rep)
END DO
DO i = 1, num_rep
     mean_m(i,:) = mean_est(i,:) - mmean_est(1,:)
END DO

Vm = dble(size_data)* MATMUL(TRANSPOSE(mean_m), mean_m)/dble(num_rep)

 DO i = 1,num_thetta
   std_mean(1,i) = sqrt(Vm(i,i))
END DO
write(*,*) 'std of mean', std_mean

write(*,*) 'count', count
  
IF ( enter_i == 1) THEN

  write(filename, '(A)') 'xmean_estk_'//trim(arg1)//'.txt'
  OPEN(unit = 22,file = filename,action = "write")
  DO j = 1, num_thetta
       DO i = 1,num_rep
          WRITE(22,*) mean_est(i,j)
      END DO
  END DO
  CLOSE(22)
  END IF

  write(filename, '(A)') 'xstd_estk_'//trim(arg1)//'.txt'
  OPEN(unit = 22,file = filename,action = "write")
    DO j = 1, num_thetta
       DO i = 1,num_rep
          WRITE(22,*) std_est(i,j)
      END DO
    END DO
    CLOSE(22)

  write(filename, '(A)') 'xcount_'//trim(arg1)//'.txt'
  OPEN(unit = 22,file = filename,action = "write")
  DO i = 1, num_thetta
      WRITE(22,*) count(i)
  END DO
  CLOSE(22)

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) function fprior(Theta, Thetai, Par1, Par2 )
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
