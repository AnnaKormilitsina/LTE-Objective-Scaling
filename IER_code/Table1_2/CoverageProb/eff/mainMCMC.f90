PROGRAM mainMCMC

USE MPARAMS
USE FUNCTIONS
USE MODEL_SOL

IMPLICIT NONE

INTEGER, PARAMETER :: Ndraws = 200000, Whentosave = Ndraws, SaveEach = 100, num_rep = 1000
REAL(8)            :: sig_xi 
REAL(8), DIMENSION(num_thetta) :: xi, randn, Par1, Par2

REAL(8) :: x(Ndraws+1,num_thetta),  Ln_x(Ndraws+1), xstart(num_thetta)
INTEGER :: i, j, k, index_acc, equilibrium, signout,  acceptance, acceptance_all(Ndraws/10000 )
REAL(8) :: Ln_xi(1,1),  f_xi, f_x, rho(1,1)
REAL    :: rand, rand_temp(Ndraws)

REAL(8) :: data_mean(size_data,1),emp_data(size_data,1), MOMdata(num_MOM,1), invSIG_data(num_MOM,num_MOM), Vmom(num_MOM,num_MOM)
REAL(8) :: D(num_thetta,num_thetta)

REAL(8) ::  grdh(num_thetta)

INTEGER :: INFO_HESS, INFO
REAL(8) :: multiplier 

REAL(8) :: Ln_xhat(1,1), meanx(num_rep,num_thetta), mmean_x(1,num_thetta), mean_m(num_rep, num_thetta) , stdx(num_rep, num_thetta),stdg(num_rep, num_thetta), Vg(num_thetta,num_thetta), x_mean(Ndraws+1,num_thetta),Vtheta(num_thetta,num_thetta), Vm(num_thetta,num_thetta), std_mean(1,num_thetta)
REAL(8) :: th_mean(num_rep,num_thetta), meanth(1,num_thetta), th_cov(num_thetta,num_thetta)

character(30) ::filename
CHARACTER(8)  :: arg, arg1, arg2, arg3, arg4
INTEGER       :: enter_i, enter_j, nsig_xi, count(2), countg(2)

! size_data
write( arg3 ,'(I8)') size_data
arg3 = adjustl(arg3)

!input 1, if 1 - starting value is Pmean and MCMC mean value is stored in mean_est.txt as a result, otherwise from mean_est.txt
CALL get_command_argument(1,arg1)
read(arg1,FMT='(I5)') enter_i
arg1 = adjustl(arg1)

 ! input 2 = value of multiplier (scaling parameter)
CALL get_command_argument(2,arg2)
read(arg2,FMT='(I5)') enter_j
multiplier = dble(enter_j)

! input 3 is the value of SIGMA, scaling of the proposal distribution to achieve 30-40% acceptance rate
CALL get_command_argument(3,arg4)
read(arg4,FMT='(I5)') nsig_xi

sig_xi = 0.0D+0
DO i = 1, nsig_xi
  sig_xi = sig_xi + 0.001D+0
END DO
write(*,*) 'sig_xi = ', sig_xi, nsig_xi

! Starting value of MCMC chain
!=========================================================================================
IF (enter_i == 1) THEN
  x(1,:) = Pmean
ELSE!IF (enter_i == 2) THEN
  x(1,:) = Pmean
  OPEN(unit = 1 ,file = 'xmean_est.txt', status = "old", action = "read")
  READ( 1, *, end = 7 ) xstart
  7 continue
  CLOSE(1)
  x(1,:) = xstart
END IF

write(*,*) 'x1 = ', x(1,:)

count  = 0
countg = 0
!==========================================================================================
DO k = 1,num_rep
  write(*,*) 'k = ',k
  IF ( 1 == 1) THEN
    ! Generate data
    CALL RANDOM_SEED()
    DO i = 1,size_data
      emp_data(i,1) = 0.5D+0*RANDOM_NORMAL() + 0.5D+0
    END DO

    MOMdata(1,1) = sum(emp_data)/dble(size_data) 
    DO i = 1,size_data
       data_mean(i,1) = (emp_data(i,1) - MOMdata(1,1))**2
    END DO
    MOMdata(2,1) = sum(data_mean)/dble(size_data) 

    invSIG_data(1,1)  = MOMdata(2,1)
    DO i = 1,size_data
       data_mean(i,1) = (emp_data(i,1) - MOMdata(1,1))**3
    END DO
    invSIG_data(1,2)  = sum(data_mean)/dble(size_data)
    invSIG_data(2,1)  = invSIG_data(1,2) 
    DO i = 1,size_data
       data_mean(i,1) = (emp_data(i,1) - MOMdata(1,1))**4
    END DO
    invSIG_data(2,2) = sum(data_mean)/dble(size_data) - MOMdata(2,1)**2

    CALL MATINV(invSIG_data, num_MOM, INFO)

  END IF

  !!! START MCMC
  IF (1 == 1) THEN
    Vmom = invSIG_data
    CALL MATINV(Vmom, num_MOM, INFO)

    ! Identity weighting matrix
!    invSIG_data(1,1)= 1.0D+0
!    invSIG_data(2,1)= 0.0D+0
!    invSIG_data(1,2)= 0.0D+0
!    invSIG_data(2,2)= 1.0D+0
    !=========================================================================================
    CALL RANDOM_SEED() 
    CALL OBJECTIVE(Ln_xi, x(1,:), MOMdata, invSIG_data)
    IF ( 1 == 1 ) THEN

      Ln_x(1)  = Ln_xi(1,1)

      acceptance = 0
      acceptance_all = 0
      index_acc  = 1

      grdh = 1.0D-04*x(1,:)
      CALL HESSIAN_SYM(D, x(1,:), grdh, emp_data, invSIG_data, INFO_HESS)
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
      xi       = x(1,:)

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
       IF (MOD(i,10000) == 0 .AND. i < 50000) THEN
              index_acc = index_acc + 1
             acceptance_all(index_acc) = acceptance
            write(*,*) 'i=', i ,'Dist', Ln_x(i),'a ac',real( acceptance)/real(i)*real(100),'m ac', real(acceptance_all(index_acc)-acceptance_all(index_acc-1))/real(100)
        END IF
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
       CALL OBJECTIVE(Ln_xi, xi, MOMdata, invSIG_data)
       f_xi = fprior(x(i,:),xi, Par1, Par2 )
       rho        = max( 0.0D+00, min( exp(multiplier*(Ln_xi(1,1) - Ln_x(i)))*f_xi, 1.0D+00 ) )
       rand = rand_temp(i)
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

    meanx(k,1) = sum(x(NDraws/2+2:Ndraws+1,1))/dble(NDraws/2)
    meanx(k,2) = sum(x(NDraws/2+2:Ndraws+1,2))/dble(NDraws/2)

    DO i = 1, NDraws + 1
       x_mean(i,:) = x(i,:) - meanx(k,:) 
    END DO 
    Vg = dble(size_data)* MATMUL(TRANSPOSE(x_mean(NDraws/2+2:Ndraws+1,:)),x_mean(NDraws/2+2:Ndraws+1,:)) /dble(NDraws/2)
!    Vtheta = dble(multiplier)**2*MATMUL(MATMUL(Vg,invSIG_data),Vg)
    Vtheta = dble(multiplier)**2*MATMUL(MATMUL(MATMUL(MATMUL(Vg,invSIG_data), Vmom), invSIG_data), Vg)

!    write(*,*) (dble(size_data))*(abs(meanx(k,1) - Pmean(1))/1.96D+0)**2 ,  (Vtheta(1,1))
!    write(*,*) (dble(size_data))*(abs(meanx(k,2) - Pmean(2))/1.96D+0)**2 ,  (Vtheta(2,2))

    DO i = 1,num_thetta
       stdx(k,i) = sqrt(Vtheta(i,i)) 
    END DO

    DO i = 1,num_thetta
       stdg(k,i) = sqrt(Vg(i,i))
    END DO

    IF ( abs(meanx(k,1) - Pmean(1))  < 1.96D+0/sqrt(dble(size_data))*stdx(k,1) ) THEN
       count(1) = count(1)+1
    END IF
    IF ( abs(meanx(k,1) - Pmean(1))  < 1.96D+0/sqrt(dble(size_data))*stdg(1,1) ) THEN
       countg(1) = countg(1)+1
    END IF

    IF (abs(meanx(k,2) - Pmean(2))< 1.96D+0/sqrt(dble(size_data))*stdx(k,2) ) THEN
       count(2) = count(2)+1
    END IF
    IF (abs(meanx(k,2) - Pmean(2))< 1.96D+0/sqrt(dble(size_data))*stdg(2,2) ) THEN
       countg(2) = countg(2)+1
    END IF
   END IF
 END IF
END DO ! k = 1,100

meanth(1,1) = sum(meanx(:,1))/dble(num_rep)
meanth(1,2) = sum(meanx(:,2))/dble(num_rep)
th_cov = 0.0D+0
DO i = 1, num_rep
   th_mean(1,:) = meanx(i,:) - meanth(1,:)
   th_cov = th_cov + MATMUL( TRANSPOSE(th_mean ) , th_mean )
END DO

th_cov = dble(size_data)* th_cov/dble(num_rep)

write(*,*) 'std1', sqrt(th_cov(1,1) ), sqrt(th_cov(2,2))

write(*,*) 'coverage probabilities', count
write(*,*) 'coverage probabilities g', countg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*) 'mean of std'
  DO j = 1,num_thetta
     write(*,*)   sum(stdx(:,j))/dble(num_rep)
  END DO

DO j = 1,num_thetta
     mmean_x(1,j) =  sum(meanx(:,j))/dble(num_rep)
END DO

DO i = 1, num_rep
     mean_m(i,:) = meanx(i,:) - mmean_x(1,:)
END DO

Vm = MATMUL(TRANSPOSE(mean_m), mean_m)/dble(num_rep)
DO i = 1,num_thetta
   std_mean(1,i) = sqrt(dble(size_data) * Vm(i,i))
END DO

write(*,*) 'std of mean', std_mean
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(filename, '(A)') 'xcovprob_'//trim(arg2)//'.txt'
OPEN(unit = 22,file = filename,action = "write")
  DO j = 1, num_thetta
     WRITE(22,*) dble(count(j))*100.0D+0/dble(num_rep)
  END DO
CLOSE(22)

write(filename, '(A)') 'xcovprobg_'//trim(arg2)//'.txt'
OPEN(unit = 22,file = filename,action = "write")
  DO j = 1, num_thetta
     WRITE(22,*) dble(countg(j))*100.0D+0/dble(num_rep)
  END DO
CLOSE(22)

write(filename, '(A)') 'xmeanx_'//trim(arg2)//'.txt'
OPEN(unit = 22,file = filename,action = "write")
DO i = 1,num_thetta
  DO j = 1, num_rep
     WRITE(22,*) meanx(j,i)
  END DO
END DO
CLOSE(22)

write(filename, '(A)') 'xstdx_'//trim(arg2)//'.txt'
OPEN(unit = 22,file = filename,action = "write")
DO i = 1,num_thetta  
   DO j = 1, num_rep
     WRITE(22,*) stdx(j,i)
  END DO
END DO
CLOSE(22)

write(filename, '(A)') 'xstdg_'//trim(arg2)//'.txt'
OPEN(unit = 22,file = filename,action = "write")
DO i = 1,num_thetta
   DO j = 1, num_rep
     WRITE(22,*) stdg(j,i)
  END DO
END DO
CLOSE(22)

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
    END IF
    IF ( SUM( v_gama) /= 0 ) THEN
       y_gama = Theta(v_gama)- Lower_Bound(v_gama)
       y_gamai = Thetai(v_gama)- Lower_Bound(v_gama)
       f     = f * PRODUCT( (y_gamai/y_gama)**( Par1(v_gama)-1.0D+0)*exp(-(y_gamai-y_gama)/Par2(v_gama)))
    END IF
    IF ( SUM(v_normal) /= 0 ) THEN
       y_normal  = Theta(v_normal)
       y_normali = Thetai(v_normal)
       f      = f* PRODUCT( exp(-y_normali+y_normal )**2 / (2.0D+0*pi*Par2(v_normal)**2)  / Par2(v_normal) )
    END IF
    IF ( SUM( v_igama) /= 0 ) THEN
       y_igama  = Theta(v_igama)- Lower_Bound(v_igama)
       y_igamai = Thetai(v_igama)- Lower_Bound(v_igama)
       f     = f * PRODUCT((y_igamai/y_igama)**( -Par1(v_igama)-1.0D+0)*exp(-Par2(v_igama)*(1.0D+0/y_igamai-1.0D+0/y_igama)) )
    END IF
    fprior  = f
  END FUNCTION fprior
END PROGRAM mainMCMC
