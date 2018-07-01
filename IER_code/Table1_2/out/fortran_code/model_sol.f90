MODULE MODEL_SOL

USE MPARAMS
USE MODEL
IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
SUBROUTINE  OBJECTIVE(obj, x, empdata, invSIG_emp)
IMPLICIT NONE
!INTEGER, INTENT(IN)  :: size_data
REAL(8), INTENT(IN)  :: empdata(num_MOM,1), invSIG_emp(num_MOM,num_MOM), x(num_thetta,1)
REAL(8), INTENT(OUT) :: obj(1,1)
REAL(8) :: emp_theory(num_MOM,1), meandata(1,1), emp_x2(size_data,1)
INTEGER              :: i, j

emp_theory(1,1) = empdata(1,1) - x(1,1)
emp_theory(2,1) = empdata(2,1) - x(2,1) 

obj = -0.5D+0*dble(size_data)* MATMUL(MATMUL(TRANSPOSE(emp_theory),invSIG_emp),emp_theory)
!obj = -MATMUL(MATMUL(TRANSPOSE(emp_theory),invSIG_emp),emp_theory)

END SUBROUTINE objective

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MATINV(A, N, INFO)
! Calculates matrix inverse for double floating type variables
! A is the input matrix of dimension N (square), INFO is not zero if there are problems with calculating the inverse.

!(c) Anna Kormilitsina
! Date: June 10, 2009


IMPLICIT NONE
INTEGER, INTENT(IN)    :: N
REAL(8), INTENT(INOUT) :: A(N,N)
INTEGER                :: IPIV(N)
INTEGER, INTENT(OUT)   :: INFO
REAL(8)                :: WORK(N)

CALL  DGETRF(N,N,A,N,IPIV,INFO)
IF (INFO == 0) THEN
   CALL DGETRI( N, A, N, IPIV, WORK, N, INFO )
   IF (INFO /=0) THEN
       write(*,*) 'ERROR: problem with DGETRI in MATINV'
       write(*,*) 'INFO', INFO
       A =  0.0D+00
   END IF 
ELSE 
   write(*,*) 'ERROR : problems with DGETRF in  MATINV'
   write(*,*) 'INFO',INFO
   A = 0.0D+00
END IF
END SUBROUTINE MATINV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SMATINV(A,N, INFO)

! Calculates matrix inverse for single floating type variables
! A is the input matrix of dimension N (square), INFO is not zero if there are problems with calculating the inverse.

!(c) Anna Kormilitsina
! Date: June 10, 2009

IMPLICIT NONE
INTEGER, INTENT(IN) :: N
INTEGER, INTENT(OUT):: INFO
REAL, INTENT(INOUT) :: A(N,N)
INTEGER             :: IPIV(N)
REAL(8)             :: WORK(N)
REAL(8)             :: AA(N,N)

AA = dble(A)
CALL  DGETRF(N,N,AA,N,IPIV,INFO)

IF (INFO /=0) THEN
!    write(*,*) 'Problem with DGETRF in SMATINV'
    AA = 0.0D+00
ELSE
    CALL DGETRI( N, AA, N, IPIV, WORK, N, INFO )
    IF (INFO /= 0) THEN
!         write(*,*) 'Problem with DGETRI in SMATINV'
         AA = 0.0D+00
    END IF
END IF
A = sngl(AA)
END SUBROUTINE SMATINV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MATDET(A, N, det)
IMPLICIT NONE
INTEGER, INTENT(IN)    :: N
REAL(8), INTENT(INOUT) :: A(N,N)
REAL(8), INTENT(OUT)   :: det
INTEGER :: IPIV(N), INFO, i , D

CALL  DGETRF(N,N,A,N,IPIV,INFO)
det  = 1.0D+00

IF (INFO < 0)THEN 
   WRITE(*,*) 'Warning: Determinant can not be computed'
   RETURN
END IF

D = 0
DO i = 1,N
    IF ( IPIV(i) == i ) THEN
       D = D + 1
    END IF
    det = det*A(i,i)
END DO
det = det*(-1.0D+00)**(N-D)
END SUBROUTINE MATDET
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HESSIAN_SYM(D, x, grdh, empirical, invSIG, output)

IMPLICIT NONE
!INTEGER, INTENT(IN)  :: size_data
INTEGER, INTENT(OUT) :: output
REAL(8), INTENT(IN)  :: x(num_thetta)
REAL(8), INTENT(IN)  :: empirical(size_data,1), invSIG(num_MOM,num_MOM),  grdh(num_thetta)
REAL(8), INTENT(OUT) :: D(num_thetta,num_thetta)
COMPLEX(16)          :: Eigenval(num_thetta,num_thetta), Eigenvec(num_thetta,num_thetta)
INTEGER              :: INFO, i , j
REAL(8)              :: WR(num_thetta), WI(num_thetta), VL(1,num_thetta), VR(num_thetta,num_thetta), WORK(4*num_thetta)

output = 1
CALL HESS(D, x, grdh, empirical, invSIG)
IF ( SUM(D) == 0.0D+00) THEN
   output = 0
   WRITE(*,*) 'Warning : Hessian was not calculated in HESS WHEN CALLING  HESSIAN_SYM '
   RETURN
END IF

CALL MATINV(D, num_thetta, INFO)
IF (INFO /=0 ) THEN
   write(*,*) 'Warning: Problem with matrix inverse in HESSIAN SYM!'
   D = 0.0D+00
   output = 0
   RETURN
END IF
CALL DGEEV( 'N', 'V', num_thetta, -D, num_thetta, WR, WI, VL, 1, VR, num_thetta, WORK, 4*num_thetta, INFO )
IF (INFO /= 0) THEN
   WRITE(*,*) 'Warning: Problem with DGGEV in HESSIAN_SYM!'
   D = 0.0D+00
   output = 0
   RETURN
END IF
        
Eigenval = cmplx(0.0D+00, 0.0D+00)
DO i = 1,num_thetta
  IF (WI(i) > 0.0D+00 ) THEN 
     Eigenvec(:,i)     = cmplx(VR(:,i),  VR(:,i+1))
     Eigenvec(:,i+1)   = cmplx(VR(:,i), -VR(:,i+1))
     Eigenval(i,i)     = cmplx(WR(i), WI(i))
     Eigenval(i+1,i+1) = cmplx(WR(i+1),WI(i+1))
  ELSEIF (WI(i) == 0.0D+00 ) THEN
     Eigenvec(:,i) = cmplx(VR(:,i), 0.0D+00)
     Eigenval(i,i) = cmplx(WR(i),0.0D+00)
  END IF
END DO

D = REAL(MATMUL(MATMUL(Eigenvec,abs(Eigenval)),TRANSPOSE(Eigenvec)))

CALL DPOTRF('U', num_thetta, D,  num_thetta , INFO ) ! Cholesky Decomposition
IF (INFO /= 0) THEN
     WRITE(*,*) 'Warning: Problem with DPOTRF in HESSIAN_SYM!,INFO=',info
     D = 0.0D+00
     output = 0
     RETURN
END IF

DO i = 2,num_thetta
  DO j = 1, i-1
    D(i,j) = 0.0D+00
  END DO
END DO

END SUBROUTINE HESSIAN_SYM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HESS(hessian, x0, grdh, empirical, invSIG)
IMPLICIT NONE
!INTEGER, INTENT(IN)  :: size_data
REAL(8), INTENT(OUT) :: hessian(num_thetta,num_thetta)
REAL(8), INTENT(IN)  :: x0(num_thetta), empirical(num_MOM,1), invSIG(num_MOM, num_MOM), grdh(num_thetta)
REAL(8) ::  fune1(1,1), fune2(1,1), fune3(1,1), fune4(1,1), eye(num_thetta,num_thetta), ones(num_thetta), stps, dh(num_thetta), ax0(num_thetta), dax0(num_thetta), hlp(num_thetta), xdh(num_thetta), ee(num_thetta,num_thetta), ee_i(num_thetta), ee_j(num_thetta)
INTEGER :: i, j, k

ones = 1.0D+00

stps = EPSILON(1.0D+00)**(1.0/3.0)
! eps: floating point relative accuracy or machine precision: 2.22e-16
! stps: step size recommended by Dennis and Schnabel: 6.006e-6

!CALL OBJECTIVE( fune1, x0, empirical, invSIG, size_data, equilibrium1, signout1 )
!CALL gx_hx(gx,hx,x0, signout1, equilibrium1)

!IF (equilibrium1 /=1 .OR. signout1 /=0 ) THEN
!   WRITE(*,*) 'Error in HESS : objective can not be calculated'
!   return
!END IF 

! ** Computation of stepsize (dh)
IF ( sum(grdh) > 0.0D+00 ) THEN
    dh = grdh
ELSE
    ax0 = ABS(x0)
    DO j = 1,num_thetta
       hlp(j) = MAX( x0(j),1.0D-02 )
    END DO
    IF ( SUM(x0) /= 0.0D+00 ) THEN
        dax0 = x0/ax0
        dh  = stps * hlp * x0/ABS(x0)
    ELSE
        dh  = stps * hlp      
    END IF
END IF

ee = 0.0D+00
DO k = 1,num_thetta
  ee(k,k) = dh(k)
END DO

hessian = 0.0D+00

DO i = 1, num_thetta
   DO j = i, num_thetta

      ee_i = ee(:,i)
      ee_j = ee(:,j)

     ! equilibrium1 = 3
     ! equilibrium2 = 3
     ! equilibrium3 = 3
     ! equilibrium4 = 3
      !!signout1     = 2
      !signout2     = 2
      !signout3     = 2
      !signout4     = 2

    !  DO WHILE ( equilibrium1 /= 1 .OR. equilibrium2 /= 1 .OR. equilibrium3 /= 1 .OR. equilibrium4 /= 1 .OR. signout1 /= 0 .OR. signout2 /= 0 .OR. signout3 /= 0 .OR. signout4 /= 0 ) 
 !        CALL OBJECTIVE(fune1, x0 + ee_i + ee_j, empirical, invSIG, size_data, equilibrium1,signout1)
  CALL OBJECTIVE(fune1,  x0 + ee_i + ee_j, empirical, invSIG)  
  CALL OBJECTIVE(fune2,  x0 - ee_i + ee_j, empirical, invSIG)  
  CALL OBJECTIVE(fune3,  x0 + ee_i - ee_j, empirical, invSIG) 
  CALL OBJECTIVE(fune4,  x0 - ee_i - ee_j, empirical, invSIG)
  !       CALL OBJECTIVE(fune2, x0 - ee_i + ee_j, empirical, invSIG, size_data, equilibrium2,signout2)
  !       CALL OBJECTIVE(fune3, x0 + ee_i - ee_j, empirical, invSIG, size_data, equilibrium3,signout3)
  !       CALL OBJECTIVE(fune4, x0 - ee_i - ee_j, empirical, invSIG, size_data, equilibrium4,signout4)
!          IF ( equilibrium1 /= 1 .OR. equilibrium2 /= 1 .OR. equilibrium3 /= 1 .OR. equilibrium4 /= 1 .OR. signout1 /= 0 .OR. signout2 /= 0 .OR. signout3 /= 0 .OR. signout4 /= 0 ) THEN
 !             ee_i = ee_i/2.0D+00
  !            ee_j = ee_j/2.0D+00
   !       END IF
!      END DO
      hessian( i,j ) = (fune1(1,1) - fune2(1,1) - fune3(1,1) + fune4(1,1))  / (4.0D+00 * dh(i) * dh(j))
      IF ( i /= j ) THEN
            hessian(j,i) = hessian(i,j)
      END IF
   END DO
END DO
END SUBROUTINE HESS

END MODULE MODEL_SOL
