MODULE MODEL_SOL

USE MPARAMS
USE MODEL
IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
SUBROUTINE  OBJECTIVE(obj, x, empirical, invSIG_emp, size_emp, equilibrium,signout)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: size_emp
INTEGER              :: equilibrium, signout
REAL(8), INTENT(IN)  :: empirical(size_emp,1), invSIG_emp(size_emp,size_emp), x(num_thetta)
REAL(8), INTENT(OUT) :: obj(1,1)
REAL(8)              :: emp_theory(size_emp,1), theory(size_emp,1)
INTEGER :: i, j

CALL MOMth(theory, size_emp, equilibrium, signout, x)

!write(*,*) 'momth',theory

IF ( signout == 1 .OR. equilibrium /= 1 ) THEN
   obj = - 10000.0D+00
   return
END IF

emp_theory = empirical - theory

IF (1==2) THEN
OPEN(unit = 22 ,file = "momth.txt", action = "write")
DO i =1, size_emp
    write(22,*) theory(i,1)
END DO
CLOSE(22)

OPEN(unit = 22 ,file = "emp_th.txt", action = "write")
DO i =1, size_emp
    write(22,*) emp_theory(i,1)
END DO
CLOSE(22)

OPEN(unit = 22 ,file = "emp.txt", action = "write")
DO i =1, size_emp
    write(22,*) empirical(i,1)
END DO
CLOSE(22)
END IF
obj = - 0.5D+0*dble(size_data)* MATMUL(MATMUL(TRANSPOSE(emp_theory),invSIG_emp),emp_theory)

END SUBROUTINE objective

SUBROUTINE MOMth(theory, size_theory, equilibrium, signout, x)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: size_theory
INTEGER, INTENT(OUT) :: equilibrium, signout
REAL(8), INTENT(IN)  :: x(num_thetta)
REAL(8), INTENT(OUT) :: theory(size_theory,1)
INTEGER              :: i, j
REAL(8)              :: init(num_shocks),  gx(sizey,sizex), hx(sizex,sizex)
REAL(8)              :: VARY(num_MOM,num_MOM, nlags+1), moments(num_MOM, num_MOM), sigx(sizex,sizex)

CALL gx_hx(gx,hx,x, signout, equilibrium)
!write(*,*) 'gx', gx
!write(*,*) 'hx',hx
IF ( 1 == 12) THEN
  OPEN(unit = 22 ,file = "gx.txt", action = "write")
  DO i =1, sizex
   DO j = 1,sizey
      write(22,*) gx(j,i)
   END DO
  END DO 
  CLOSE(22)
  OPEN(unit = 22 ,file = "hx.txt", action = "write")
  DO i =1, sizex
    DO j = 1,sizex
      write(22,*) hx(j,i)
    END DO
  END DO
  CLOSE(22)
END IF

init(1) = x(5) * 0.01D+0
init(2) = x(6) * 0.01D+0
init(3) = x(7)* 0.01D+0

CALL MOM( moments, sigx, gx, hx, init )

DO i = 1, nlags+1
   VARY(:,:,i) = moments
   !%Get E{x(t)*x(t-1)'}
   sigx = MATMUL( hx, sigx )!, sigx)
   !%Get E{y(t)*y(t+J)'}
   moments = ( MATMUL(MATMUL(gx,sigx),TRANSPOSE(gx)))
END DO

DO i = 1, nlags + 1
    j = i
    theory(i,1) = VARY(1,1,j)
END DO
DO i = nlags + 2 , 2*(nlags+1)
    j = i - nlags - 1 
    theory(i,1) = VARY(2,1,j)
END DO
DO i = 2*(nlags + 1)+1 , 3*(nlags+1)
    j = i - 2*(nlags +1) 
    theory(i,1) = VARY(3,1,j)
END DO
DO i = 3*(nlags+1)+1 , 4*(nlags + 1)-1
    j = i-3*(nlags+1)+1
    theory(i,1) = VARY(1,2,j)
END DO
DO i = 4*(nlags+1) , 5*(nlags+1) - 1
    j = i - 4*(nlags+1)+1
    theory(i,1) = VARY(2,2,j)
END DO
DO i = 5*(nlags + 1) , 6*(nlags+1)-1
    j = i - 5*(nlags +1)+1
    theory(i,1) = VARY(3,2,j)
END DO
DO i = 6*(nlags+1) , 7*(nlags + 1)-2
    j = i-6*(nlags+1)+2
    theory(i,1) = VARY(1,3,j)
END DO
DO i = 7*(nlags+1)-1 , 8*(nlags+1) - 3
    j = i - 7*(nlags+1)+3 
    theory(i,1) = VARY(2,3,j)
END DO
DO i = 8*(nlags + 1)-2 , 9*(nlags+1)-3
    j = i - 8*(nlags + 1) + 3
    theory(i,1) = VARY(3,3,j)
END DO

END SUBROUTINE MOMth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE solab( x, a, b, gx, lgx, hx, lhx, equilibrium )
IMPLICIT NONE
LOGICAL, PARAMETER     :: tf   = .TRUE.
REAL(8), PARAMETER     :: eps1 = 0.0D+00
INTEGER, INTENT(IN)    :: lgx, lhx
REAL(8), INTENT(INOUT) :: a(sizef,sizef),b(sizef,sizef)
REAL(8), INTENT(OUT)   :: gx(lgx,lhx), hx( lhx, lhx)
INTEGER, INTENT(OUT)   :: equilibrium
INTEGER                :: j, k, ierr,i
INTEGER                :: INFO
INTEGER, PARAMETER     :: LWORK = 10*sizef
REAL(8)                :: ALPHAR(sizef), ALPHAI(sizef), BETA(sizef), WORK(LWORK)
EXTERNAL DELZTG
COMPLEX(8)             :: V(sizef,sizef)
COMPLEX(8)             :: invV11( lhx, lhx)
COMPLEX(8)             :: Eigenvalue(sizef,sizef)
INTEGER                :: ipvt(lhx)
REAL(8)                :: rcond
COMPLEX(8)             :: work_zgeco(lhx), work_zgedi(lhx)
INTEGER,PARAMETER      :: JOB = 01
REAL(8)                :: VR(sizef,sizef), VL(1,sizef)

REAL(8), INTENT(IN)    :: x(num_thetta)

Eigenvalue = cmplx(0.0D+00,0.0D+00)
equilibrium = 3
gx          = 0.0D+00 
hx          = 0.0D+00

CALL DGGEV( 'N', 'V', sizef, a, sizef, b, sizef, ALPHAR, ALPHAI, BETA, VL, 1, VR, sizef, WORK, LWORK, INFO )
IF (INFO /= 0 ) THEN
    ! write(*,*) 'Warning: DGGEV did not calculate eigenvectors corectly in solab'
    gx = 0.0D+00
    hx = 0.0D+00
    equilibrium = 6
ELSE
    k = 0
    DO j = 1, sizef
      IF ( BETA(j) /= 0.0D+00 ) THEN
        IF ( ALPHAI(j) == 0.0D+00 .AND.   ABS(ALPHAR(j))< ABS(BETA(j)) ) THEN
          Eigenvalue( k+1 ,k+1) = cmplx( ALPHAR(j)/BETA(j), 0.0D+00 )
          V(:,k+1) = cmplx(VR(:, j),0.0D+00)
          k = k+1
         ELSEIF ( ALPHAI(j) > 0.0D+00 .AND. SQRT(ALPHAR(j)**2 + ALPHAI(j)**2) < ABS( BETA(j)) ) THEN
          Eigenvalue(k+1,k+1) = cmplx(ALPHAR(j)/BETA(j), ALPHAI(j)/BETA(j))
          Eigenvalue(k+2,k+2) = cmplx(ALPHAR(j)/BETA(j),-ALPHAI(j)/BETA(j))
          V(:,k+1)   = cmplx( VR(:,j),  VR(:,j+1))
          V(:,k+2)   = cmplx( VR(:,j), -VR(:,j+1))
          k = k + 2
        END IF
      END IF
    END DO
    IF ( k < lhx ) THEN
      gx = 0.0D+00
      hx = 0.0D+00
      equilibrium = 0
    ELSEIF ( k > lhx ) THEN
      gx = 0.0D+00
      hx = 0.0D+00
      equilibrium = 2
    ELSEIF ( k == lhx ) THEN
      invV11  = V(1:lhx,1:lhx)
      CALL zgeco ( invV11 , lhx, lhx, ipvt, rcond, work_zgeco )
      IF (rcond == 0.0D+00) THEN
           WRITE(*,*) 'Invertibility condition violated in solab'
           gx   = 0.0D+00
           hx   = 0.0D+00
           equilibrium = 4
       ELSE
           CALL zgedi( invV11, lhx, lhx, ipvt,  , work_zgedi, job )
           hx = real( MATMUL(MATMUL(V(1:lhx,1:lhx),Eigenvalue(1:k,1:k)),invV11),8)
           gx = real(MATMUL(V(lhx+1:lhx+lgx,1:lhx),invV11),8)
           equilibrium = 1
       END IF
    END IF
END IF
END SUBROUTINE solab

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

SUBROUTINE MOM( momentsy, momentsx , gx, hx, init)

IMPLICIT NONE
 
REAL(8), INTENT(OUT) :: momentsy(sizey,sizey), momentsx(sizex,sizex)
REAL(8), INTENT(IN)  :: init(num_shocks),  gx(sizey,sizex), hx(sizex,sizex)
REAL(8)              :: var_shock(sizex,sizex), err(sizex,num_shocks)
REAL(8)              :: hx_old(sizex,sizex), sig_old(sizex,sizex), sigx_old(sizex,sizex), diferenz, diff(sizex*sizex), sigx(sizex,sizex), sigy(sizey,sizey)
INTEGER:: i, max_element, j

err             = 0.0D+00
err( sth  , 1 ) = init(1)
err( seps , 2 ) = init(2)
err( sgam , 3 ) = init(3)

var_shock     = MATMUL(err,TRANSPOSE(err))

!Doubling algorithm
hx_old   = hx
sig_old  = var_shock
sigx_old = 0.0D+00
DO i = 1,sizex
    sigx_old(i,i) = 1.0D+00
END DO

diferenz = 0.1D+00
DO WHILE (diferenz > 1.0D-25) 
    sigx     = MATMUL(MATMUL( hx_old, sigx_old ), TRANSPOSE(hx_old) ) + sig_old
    diferenz = MAXVAL(ABS(sigx - sigx_old))
    sig_old  = MATMUL( MATMUL( hx_old,sig_old ), TRANSPOSE(hx_old) ) + sig_old
    hx_old   = MATMUL(hx_old,hx_old)
    sigx_old = sigx
END DO

!%Get E{x(t)*x(t+J)'}
!sigx = MATMUL( MATMUL( hx**(-min(0,lag)), sigx), (TRANSPOSE(hx))**(max(0,lag)))
!%Get E{y(t)*y(t+J)'}
momentsy = ( MATMUL(MATMUL(gx,sigx),TRANSPOSE(gx)))
momentsx = (sigx)

END SUBROUTINE MOM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HESSIAN_SYM(D, x, grdh, empirical, invSIG, size_emp, output)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: size_emp
INTEGER, INTENT(OUT) :: output
REAL(8), INTENT(IN)  :: x(num_thetta)
REAL(8), INTENT(IN)  :: empirical(size_emp), invSIG(size_emp),  grdh(num_thetta)
REAL(8), INTENT(OUT) :: D(num_thetta,num_thetta)
COMPLEX(16)          :: Eigenval(num_thetta,num_thetta), Eigenvec(num_thetta,num_thetta)
INTEGER              :: INFO, i , j
REAL(8)              :: WR(num_thetta), WI(num_thetta), VL(1,num_thetta), VR(num_thetta,num_thetta), WORK(4*num_thetta)

output = 1
CALL HESS(D, x, grdh, empirical, invSIG, size_emp)
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
SUBROUTINE HESS(hessian, x0, grdh, empirical, invSIG,  size_emp)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: size_emp
REAL(8), INTENT(OUT) :: hessian(num_thetta,num_thetta)
REAL(8), INTENT(IN)  :: x0(num_thetta), empirical(size_emp), invSIG(size_emp), grdh(num_thetta)
REAL(8) ::  fune1(1,1), fune2(1,1), fune3(1,1), fune4(1,1), eye(num_thetta,num_thetta), ones(num_thetta), stps, dh(num_thetta), ax0(num_thetta), dax0(num_thetta), hlp(num_thetta), xdh(num_thetta), ee(num_thetta,num_thetta), ee_i(num_thetta), ee_j(num_thetta)
REAL(8) :: gx(sizey,sizex), hx(sizex,sizex)
INTEGER :: i, j, k, equilibrium1, equilibrium2, equilibrium3, equilibrium4, signout1, signout2, signout3, signout4

ones = 1.0D+00

stps = EPSILON(1.0D+00)**(1.0/3.0)
! eps: floating point relative accuracy or machine precision: 2.22e-16
! stps: step size recommended by Dennis and Schnabel: 6.006e-6

!CALL OBJECTIVE( fune1, x0, empirical, invSIG, size_data, equilibrium1, signout1 )
CALL gx_hx(gx,hx,x0, signout1, equilibrium1)

IF (equilibrium1 /=1 .OR. signout1 /=0 ) THEN
   WRITE(*,*) 'Error in HESS : objective can not be calculated'
   return
END IF 

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

      equilibrium1 = 3
      equilibrium2 = 3
      equilibrium3 = 3
      equilibrium4 = 3
      signout1     = 2
      signout2     = 2
      signout3     = 2
      signout4     = 2

      DO WHILE ( equilibrium1 /= 1 .OR. equilibrium2 /= 1 .OR. equilibrium3 /= 1 .OR. equilibrium4 /= 1 .OR. signout1 /= 0 .OR. signout2 /= 0 .OR. signout3 /= 0 .OR. signout4 /= 0 ) 
         CALL OBJECTIVE(fune1, x0 + ee_i + ee_j, empirical, invSIG, size_emp, equilibrium1,signout1)
         CALL OBJECTIVE(fune2, x0 - ee_i + ee_j, empirical, invSIG, size_emp, equilibrium2,signout2)
         CALL OBJECTIVE(fune3, x0 + ee_i - ee_j, empirical, invSIG, size_emp, equilibrium3,signout3)
         CALL OBJECTIVE(fune4, x0 - ee_i - ee_j, empirical, invSIG, size_emp, equilibrium4,signout4)
          IF ( equilibrium1 /= 1 .OR. equilibrium2 /= 1 .OR. equilibrium3 /= 1 .OR. equilibrium4 /= 1 .OR. signout1 /= 0 .OR. signout2 /= 0 .OR. signout3 /= 0 .OR. signout4 /= 0 ) THEN
              ee_i = ee_i/2.0D+00
              ee_j = ee_j/2.0D+00
          END IF
      END DO
      hessian( i,j ) = (fune1(1,1) - fune2(1,1) - fune3(1,1) + fune4(1,1))  / (4.0D+00 * dh(i) * dh(j))
      IF ( i /= j ) THEN
            hessian(j,i) = hessian(i,j)
      END IF
   END DO
END DO
END SUBROUTINE HESS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GRADIENT_MOM(grad, x0, size_mom,  grdh)

IMPLICIT NONE

INTEGER, INTENT(IN)  :: size_mom
REAL(8), INTENT(OUT) :: grad(size_mom, num_thetta)
REAL(8), INTENT(IN)  :: x0(num_thetta), grdh(num_thetta)
INTEGER :: equilibrium, signout
REAL(8) ::  fune1(size_mom), fune2(size_mom), eye(num_thetta,num_thetta), ones(num_thetta), stps, dh(num_thetta), ax0(num_thetta), dax0(num_thetta), hlp(num_thetta), xdh(num_thetta),dhm(num_thetta,num_thetta), ee(num_thetta,num_thetta), ee_i(num_thetta), gx(sizey,sizex), hx(sizex,sizex), init(num_shocks)
INTEGER :: i, j

eye = 0.0D+00
DO j = 1, num_thetta
    eye(j,j) = 1.0D+00
END DO
ones = 1.0D+00!

stps = EPSILON(1.0D+00)**(1.0/3.0)

CALL gx_hx(gx,hx,x0, signout, equilibrium) 
IF ( SUM(gx) == 0.0D+00) THEN
   WRITE(*,*) 'Error in GRADIENT : objective can not be calculated '
   return
END IF 

!! ** Computation of stepsize (dh)
IF ( sum(grdh) > 0.0D+00 ) THEN
    dh = grdh
ELSE
    ax0 = ABS(x0)
    DO j = 1,num_thetta
       hlp(j) = MAX( x0(j),1.0D-02 )
    END DO
    IF ( SUM(x0) == 0.0D+00 ) THEN
        dax0 = x0/ax0
        dh  = stps * hlp * x0/ABS(x0)
    ELSE
        dh  = stps * hlp
    END IF
END IF

DO j = 1,num_thetta
    dhm(:,j) = dh
END DO

ee = eye*dhm

DO i = 1, num_thetta
      ee_i = ee(:,i)
      fune1 = 0.0D+00
      fune2 = 0.0D+00
      DO WHILE ( SUM(fune1) == 0.0D+00 .AND. SUM(fune2) == 0.0D+00 ) 
         CALL  MOMth(fune1, size_mom, equilibrium, signout, x0 + ee_i )
         CALL  MOMth(fune2, size_mom, equilibrium, signout, x0 - ee_i )
           ee_i = ee_i/2.0D+00
      END DO
      grad(:,i) = (fune1 - fune2)  / 2.0D+0/dh(i)
END DO
END SUBROUTINE GRADIENT_MOM

!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gx_hx(gx_part, hx_part, x, signout, equilibrium)
IMPLICIT NONE

REAL(8), INTENT(IN)  :: x(num_thetta)
REAL(8), INTENT(OUT) :: gx_part(sizey,sizex), hx_part(sizex,sizex)
REAL(8)              :: nf(sizef),nfx(sizef*(sizex-sizetheta)),nfxp(sizef*(sizex-sizetheta)),nfy(sizef*sizey),nfyp(sizef*sizey), A_cu(sizef,sizef), A_cup(sizef,sizef)
INTEGER, INTENT(OUT) :: equilibrium, signout
INTEGER              :: j, i           
REAL(8) :: params(33), gx(sizey,sizex-sizetheta), hx(sizex-sizetheta, sizex - sizetheta)

CALL fxfyfxpfypf(x,  nf, nfx, nfxp, nfy, nfyp, signout )

IF (1 == 2) THEN
OPEN(unit = 22 ,file = "NFX.txt", action = "write")
DO i = 1, sizef*( sizex-sizetheta)
  write(22,*) nfx(i)
END DO
CLOSE(22)

OPEN(unit = 22 ,file = "NFXP.txt", action = "write")
DO i =1, sizef*(sizex-sizetheta)
  write(22,*) nfxp(i)
END DO
CLOSE(22)

OPEN(unit = 22 ,file = "NFY.txt", action = "write")
DO i =1, sizef*sizey
  write(22,*) nfy(i)
END DO
CLOSE(22)

OPEN(unit = 22 ,file = "NFYP.txt", action = "write")
DO i =1, sizef*sizey
  write(22,*) nfyp(i)
END DO
CLOSE(22)

END IF

IF (signout == 1) THEN
   equilibrium = 0
   gx = 0.0D+00
   hx = 0.0D+00
   return
END IF

A_cup(: ,         1 : sizex - sizetheta ) = -RESHAPE(nfxp,(/sizef,sizex-sizetheta/))
A_cup(: , sizex - sizetheta + 1 : sizef ) = -RESHAPE(nfyp,(/sizef,sizey/))
A_cu( : ,         1 : sizex - sizetheta ) =  RESHAPE(nfx, (/sizef,sizex-sizetheta/))
A_cu( : , sizex -sizetheta + 1 : sizef  ) =  RESHAPE(nfy, (/sizef,sizey/))

CALL solab(x, A_cu, A_cup, gx, sizey, hx, sizex-sizetheta,equilibrium )

IF (1 == 2) THEN
OPEN(unit = 22 ,file = "gx.txt", action = "write")
DO i =1, sizex-sizetheta
 DO j = 1,sizey
  write(22,*) gx(j,i)
 END DO
END DO
CLOSE(22)
OPEN(unit = 22 ,file = "hx.txt", action = "write")
DO i =1, sizex-sizetheta
 DO j = 1,sizex-sizetheta
  write(22,*) hx(j,i)
 END DO
END DO
CLOSE(22)
END IF

IF (sizetheta == 0) THEN
   gx_part = gx
   hx_part = hx
ELSE 
   CALL gx_hx_partstate(gx_part, hx_part, gx, hx, RESHAPE(nfxp,(/sizef,sizex-sizetheta/)), RESHAPE(nfy, (/sizef,sizey/) ), RESHAPE(nfyp, (/sizef,sizey/) ) )
END IF


END SUBROUTINE gx_hx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gx_hx_partstate(gX_part, hX_part, GX, HX, nfXP, nfY, nfYP)

REAL(8), INTENT(IN)  :: GX(sizey,sizex-sizetheta), HX(sizex-sizetheta,sizex-sizetheta)
REAL(8), INTENT(OUT) :: gX_part(sizey,sizex), hX_part(sizex,sizex)
REAL(8)              :: P(sizetheta,sizetheta), nf(sizef) , nfXP(sizef,sizex-sizetheta), nfY(sizef,sizey), nfYP(sizef,sizey)
INTEGER              :: j, i, sizeym, sizezm, sizethetam, sizethetam_1, nthetam, m, n, INFO
REAL(8)              :: a( sizef, sizex-sizetheta  ), inva(sizef-sizetheta - sizez,sizef-sizetheta-sizez ), hgtheta(sizef-sizetheta,sizetheta), hgtheta_ba1(sizef-sizetheta,sizetheta)

hX_part = 0.0D+0
gX_part = 0.0D+0

P = HX(sizex-2*sizetheta+1: sizex-sizetheta,sizex-2*sizetheta+1:sizex-sizetheta)

sizeym       = sizey - sum(nz)
sizezm       = sum(nz)
sizethetam   = sizetheta
sizethetam_1 = sizetheta - ntheta(nSUBPER)

a = MATMUL( nfYP, GX ) + nfXP

DO n = 0, nSUBPER-1
 m = nSUBPER - n
    gX_part( sizeym+1 : sizey, sizex-sizetheta + sizethetam_1 + 1 : sizex-sizetheta + sizethetam ) = MATMUL( GX( sizeym + 1: sizey, sizex-2*sizetheta + sizethetam_1 + 1 : sizex- 2*sizetheta + sizethetam), P( sizethetam_1 + 1 : sizethetam, sizethetam_1 + 1 : sizethetam ) )
    
inva = 0.0D+0
inva(1:sizef-sizetheta-sizezm,1:sizex-2*sizetheta) = a(sizezm+1:sizef-sizetheta,1:sizex-2*sizetheta)
inva(1:sizef-sizetheta-sizezm,sizex-2*sizetheta+1:sizex-2*sizetheta+sizeym) = nfY(sizezm+1 : sizef-sizetheta, 1:sizeym)

DO i = sizef-sizetheta-sizezm+1, sizef-sizetheta- sizez
    inva(i,i) = 1.0D+0
END DO

CALL MATINV( inva , sizef-sizetheta-sizez , INFO)

IF (INFO /= 0 ) THEN
write(*,*) 'problem occured in gx hx partstate'
    gX_part = 0.0D+0
    hX_part = 0.0D+0
    return
END IF
    
     hgtheta =  MATMUL( MATMUL( inva(1:sizex-2*sizetheta+sizeym,1:sizex-2*sizetheta+sizeym), nfY(sizezm+1:sizef-sizetheta,sizeym+1: sizey ) ) , GX( sizeym + 1:sizey, sizex-2*sizetheta + sizethetam_1 + 1 : sizex-2*sizetheta + sizethetam) )
     hgtheta_ba1 = - MATMUL( hgtheta, P( sizethetam_1 + 1 : sizethetam, sizethetam_1 + 1 : sizethetam ) )
     hX_part(1:sizex-2*sizetheta, sizex-sizetheta + sizethetam_1 + 1 : sizex - sizetheta + sizethetam) = hgtheta_ba1(1:sizex-2*sizetheta,:)
     gX_part(1:sizeym           , sizex-sizetheta + sizethetam_1 + 1 : sizex - sizetheta + sizethetam) = hgtheta_ba1(sizex-2*sizetheta+1:sizex-2*sizetheta+sizeym,:)

    hX_part( 1:sizex-2*sizetheta , sizex-2*sizetheta + 1 : sizex-2*sizetheta + sizethetam) = hgtheta(1:sizex-2*sizetheta,:) +  HX(1:sizex-2*sizetheta ,sizex - 2*sizetheta + sizethetam_1 + 1 : sizex - 2*sizetheta + sizethetam)
    gX_part( 1:sizeym, sizex-2*sizetheta + sizethetam_1 + 1 : sizex-2*sizetheta + sizethetam ) = hgtheta(sizex-2*sizetheta+1:sizex-2*sizetheta+sizeym,:) +  GX(1:sizeym,sizex - 2*sizetheta + sizethetam_1 + 1 : sizex - 2*sizetheta + sizethetam)

    ! update nfy and nfz and the size of y vector
    sizeym     = sizeym + nz(m)
    sizezm     = sizezm - nz(m)
    sizethetam = sizethetam - ntheta(m)
    IF ( m == 1) THEN
        sizethetam_1 = 0  
    ELSE 
        sizethetam_1 =  sizethetam - ntheta(m-1)
    END IF
END DO

GX_part( :     , 1 : sizex - 2*sizetheta) = GX(:         , 1 : sizex - 2*sizetheta)
hX_part(1:sizex - 2*sizetheta, 1 : sizex - 2*sizetheta) = HX( 1 : sizex- 2*sizetheta, 1 : sizex - 2* sizetheta)
hX_part(sizex-2*sizetheta+1       : sizex-sizetheta,sizex-2*sizetheta+1:sizex - sizetheta) =  P
DO i = 1, sizetheta
   hX_part( sizex-sizetheta + i ,  sizex-2*sizetheta +i) = 1.0D+0
END DO

END SUBROUTINE gx_hx_partstate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE covM(covdata, data, sizedata1, sizedata2)
IMPLICIT NONE

INTEGER, INTENT(IN)  :: sizedata1, sizedata2
REAL(8), INTENT(IN)  :: data(sizedata1,sizedata2)
REAL(8), INTENT(OUT) :: covdata(sizedata2,sizedata2)
INTEGER              :: i
REAL(8) :: meandata(1,sizedata2), data_mean(sizedata1,sizedata2)


DO i = 1, sizedata2
   meandata(1,i) = sum(data(:,i))/dble(sizedata1)
END DO

DO i = 1,sizedata1
   data_mean(i,:) = data(i,:) - meandata(1,:)
END DO
covdata = MATMUL(TRANSPOSE(data_mean), data_mean)/dble(sizedata1)

END SUBROUTINE covM
END MODULE MODEL_SOL
