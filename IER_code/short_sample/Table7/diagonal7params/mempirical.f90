MODULE EMPIRICAL

USE FUNCTIONS
USE MPARAMS
USE MODEL
USE MODEL_SOL

IMPLICIT NONE

CONTAINS

SUBROUTINE sample_mom( hatMOM, invSIG, V, nMOM, data_sample )

IMPLICIT NONE

INTEGER, INTENT(IN)  :: nMOM
REAL(8), INTENT(IN)  :: data_sample(size_data,num_MOM)
REAL(8), INTENT(OUT) :: hatMOM(nMOM,1), invSIG(nMOM,nMOM), V(nMOM,nMOM)
REAL(8)              :: ap(num_MOM)

REAL(8)  :: invSIG_dble(nMOM,nMOM)
INTEGER :: i, INFO, j, p, h

REAL(8) :: Ghat( nMOM,nMOM+num_MOM), fz_a(nMOM +num_MOM, size_data-nlags) , fz(nMOM+num_MOM , size_data-nlags) , ahat(nMOM+num_MOM), sigA(nMOM+num_MOM,nMOM+num_MOM), Ch(nMOM+num_MOM,nMOM+num_MOM),  eye(nMOM, nMOM)

eye = 0.0D+0
DO i = 1,nMOM
   eye(i,i) = 1.0D+0 
END DO

ap = 0.0D+0
DO i = 1, num_MOM
  DO j = 1, size_data 
      ap(i)    = ap(i) + data_sample(j,i)
  END DO
END DO
ap = ap/dble(size_data)
! 1, 2   (r,r)   (r-1, r)
DO i = 1, nlags + 1
    h = i-1
    fz( i , : ) = data_sample(nlags + 1:size_data, nr)*data_sample(nlags + 1 - h:size_data - h, nr)
    ahat(i) = sum( fz( i , : ))/dble(size_data-nlags)
    DO j = 1,size_data-nlags
       fz_a( i , j ) = fz( i , j ) - ahat(i)
    END DO
    hatMOM(i,1) = ahat(i) - ap(nr)**2
END DO

! 3, 4 (r , y)   (y-1, r)
DO i = nlags + 2 , 2*(nlags+1)
!    write(*,*) 'nr, noutput', nr, noutput
    h = i - nlags - 2
    fz( i , : ) = data_sample(nlags + 1:size_data, noutput)*data_sample(nlags + 1 - h:size_data - h, nr)

    ahat(i) = sum( fz( i , : ))/dble(size_data-nlags)
    DO j = 1,size_data-nlags
       fz_a( i , j ) = fz( i , j ) - ahat(i)
    END DO
    hatMOM(i,1) = ahat(i) - ap(nr)*ap(noutput)
END DO

! 5, 6 (pai,r )   (pai-1, r)
DO i = 2*(nlags + 1)+1 , 3*(nlags+1)
    h = i - 2*(nlags +1) - 1
    fz( i , : ) = data_sample(nlags + 1:size_data, npai)*data_sample(nlags + 1 - h:size_data - h, nr)
    ahat(i) = sum( fz( i , : ))/dble(size_data-nlags)
    DO j = 1,size_data-nlags
       fz_a( i , j ) = fz( i , j ) - ahat(i)
    END DO
    hatMOM(i,1) = ahat(i) - ap(nr)*ap(npai)
END DO

! 7    (r-1,y)
DO  i = 3*(nlags+1)+1 , 4*(nlags + 1)-1
    h = i-3*(nlags+1)
    fz( i , : ) = data_sample(nlags + 1:size_data, nr)*data_sample(nlags + 1 - h:size_data - h, noutput)
    ahat(i) = sum( fz( i , : ))/dble(size_data-nlags)
    DO j = 1,size_data-nlags
       fz_a( i , j ) = fz( i , j ) - ahat(i)
    END DO
    hatMOM(i,1) = ahat(i) - ap(nr)*ap(noutput)
END DO

! 8 , 9  (y, y)  (y, y-1)
DO i = 4*(nlags+1) , 5*(nlags+1) - 1
    h = i - 4*(nlags+1) 
    fz( i , : ) = data_sample(nlags + 1:size_data, noutput)*data_sample(nlags + 1 - h:size_data - h, noutput)
    ahat(i) = sum( fz( i , : ))/dble(size_data-nlags)
    fz_a( i , : ) = fz( i , : ) - ahat(i)
    DO j = 1,size_data-nlags
       fz_a( i , j ) = fz( i , j ) - ahat(i)
    END DO
    hatMOM(i,1) = ahat(i) - ap(noutput)**2
END DO

! 10 , 11 ( pai, y) (pai-1, y)
DO i = 5*(nlags + 1) , 6*(nlags+1)-1
    h = i - 5*(nlags +1)
    fz( i , : ) = data_sample(nlags + 1:size_data, npai)*data_sample(nlags + 1 - h:size_data - h, noutput)
    ahat(i) = sum( fz( i , : ))/dble(size_data-nlags)
    DO j = 1,size_data-nlags
       fz_a( i , j ) = fz( i , j ) - ahat(i)
    END DO
    hatMOM(i,1) = ahat(i) - ap(npai)*ap(noutput)
END DO

! 12 (r-1, pai) 
DO i = 6*(nlags+1) , 7*(nlags + 1)-2
    h = i-6*(nlags+1)+1
    fz( i , : ) = data_sample(nlags + 1:size_data, nr)*data_sample(nlags + 1 - h:size_data - h, npai)
    ahat(i) = sum( fz( i , : ))/dble(size_data-nlags)
    DO j = 1,size_data-nlags
       fz_a( i , j ) = fz( i , j ) - ahat(i)
    END DO
    hatMOM(i,1) = ahat(i) - ap(npai)*ap(nr)
END DO

! 13 (y-1, pai) 
DO  i = 7*(nlags+1)-1 , 8*(nlags+1) - 3
    h = i - 7*(nlags+1)+2
    fz( i , : ) = data_sample(nlags + 1:size_data, noutput)*data_sample(nlags + 1 - h:size_data - h, npai)
    ahat(i) = sum( fz( i , : ))/dble(size_data-nlags)
    DO j = 1,size_data-nlags
       fz_a( i , j ) = fz( i , j ) - ahat(i)
    END DO
    hatMOM(i,1) = ahat(i) - ap(npai)*ap(noutput)
END DO
! 14 15 (pai pai) (pai pai-1)
DO i = 8*(nlags + 1)-2 , 9*(nlags+1)-3
    h = i - 8*(nlags + 1) + 2
    fz( i , : ) = data_sample(nlags + 1:size_data, npai)*data_sample(nlags + 1 - h:size_data - h, npai)
    ahat(i) = sum( fz( i , : ))/dble(size_data-nlags)
    DO j = 1,size_data-nlags
       fz_a( i , j ) = fz( i , j ) - ahat(i)
    END DO
    hatMOM(i,1) = ahat(i) - ap(npai)**2
END DO

fz(  9*(H+1)-2,  : ) = data_sample(nlags+1:size_data,nr)
ahat(9*(H+1)-2     ) = sum( fz( 9*(H+1)-2 , : ))/dble(size_data-nlags)
fz_a( 9*(H+1)-2, : ) = fz( 9*(H+1)-2 , : ) - ahat(9*(H+1)-2)
fz( 9*(H+1)-1,   : ) = data_sample(H+1:size_data,noutput)
ahat(9*(H+1)-1) = sum( fz( 9*(H+1)-1 , : ))/dble(size_data-nlags)
fz_a( 9*(H+1)-1 , : )= fz( 9*(H+1)-1 , : ) - ahat(9*(H+1)-1)
fz( 9*(H+1)  , : )   = data_sample(H+1:size_data,npai)
ahat(9*(H+1)) = sum( fz( 9*(H+1), : ))/dble(size_data - nlags)
fz_a( 9*(H+1) , : ) = fz( 9*(H+1) , : ) - ahat(9*(H+1))

p = (size_data - nlags)**(1.0D+0/3.0D+0)
    sigA = MATMUL(fz_a, TRANSPOSE(fz_a))/dble(size_data-nlags)
    DO h = 1,p
        Ch =  MATMUL(fz_a(:,h+1:size_data-nlags), TRANSPOSE(fz_a(:,1:size_data-nlags-h)))/dble(size_data-nlags) 
        sigA = sigA + (1.0D+0-dble(h)/dble(p+1))*(Ch + TRANSPOSE(Ch))
    END DO
    !Ghat = zeros(num_MOM**2*(1+H)- num_MOM,num_MOM**2*(1+H))
    Ghat = 0.0D+0
!    Ghat(1:end,1:num_MOM**2*(1+H)- num_MOM) = eye(num_MOM**2*(1+H)- num_MOM)
DO i = 1,num_MOM**2*(1+H)- num_MOM
   Ghat(i,i) = 1.0D+0 
END DO 

Ghat(1          : H+1         , num_MOM**2*(1+H)- num_MOM + 1)  = -2*ap(nr)
Ghat(H+2         : 2*(H+1)     , num_MOM**2*(1+H)- num_MOM + 1) = -ap(noutput)
Ghat(2*(H+1) + 1 : 3*(H+1)     , num_MOM**2*(1+H)- num_MOM + 1) = -ap(npai)
Ghat(3*(H+1) + 1 : 4*(H+1) - 1 , num_MOM**2*(1+H)- num_MOM + 1) = -ap(noutput)
Ghat(6*(H+1)     : 7*(H+1) - 2 , num_MOM**2*(1+H)- num_MOM + 1) = -2*ap(npai)

Ghat(  (H+1) + 1 : 2*(H+1)     , num_MOM**2*(1+H)- num_MOM + 2) = -ap(nr)
Ghat(3*(H+1) + 1 : 4*(H+1) - 1 , num_MOM**2*(1+H)- num_MOM + 2) = -ap(nr)
Ghat(4*(H+1) + 1 : 5*(H+1) - 1 , num_MOM**2*(1+H)- num_MOM + 2) = -2*ap(noutput)
Ghat(5*(H+1) + 1 : 6*(H+1) - 1 , num_MOM**2*(1+H)- num_MOM + 2) = -ap(npai)
Ghat(7*(H+1) - 1 : 8*(H+1) - 3 , num_MOM**2*(1+H)- num_MOM + 2) = -ap(npai)

Ghat(2*(H+1) + 1 : 3*(H+1)     , num_MOM**2*(1+H)- num_MOM + 3) = -ap(nr)
Ghat(5*(H+1)     : 6*(H+1) - 1 , num_MOM**2*(1+H)- num_MOM + 3) = -ap(noutput)
Ghat(6*(H+1)     : 7*(H+1) - 2 , num_MOM**2*(1+H)- num_MOM + 3) = -ap(nr)
Ghat(7*(H+1) - 1 : 8*(H+1) - 3 , num_MOM**2*(1+H)- num_MOM + 3) = -ap(noutput)
Ghat(8*(H+1) - 2 : 9*(H+1) - 3 , num_MOM**2*(1+H)- num_MOM + 3) = -2*ap(npai)

V = MATMUL(Ghat, MATMUL(sigA, TRANSPOSE(Ghat)))
invSIG = V
CALL  MATINV(invSIG, nMOM , INFO)

END SUBROUTINE sample_mom
END MODULE EMPIRICAL
