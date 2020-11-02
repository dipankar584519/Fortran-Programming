!-------------------------------------------------------!
!Lab PH707 (LG1)
!Title: Polynomial Fitting
!Date: 12/10/2020
!Name: Dipankar Barman
!Roll No: 206121013
!Email: dipankar1998@iitg.ac.in
!-------------------------------------------------------!
        PROGRAM Polynomial_Alternative_Fitting
        IMPLICIT NONE
        REAL::h,avg,sy
        INTEGER::i,j,k,m,n
        REAL,ALLOCATABLE,DIMENSION(:)::x,y,v,c
        REAL,ALLOCATABLE,DIMENSION(:,:)::Z,Z_inv
        OPEN(UNIT=11,FILE="PolynomialFit.inp")
        
        WRITE(*,*)"GIVE THE NUMBER OF COEFFICIENTS:"
        READ(*,*)m     !oredr of the polynomial
        WRITE(*,*)"GIVE THE NUMBER OF POINTS:"
        READ(*,*)n     !number of points
        ALLOCATE(x(n),y(n),c(m),Z(m,m),v(m),Z_inv(m,m))
        WRITE(*,*)"    x-value         y-value"
        DO i=1,n
            READ(11,*)x(i),y(i)
            WRITE(*,*)x(i),y(i)
        ENDDO
        WRITE(*,*)"Z is:"
        DO i=1,m
           DO j=1,m
              Z(i,j)=avg(x,y,i+j-2,0,n)
           ENDDO
           WRITE(*,*)(Z(i,j),j=1,m)
        ENDDO
        WRITE(*,*)"b is:"
        DO i=1,m
        v(i)=avg(x,y,i-1,1,n)
        WRITE(*,*)v(i)
        ENDDO
        WRITE(*,*)"SOLUTIONS ARE:" 
        CALL MatrixGaussJordon(Z,v,c,m)
        DO i=1,m
        j=i-1
        WRITE(*,*)"a[",j,"]=",c(i)
        ENDDO
        
        CALL least_square_fit(c,x,y,m,n,h,sy)
        WRITE(*,*)"Correlation Coefficient:C_r=",h
        WRITE(*,*)"S_{y/x}=",sy
     
        CLOSE(11)
        END PROGRAM Polynomial_Alternative_Fitting
!--------------FUNCTION FOR AVERAGE-------------------!
    REAL FUNCTION avg(x,y,a,b,n)
        REAL,INTENT(IN),DIMENSION(:)::x(n),y(n)
        INTEGER,INTENT(IN)::a,b,n 
        INTEGER::i 
        REAL::h
        h=0
        DO i=1,n
        h=h+x(i)**a*y(i)**b
        ENDDO
        avg=h/n
       END FUNCTION avg
!-------SUBROUTINE OF SOLVING LINEAR EQUATIONS---------!        
        SUBROUTINE MatrixGaussJordon(B,v,x,m)
        IMPLICIT NONE
        INTEGER::i,j,k,m,n,MaxRow
        DOUBLE PRECISION::c,store,MaxEl
        REAL,DIMENSION(:,:)::B(m,m),t(m,m)
        REAL,DIMENSION(:)::v(m),x(m),s(m)
        DO k=1,m
        !performing the row exchange operation 
          MaxEl=abs(B(k,k))
          MaxRow=k
         DO i=k+1,m
           IF(abs(B(i,k))>MaxEl) THEN
             MaxEL=abs(B(i,k))
             MaxRow=i
           END IF
          ENDDO
         IF(MaxRow/=k)THEN 
          DO i=k,m
           t(k,i)=B(maxRow,i)
           B(maxRow,i) = B(k,i)
           B(k,i) =t(k,i)
          ENDDO     
           s(k)=v(maxRow)
           v(maxRow)=v(k)
           v(k)=s(k)
         ENDIF 
         !performing elimination process
          DO i=1,m
            IF(i/=k)THEN
              c=B(i,k)/B(k,k)
              DO j=1,m
                B(i,j)=B(i,j)-c*B(k,j)
              ENDDO
              v(i)=v(i)-c*v(k)
            ENDIF
          ENDDO
        ENDDO
        !solutions of the matrix equations
         DO i=1,m
           x(i)=v(i)/B(i,i)
         ENDDO
        END SUBROUTINE MatrixGaussJordon
         
!-------------SUBROUTINE CORRELATION FUNCTION-----------!   

 SUBROUTINE least_square_fit(c,x,y,m,n,h,sy)
        IMPLICIT NONE
        INTEGER::i,j,k,m,n
        !DOUBLE PRECISION::
        REAL::d,t,e,h,sy
        !REAL,DIMENSION(:,:)::B(m,m),t(m,m)
        REAL,DIMENSION(:)::c(m),x(n),y(n),s(n),l(n),q(n)
                        
        d=SUM(y)/n
        l=0
        DO i=1,n
        DO j=1,m
        l(i)=l(i)+c(j)*x(i)**(j-1)
        ENDDO
        s(i)=(y(i)-l(i))**2
        q(i)=(y(i)-d)**2
        ENDDO
        t=SUM(s)
        e=SUM(q)
        h=SQRT((e-t)/e)
        sy=SQRT(t/(n-m))
        END SUBROUTINE least_square_fit
