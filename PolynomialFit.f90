!-------------------------------------------------------!
!Lab PH707 (LG1)
!Title: Polynomial Fitting
!Date: 13/10/2020
!Name: Dipankar Barman
!Roll No: 206121013
!Email: dipankar1998@iitg.ac.in
!-------------------------------------------------------!
!This program has been modified in Lab to show S, S_y explicitly.
!-------------------------------------------------------!
        PROGRAM Polynomial_Fitting
        IMPLICIT NONE
        REAL::f,e,t,h,sy,ss
        INTEGER::i,j,k,m,m1,n
        REAL,ALLOCATABLE,DIMENSION(:)::x,y,u,v,c,s
        REAL,ALLOCATABLE,DIMENSION(:,:)::B,B_inv,A,Z,ZT
        OPEN(UNIT=11,FILE="xy.txt")
        OPEN(UNIT=19,FILE="fit.out")
        WRITE(*,*)"GIVE THE DEGREE OF THE POLYNOIAL:"
        READ(*,*)m1     !oredr of the polynomial
        m=m1+1
        WRITE(*,*)"GIVE THE NUMBER OF POINTS:"
        READ(*,*)n     !number of points
        ALLOCATE(x(n),y(n),s(n),c(m),A(n,2),Z(n,m),ZT(m,n),u(m),v(m),B(m,m),B_inv(m,m))
        DO i=1,n
            READ(11,*)(A(i,j),j=1,2)
            x(i)=A(i,1)
            y(i)=A(i,2)
        ENDDO
        !WRITE(*,*)"Z is:"
        !WRITE(19,*)"Z is:"
        DO i=1,n
           DO j=1,m
              Z(i,j)=x(i)**(j-1)
           ENDDO
           !WRITE(19,*)(Z(i,j),j=1,m)
           !WRITE(*,*)(Z(i,j),j=1,m)
        ENDDO
        DO i=1,m
          DO j=1,n
           ZT(i,j)=Z(j,i)
          ENDDO
        ENDDO
        WRITE(19,*)"ZTZ is:"
        WRITE(*,*)"ZTZ is:"
        DO i=1,m
           DO j=1,m
              B(i,j)=0
              v(i)=0
              DO k=1,n
                 B(i,j)=B(i,j)+ZT(i,k)*Z(k,j)
                 v(i)=v(i)+ZT(i,k)*y(k)
              ENDDO
           ENDDO
           WRITE(19,*)(B(i,j),j=1,m)
           WRITE(*,*)(B(i,j),j=1,m)
        ENDDO
        WRITE(19,*)"ZTY is:"
        DO i=1,m
        WRITE(19,*)v(i)
        ENDDO
        WRITE(19,*)"SOLUTIONS ARE:" 
        CALL MatrixGaussJordon(B,v,c,m)
        DO i=1,m
        j=i-1
        WRITE(*,*)"a[",j,"]=",c(i)
        WRITE(19,*)"a[",j,"]=",c(i)
        ENDDO
        WRITE(19,*) 
        CALL least_square_fit(c,x,y,m,n,h,sy,ss)
        WRITE(*,*)"Correlation Coefficient:C_r=",h
        WRITE(19,*)"Correlation Coefficient:C_r=",h
        WRITE(*,*)"S=",sy
        WRITE(19,*)"S=",sy
        WRITE(*,*)"S_y=",ss
        WRITE(19,*)"S_y=",ss
        
        CLOSE(11)
        CLOSE(19)
        
        END PROGRAM Polynomial_Fitting

!------------SUBROUTINE OF SOLVING LINEAR EQUATIONS-------------!        
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
         
!-----------------SUBROUTINE CORRELATION FUNCTION----------------!   

 SUBROUTINE least_square_fit(c,x,y,m,n,h,t,e)
        IMPLICIT NONE
        INTEGER::i,j,k,m,n
        !DOUBLE PRECISION::
        REAL::d,t,e,h,sy,ss
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
        END SUBROUTINE least_square_fit
