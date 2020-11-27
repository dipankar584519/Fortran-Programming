!---------------------------------------------------------!
!-y''[x]+pV[x]y[x]=Lambda y[x]; V[x]=x^2/2; Lambda=Pi^2 n^2
!p=2m/hbar^2; Lambda=p E; y[\pm\infty]=0
!-q(y(i+1)+y(i-1))+(2q+pV(x))y(i)=Lambda y(i)
!---------------------------------------------------------!
        PROGRAM HEATEQN
        implicit NONE
        DOUBLE PRECISION::x0,xf,dx,p,q,Alpha,Beta,Norm2,V,ea1,stopping,GUESS
        INTEGER::i,j,k,h,n,m
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::x,y1,p1,q1,y3,w1,x1
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::A,Mm,A_inv
        OPEN(UNIT=19,FILE="heat.txt")
        stopping=0.001;x0=-2.0; xf=2.0; n=20
        ALLOCATE(x(n),A(n-1,n-1),A_inv(n-1,n-1),p1(100),x1(n-1),y1(n-1),q1(100),y3(n-1),w1(100),Mm(n-1,n-1))
        x(0)=x0; x(n)=xf
        dx=(xf-x0)/(n)
        DO i=1,n-1
        x(i)=x(i-1)+dx
        ENDDO
        
!setting LHS matrix
        DO i=1,n-1
            DO j=1,n-1
                IF(j==i)THEN
                Mm(i,i)=1/dx**2+1*V(x(i))
                ELSEIF((j==(i-1).OR.(j==(i+1))))THEN
                Mm(i,j)=-0.5/dx**2
                ELSE 
                Mm(i,j)=0
                ENDIF
            ENDDO    
        ENDDO
!writing lHS matrix and RHS coloumn vector
       ! WRITE(*,*)"LHS matrix"
       ! DO i=1,n-1
       ! WRITE(*,*)(Mm(i,j),j=1,n-1)
       ! ENDDO
        WRITE(*,*)"ENTER THE GUESS VALUE:"
        READ(*,*)GUESS

       ! WRITE(*,*)"A=M-(Identity*GUESS) is:"
        DO i=1,n-1
            DO j=1,n-1
                IF(j==i)THEN
                    A(i,i)=Mm(i,i)-GUESS
                ELSE 
                    A(i,j)=Mm(i,j)
                ENDIF
            ENDDO
        ENDDO
       ! DO i=1,n-1
       !     WRITE(*,*)(A(i,j),j=1,n-1)
       ! ENDDO
       ! WRITE(*,*)"A_inv is:"
         CALL Gauss_Jordon_Inversion(A,n-1,A_inv)       
        
       ! DO i=1,n-1
       !     WRITE(*,*)(A_inv(i,j),j=1,n-1)
       ! ENDDO
       ! WRITE(*,*)"Initial guess:"
        DO i=1,n-1
            x1(i)=1
       !     WRITE(*,*)x1(i)
        ENDDO
 
        y1=matmul(A_inv,x1)/maxval(matmul(A_inv,x1),n-1)
      !  WRITE(*,*)"without Coefficient A_inv.X is:"
      !  DO i=1,n-1
      !   WRITE(*,*)y1(i)
      !  ENDDO
        
    Loop2:DO k=2,100
        y1=matmul(A_inv,y1)/maxval(matmul(A_inv,y1),n-1)
        y3=matmul(A_inv,y1)
        q1(k)=0
        DO h=1,n-1
        q1(k)=q1(k)+y3(h)*y1(h)
        ENDDO
        w1(k)=0
        DO h=1,n-1
        w1(k)=w1(k)+y1(h)*y1(h)
        ENDDO
        p1(k)=GUESS+w1(k)/q1(k)
       ! v=p1(k)
        ea1=abs((p1(k)-p1(k-1))/p1(k))*100
     !   WRITE(*,*)"Eigenvalue=",v,"Loop N0=",k
      !  WRITE(*,*)"Del(Eigenvalue)=",ea1
         Norm2=sqrt(sum(y1**2))
     !    WRITE(*,*)"Normalization is:",Norm2
      !   WRITE(*,*)"Normalized Eigenvector is:"
     !    DO h=0,n
      !  WRITE(*,*)x(h),y1(h)/Norm2
      !  ENDDO      
        !WRITE(*,*)0.000000000
        IF(ea1<stopping)EXIT Loop2
        ENDDO Loop2
        
        !WRITE(19,*)"Eigenvalue in Loop N0=",k
        !WRITE(19,'(f10.4)')p1(k)
        !WRITE(19,*)"        x               y           "
        DO h=0,n
        WRITE(19,'(2(f10.4))')x(h),y1(h)/Norm2
        ENDDO      
        WRITE(*,*)"Eigenvalue in Loop N0=",k
        WRITE(*,'(f10.4)')p1(k)
        WRITE(*,*)"        x               y           "
        DO h=0,n
        WRITE(*,'(2(f10.4))')x(h),y1(h)/Norm2
        ENDDO      
        
        CLOSE(19)
        END PROGRAM HEATEQN
!-------------SUBROUTINE FOR INVERSE MATRIX-----------!
       
        SUBROUTINE Gauss_Jordon_Inversion(A,n,A_inv)
        IMPLICIT NONE
        INTEGER,INTENT(IN)::n
        INTEGER::i,j,k,it,n2,MaxRow
        DOUBLE PRECISION,INTENT(IN),DIMENSION(:,:)::A(n,n)
        DOUBLE PRECISION,INTENT(OUT),DIMENSION(:,:)::A_inv(n,n)
        DOUBLE PRECISION,DIMENSION(:,:)::IM(n,n)
        DOUBLE PRECISION,DIMENSION(:,:)::Aug(n,2*n),t(n,2*n)
        DOUBLE PRECISION::MaxEl,c
        
        n2=2*n
!Generating Identity Matrix
        DO i=1,n
          DO j=1,n
            IF(i/=j)THEN
            IM(i,j)=0
            ELSE 
            IM(i,j)=1
            ENDIF
            k=n+j
            Aug(i,j)=A(i,j)
            Aug(i,k)=IM(i,j)
          ENDDO
        ENDDO
        
        DO k=1,n
!performing the row exchange operation 
          MaxEl=abs(Aug(k,k))
          MaxRow=k
         DO i=k+1,n
           IF(abs(Aug(i,k))>MaxEl) THEN
             MaxEL=abs(Aug(i,k))
             MaxRow=i
           END IF
          ENDDO
         IF(MaxRow/=k)THEN 
          DO i=k,n2
           t(k,i)=Aug(maxRow,i)
           Aug(maxRow,i) = Aug(k,i)
           Aug(k,i) =t(k,i)
          ENDDO     
        ENDIF 
!performing elimination process
          DO i=1,n
            IF(i/=k)THEN
              c=Aug(i,k)/Aug(k,k)
              DO j=1,n2
                Aug(i,j)=Aug(i,j)-c*Aug(k,j)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
!Inverse matrix is the right half of the reduced merged matrix!
         DO i=1,n
          DO j=1,n
            k=n+j
            A_inv(i,j)=Aug(i,k)/Aug(i,i)
          ENDDO
        ENDDO

        END SUBROUTINE Gauss_Jordon_Inversion
        
!----------------V(x) FUNCTION-------------------------!
    DOUBLE PRECISION FUNCTION V(x)
        DOUBLE PRECISION,INTENT(IN)::x
        DOUBLE PRECISION::b0,b1,b2,b3,b4,x0,x1,x2,x3,x4
        b0=0; b1=-5; b2=5; b3=-3; b4=1.3333
        x0=-2; x1=-1; x2=0; x3=1; x4=2
        V=b0+b1*(x-x0)+b2*(x-x0)*(x-x1)+b3*(x-x0)*(x-x1)*(x-x2)+b4*(x-x0)*(x-x1)*(x-x2)*(x-x3)
        END FUNCTION V
