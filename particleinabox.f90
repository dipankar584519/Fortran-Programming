!---------------------------------------------------------!
!-y''[x]+pV[x]y[x]=Lambda y[x]; V[x]=0; Lambda=Pi^2 n^2
!p=2m/hbar^2; Lambda=p E; y[0]=y[L]=0
!-q(y(i+1)+y(i-1))+(2q+pV(x))y(i)=Lambda y(i)
!---------------------------------------------------------!
        PROGRAM HEATEQN
        implicit NONE
        REAL::x0,xf,dx,p,q,Alpha,Beta,Norm2,v,ea1,stopping,GUESS
        INTEGER::i,j,k,h,n,m
        REAL,ALLOCATABLE,DIMENSION(:)::x,y1,p1,q1,y3,w1,x1
        REAL,ALLOCATABLE,DIMENSION(:,:)::A,Mm,A_inv
        OPEN(UNIT=19,FILE="heat.txt")
        ! REAL::
        !INTEGER::i,j,k,m,n,h
        !REAL,ALLOCATABLE,DIMENSION(:)::x,p1,y1,q1,y3,w1
        !REAL,ALLOCATABLE,DIMENSION(:,:)::Mm,A,A_inv
        !OPEN(UNIT=11,FILE="eigenvalue.inp")
        !WRITE(*,*)"GIVE THE ORDER OF THE MATRIX:"
        !READ(*,*)n     !oredr of the matrix
        stopping=0.01;x0=0; xf=1; n=10
        ALLOCATE(x(n),A(n-1,n-1),A_inv(n-1,n-1),p1(100),x1(n-1),y1(n-1),q1(100),y3(n-1),w1(100),Mm(n-1,n-1))
        !n is # for which x is max
        x(0)=x0; x(n)=xf
        dx=(xf-x0)/(n)
        DO i=1,n-1
        x(i)=x(i-1)+dx
        !WRITE(*,*)x(i),dx
        ENDDO
        
!setting LHS matrix
        DO i=1,n-1
            DO j=1,n-1
                IF(j==i)THEN
                Mm(i,i)=2/dx**2
                ELSEIF((j==(i-1).OR.(j==(i+1))))THEN
                Mm(i,j)=-1/dx**2
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
        v=p1(k)
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
        
        WRITE(19,*)"Eigenvalue in Loop N0=",k
        WRITE(19,'(f10.4)')v
        WRITE(19,*)"        x               y           "
        DO h=0,n
        WRITE(19,'(2(f10.4))')x(h),y1(h)/Norm2
        ENDDO      
        
        
        CLOSE(19)
        END PROGRAM HEATEQN
!-------------SUBROUTINE FOR INVERSE MATRIX-----------!
       
        SUBROUTINE Gauss_Jordon_Inversion(A,n,A_inv)
        IMPLICIT NONE
        INTEGER,INTENT(IN)::n
        INTEGER::i,j,k,it,n2,MaxRow
        REAL,INTENT(IN),DIMENSION(:,:)::A(n,n)
        REAL,INTENT(OUT),DIMENSION(:,:)::A_inv(n,n)
        REAL,DIMENSION(:,:)::IM(n,n)
        REAL,DIMENSION(:,:)::Aug(n,2*n),t(n,2*n)
        REAL::MaxEl,c
        
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
    REAL FUNCTION V(x)
        REAL,INTENT(IN)::x
        V=-4/x**2
        END FUNCTION V
