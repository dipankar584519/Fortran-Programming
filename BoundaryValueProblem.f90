!---------------------------------------------------------!
!SOLVING: y''[x]+y[x]+xy=0,y[-1]=1,y'[1]=2,h=0.5
!---------------------------------------------------------!
        PROGRAM HEATEQN
        implicit NONE
        REAL::x0,xf,T0,dx,h,p,q,Alpha,Beta
        INTEGER::i,j,k,n 
        REAL,ALLOCATABLE,DIMENSION(:)::x,y,T
        REAL,ALLOCATABLE,DIMENSION(:,:)::A
        OPEN(UNIT=19,FILE="heat.txt")
        !WRITE(*,*)"n"
        !READ(*,*)n
        !Alpha=1;Beta=64;
        x0=-1; xf=1.5; n=5
        ALLOCATE(x(n),y(n-1),A(n-1,n-1),T(n))
        !n is # for which x is max
        x(0)=x0; x(n)=xf;T(0)=1
        dx=(xf-x0)/(n)
        DO i=1,n-1
        x(i)=x(i-1)+dx
        WRITE(*,*)x(i),dx
        ENDDO
        
        
!setting RHS coloumn vector
        y(1)=-(2-dx)*1
        y(n-1)=-(2+dx)*4*dx
        DO i=2,n-2
        y(i)=0
        ENDDO
!setting LHS matrix
        A(1,1)=-4+2*dx*dx*x(1); A(1,2)=2+dx      
        A(n-1,n-2)=4; A(n-1,n-1)=-4+2*dx*dx*x(n-1)
        DO i=2,n-2
            DO j=1,n-1
                IF(j==i)THEN
                 A(i,i)=-4+2*dx*dx*x(i)
                ELSEIF(j==(i-1))THEN
                A(i,j)=2-dx
                ELSEIF(j==(i+1))THEN
                A(i,j)=2+dx
                ELSE 
                A(i,j)=0
                ENDIF
            ENDDO    
        ENDDO
!writing lHS matrix and RHS coloumn vector
        WRITE(*,*)"LHS matrix"
        DO i=1,n-1
        WRITE(*,*)(A(i,j),j=1,n-1)
        ENDDO
        WRITE(*,*)"RHS"
        DO i=1,n-1
        WRITE(*,*)y(i)
        ENDDO
!solution of the linear equations
       WRITE(*,*)"     i            x(i)            T(i)"
        CALL MatrixGaussJordon(A,y,T,n-1)
        T(n)=4*dx+T(n-2)
        
        DO i=0,n
        WRITE(*,*)i,x(i),T(i)
        WRITE(19,*)i,x(i),T(i)
        ENDDO
       
        CLOSE(19)
        END PROGRAM HEATEQN
        
!----------------Q(x) FUNCTION-------------------------!
    REAL FUNCTION Q(x)
        REAL,INTENT(IN)::x
        Q=-4/x**2
        END FUNCTION Q

!-------SUBROUTINE OF SOLVING LINEAR EQUATIONS---------!        
        SUBROUTINE MatrixGaussJordon(B,v,x,m)
        IMPLICIT NONE
        INTEGER::i,j,k,m,n,MaxRow
        DOUBLE PRECISION::c,store,MaxEl
        REAL,DIMENSION(:,:)::B(m,m),t(m,m)
        REAL,DIMENSION(:)::v(m),x(m),s(m)
        DO k=1,m
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
