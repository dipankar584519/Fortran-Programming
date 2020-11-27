        PROGRAM HEATEQN
        implicit NONE
        REAL::x0,xf,T0,dx,h,p,q
        INTEGER::i,j,k,n 
        REAL,ALLOCATABLE,DIMENSION(:)::x,y,T
        REAL,ALLOCATABLE,DIMENSION(:,:)::A
        OPEN(UNIT=19,FILE="heat.txt")
        !WRITE(*,*)"n, xi,xf, y0"
        !READ(*,*)n,xi,xf,y0
        T0=20; x0=0; xf=10; n=9;h=0.01
        ALLOCATE(x(n),y(n),A(n,n),T(n))
        !n is #values we want (excluding 0 and n+1)
        x(0)=x0; x(n+1)=xf; T(0)=40; T(n+1)=200
        dx=(xf-x0)/(n+1)
        q=-(2+h*dx**2)
        p=-h*dx*dx*T0
        DO i=1,n
            DO j=1,n
                IF(j==i)THEN
                A(j,j)=q
                ELSEIF((j==(i+1)).OR.(j==(i-1)))THEN
                A(i,j)=1
                ELSE 
                A(i,j)=0
                ENDIF
            ENDDO
        ENDDO
        WRITE(*,*)"LHS matrix",q,p
        DO i=1,n
        WRITE(*,*)(A(i,j),j=1,n)
        ENDDO
        DO i=1,n
        x(i)=x(i-1)+dx
        y(i)=p
        ENDDO
        y(1)=p-T(0)
        y(n)=p-T(n+1)
        WRITE(*,*)"RHS"
        DO i=1,n
        WRITE(*,*)y(i)
        ENDDO
        WRITE(*,*)"     x(i)            T(i)"
        CALL MatrixGaussJordon(A,y,T,n)
        DO i=0,n+1
        WRITE(*,*)x(i),T(i)
        WRITE(19,*)x(i),T(i)
        ENDDO
        CLOSE(19)
        END PROGRAM HEATEQN
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

