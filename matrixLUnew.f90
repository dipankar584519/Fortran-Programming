 !------------------------------------------------------!
 !A=LU; Ly=LUx=B
 !WITHOUT PIVOTING, MAY BLOW UP FOR 0 PIVOT
 !------------------------------------------------------!
PROGRAM MatrixLU 
        IMPLICIT NONE
        INTEGER::i,j,k,m,n,MaxRow
        REAL::pivot,c,store,MaxEl
        REAL,ALLOCATABLE,DIMENSION(:)::B,x,y
        REAL,ALLOCATABLE,DIMENSION(:,:)::A,U,L,t !MATRIX OF UNKNOWN DIMENSION
        REAL::s
        OPEN(UNIT=11,FILE="MatrixInput.inp")
        READ(11,*)n

        ALLOCATE(A(n,n+1),t(n,n+1),U(n,n),L(n,n),B(n),x(n),y(n))
        DO i=1,n
          READ(11,*)(A(i,j), j=1,n)
        ENDDO
        DO i=1,n
          READ(11,*)B(i)
          A(i,n+1)=B(i)
        ENDDO
        WRITE(*,*)"The Augmented matrix is:" 
        DO i=1,n
          WRITE(*,'(4(f8.3))')(A(i,j), j=1,n+1)
        ENDDO
!LU Decomposition 
    DO k=1,n    
        !Getting the k th row of U matrix
        DO i=k,n
          s=0.0
          DO j=1,k     
            s=s+L(k,j)*U(j,i)
          ENDDO
          U(k,i)=A(k,i)-s
        ENDDO
        !Getting the k th coloumn of L matrix
        DO i=k,n
          IF(i==k)THEN
            L(k,k)=1
          ELSE
            s=0
            DO j=1,k
              s=s+L(i,j)*U(j,k)
            ENDDO
            L(i,k)=(A(i,k)-s)/U(k,k)
          ENDIF
        ENDDO
    ENDDO
!Forward Substitution to get y
    y(1)=B(1)
    DO i=2,n
      s=0
        DO j=1,i-1
        s=s+L(i,j)*B(j)
        ENDDO
      y(i)=B(i)-s 
    ENDDO
!Backward Substitution to get x
    x(n)=y(n)/U(n,n)
    DO i=n-1,1,-1
      s=0
      DO j=i+1,n
      s=s+u(i,j)*x(j)
      ENDDO
      x(i)=(y(i)-s)/u(i,i)
    ENDDO
    

        WRITE(*,*)"Upper Triangular Matrix U is:"
        DO i=1,n
           WRITE(*,'(4(f8.3))')(U(i,j),j=1,n)
        ENDDO
        WRITE(*,*)"Lower Triangular Matrix L is:"
        DO i=1,n
           WRITE(*,'(4(f8.3))')(L(i,j),j=1,n)
        ENDDO
      t=MATMUL(L,U)

      WRITE(*,*)"Multiplication Of L and U is:"
      DO i=1,n
       WRITE(*,'(4(f8.3))')(t(i,j), j=1,n) 
      ENDDO
      WRITE(*,*)"y is:"
      DO i=1,n
      WRITE(*,*)y(i)
      ENDDO
      WRITE(*,*)"x is:"
      DO i=1,n
      WRITE(*,*)x(i)
      ENDDO
   
        CLOSE(11)
        END PROGRAM MatrixLU
