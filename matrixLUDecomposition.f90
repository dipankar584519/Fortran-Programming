        PROGRAM MatrixLU !WITHOUT PIVOTING, MAY BLOW UP FOR 0 PIVOT
        IMPLICIT NONE
        INTEGER::i,j,k,m,n,MaxRow
        REAL::pivot,c,store,MaxEl
        REAL,ALLOCATABLE,DIMENSION(:)::B,x,y
        REAL,ALLOCATABLE,DIMENSION(:,:)::A,U,L,t !MATRIX OF UNKNOWN DIMENSION
        REAL::s
        OPEN(UNIT=11,FILE="MatrixInput.inp")
        OPEN(UNIT=19,FILE="matrix1.out")
        READ(11,*)n

        ALLOCATE(A(n,n+1),t(n,n+1),U(n,n),L(n,n),B(n),x(n),y(n))
        DO i=1,n
          READ(11,*)(A(i,j), j=1,n)
        ENDDO
        DO i=1,n
          READ(11,*)B(i)
          A(i,n+1)=B(i)
        ENDDO
        WRITE(19,*)"The Augmented matrix is:" 
        DO i=1,n
          WRITE(19,'(4(f6.2))')(A(i,j), j=1,n+1)
        ENDDO
        
    DO k=1,n
      DO j=1,n
        U(k,j)=0
        L(k,j)=0
      ENDDO
    ENDDO    

    DO k=1,n           
        DO i=k,n
          s=0.0
          DO j=1,k     
            s=s+L(k,j)*U(j,i)
          ENDDO
          U(k,i)=A(k,i)-s
        ENDDO
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

        WRITE(19,*)"Upper Triangular Matrix U is:"
        DO i=1,n
           WRITE(19,'(4(f10.5))')(U(i,j),j=1,n)
        ENDDO
        WRITE(19,*)"Lower Triangular Matrix L is:"
        DO i=1,n
           WRITE(19,'(4(f10.5))')(L(i,j),j=1,n)
        ENDDO
      t=MATMUL(L,U)

      WRITE(19,*)"Multiplication Of L and U is:"
      DO i=1,n
       WRITE(19,*)(t(i,j), j=1,n)     
      ENDDO
        CLOSE(11)
        CLOSE(19)
        END PROGRAM MatrixLU

