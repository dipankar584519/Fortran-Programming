        PROGRAM matrixMultiplication
        IMPLICIT NONE
        INTEGER::i,j,k,m,n
        REAL,ALLOCATABLE,DIMENSION(:,:)::a,b,t!MATRIX OF UNKNOWN DIMENSION
        OPEN(UNIT=11,FILE="matrix.inp")
        OPEN(UNIT=19,FILE="matrix.out")
        WRITE(*,*)"Give dimensions of the matrix A,B:"
        READ(*,*)m,n

        ALLOCATE(a(m,n),b(n,m),t(m,m))

        DO i=1,m
          READ(11,*)(a(i,j), j=1,n)
        ENDDO

        DO i=1,n
          READ(11,*)(b(i,j),j=1,m)
        ENDDO
        
        DO i=1,m
         DO j=1,m
         t(i,j)=0
          DO k=1,n
          t(i,j)=t(i,j)+a(i,k)*b(k,j)
          ENDDO
         ENDDO
        ENDDO

        DO i=1,m
          WRITE(19,*)(a(i,j), j=1,n)
        ENDDO
        WRITE(19,*)
        DO i=1,n
          WRITE(19,*)(b(i,j),j=1,m)
        ENDDO
        WRITE(19,*)
        DO i=1,m
          WRITE(19,*)(t(i,j), j=1,m)
        ENDDO 
        
        CLOSE(11)
        CLOSE(19)
        END PROGRAM matrixMultiplication
