        PROGRAM MatrixInversion
        IMPLICIT NONE
        INTEGER::i,j,k,m,n,MaxRow
        REAL::pivot,c,store,MaxEl
        REAL,ALLOCATABLE,DIMENSION(:,:)::A,B,Id,t,INV !MATRIX OF UNKNOWN DIMENSION
        REAL::s
        OPEN(UNIT=11,FILE="MatrixInput.inp")
        OPEN(UNIT=19,FILE="matrix1.out")
        READ(11,*)n
        m=2*n
        ALLOCATE(A(n,m),B(n,n),t(n,m),Id(n,n),INV(n,n))
        DO i=1,n
          READ(11,*)(B(i,j), j=1,n)
        ENDDO
        DO i=1,n
          DO j=1,n
            IF(i/=j)THEN
            Id(i,j)=0
            ELSE 
            Id(i,j)=1
            ENDIF
            k=n+j
            A(i,j)=B(i,j)
            A(i,k)=Id(i,j)
          ENDDO
        ENDDO
        
        WRITE(19,*)"The matrix is:" 
        DO i=1,n
          WRITE(19,'(12(f6.2))')(A(i,j), j=1,n)
        ENDDO
        !WRITE(19,*)"After Row Reductio Operations:"
 
         DO k=1,n
          MaxEl=abs(A(k,k))
          MaxRow=k
          !WRITE(19,*)MaxEl
          DO i=k+1,n
           IF(abs(A(i,k))>MaxEl) THEN
             MaxEL=abs(A(i,k))
             MaxRow=i
            !write(19,*)MaxEl,MaxRow
           END IF
          ENDDO
        IF(MaxRow/=k)THEN 
          DO i=k,m
            t(k,i)=A(maxRow,i)
            A(maxRow,i) = A(k,i)
            A(k,i) =t(k,i)
          ENDDO     
         ENDIF 
         !DO i=1,n
          ! WRITE(19,'(6(f6.2))')(A(i,j), j=1,m)
         !ENDDO 
         DO i=1,n
            IF(i/=k)THEN
              c=A(i,k)/A(k,k)
              DO j=1,m
                A(i,j)=A(i,j)-c*A(k,j)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        !  DO i=1,n
         !  WRITE(19,'(6(f10.6))')(A(i,j)/A(i,i),j=1,m)
          !ENDDO
        WRITE(19,*)"THE INVERSE OF THE GIVEN MATRIX IS:"
        DO i=1,n
          DO j=1,n
            k=n+j
            INV(i,j)=A(i,k)/A(i,i)
          ENDDO
        ENDDO


          DO i=1,n
           WRITE(19,'(6(f10.6))')(INV(i,j),j=1,n)
          ENDDO
          t=MATMUL(INV,B)
          
        WRITE(19,*)"THE MULTIPLICATION OF INV(A) AND A IS:"
          DO i=1,n
           WRITE(19,'(6(f10.6))')(t(i,j),j=1,n)
          ENDDO 



        CLOSE(11)
        CLOSE(19)
        END PROGRAM MatrixInversion
