!----------------------------------------------------!
!Lab PH707
!Title: Gauss-Jordon Method
!Date: 06/10/2020
!Name: Dipankar Barman
!Roll No: 206121013
!Email: dipankar1998@iitg.ac.in
!----------------------------------------------------!

        PROGRAM MatrixGaussJordon
        IMPLICIT NONE
        INTEGER::i,j,k,m,n,MaxRow
        REAL::pivot,c,store,MaxEl
        REAL,ALLOCATABLE,DIMENSION(:)::B,x
        !matrix of unknown dimension
        REAL,ALLOCATABLE,DIMENSION(:,:)::A,t
        REAL::s
        !opening input and output file
        OPEN(UNIT=11,FILE="MatrixInput1.inp")
        OPEN(UNIT=19,FILE="matrixGJ.out")
        !reading dimension of the matrix from input file
        READ(11,*)n
        !allocating dimension of the matrices
        ALLOCATE(A(n,n+1),t(n,n+1),B(n),x(n))
        !reading matrix A from input file
        DO i=1,n
          READ(11,*)(A(i,j), j=1,n)
        ENDDO
        !reading matrix B, and merging with A matrix
        DO i=1,n
          READ(11,*)B(i)
          A(i,n+1)=B(i)
        ENDDO
        !writing the merged matrix in output file
        WRITE(19,*)"The Augmented matrix is:" 
        DO i=1,n
          WRITE(19,*)(A(i,j), j=1,n+1)
        ENDDO
        WRITE(19,*) "The reduced matrix is:" 
        DO k=1,n
        !performing the row exchange operation 
          MaxEl=abs(A(k,k))
          MaxRow=k
          DO i=k+1,n
           IF(abs(A(i,k))>MaxEl) THEN
             MaxEL=abs(A(i,k))
             MaxRow=i
           END IF
          ENDDO
        IF(MaxRow/=k)THEN 
          DO i=k,n+1
           t(k,i)=A(maxRow,i)
           A(maxRow,i) = A(k,i)
           A(k,i) =t(k,i)
          ENDDO     
         ENDIF 
         !performing elimination process
          DO i=1,n
            IF(i/=k)THEN
              c=A(i,k)/A(k,k)
              DO j=1,n+1
                A(i,j)=A(i,j)-c*A(k,j)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        !writing the final matrix in output file
          DO i=1,n
           WRITE(19,*)(A(i,j)/A(i,i),j=1,n+1)
          ENDDO
        !calculating and writing results
        WRITE(19,*)"Result:"
        DO i=1,n
         WRITE(19,*)"x[",i,"]=",A(i,n+1)/A(i,i)
        ENDDO
        !closing input and output file
        CLOSE(11)
        CLOSE(19)
        END PROGRAM MatrixGaussJordon
!------------------------------------------------------!
!RESULT 
!------------------------------------------------------!
! The Augmented matrix is:
!  0.10  7.00 -0.30-19.30
!  3.00 -0.10 -0.20  7.85
!  0.30 -0.20 10.00 71.40
!The reduced matrix is:
!  1.00  0.00  0.00  3.00
! -0.00  1.00  0.00 -2.50
! -0.00  0.00  1.00  7.00
!Result:
! x[           1 ]=   3.00000000    
! x[           2 ]=  -2.50000000    
! x[           3 ]=   6.99999952    
!-----------------------------------------------------!
