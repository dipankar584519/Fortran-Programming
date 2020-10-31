!-----------------------------------------------------!
!Lab PH707
!Title: Gauss Elimination Method
!Date: 06/10/2020
!Name: Dipankar Barman
!Roll No: 206121013
!Email: dipankar1998@iitg.ac.in
!----------------------------------------------------!

        PROGRAM matrixGaussElimination
        IMPLICIT NONE
        INTEGER::i,j,k,m,n,MaxRow
        REAL::pivot,c,store,MaxEl
        REAL,ALLOCATABLE,DIMENSION(:)::B,x
        !matrix of unknown dimension
        REAL,ALLOCATABLE,DIMENSION(:,:)::A,t 
        REAL::s
        !opening input and output file
        OPEN(UNIT=11,FILE="MatrixInput1.inp")
        OPEN(UNIT=19,FILE="matrixGE.out")   
        READ(11,*)n
        !allocating dimension of matrices
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
        WRITE(19,*)"The Augmented matrix is:" 
        !writing the merged matrix in output file
        DO i=1,n
          WRITE(19,*)(A(i,j), j=1,n+1)
        ENDDO
        !recognising max value from coloums 
        DO k=1,n
          MaxEl=abs(A(k,k))               
          MaxRow=k 
          DO i=k+1,n
           IF(abs(A(i,k))>MaxEl) THEN
             MaxEL=abs(A(i,k))
             MaxRow=i
           END IF
          ENDDO
          !row exchange operation
        IF(MaxRow/=k)THEN 
          DO i=k,n+1
           t(k,i)=A(maxRow,i)
           A(maxRow,i) = A(k,i)
           A(k,i) =t(k,i)
          ENDDO     
         ENDIF 
         !performing forward (Gauss) elimination
          DO i=k+1,n
            IF(i/=k)THEN
              c=A(i,k)/A(k,k)
              DO j=1,n+1
                A(i,j)=A(i,j)-c*A(k,j)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        !writing the reduced matrix in output file
           WRITE(19,*) "The reduced matrix is:" 
          DO i=1,n
           WRITE(19,*)(A(i,j),j=1,n+1)
          ENDDO
        !back substitution for solutions, and writing on output file
        x(n)=A(n,n+1)/A(n,n)
        WRITE(19,*)"Result:"
        WRITE(19,*)"x[",n,"]=",x(n) 
        DO i=n-1,1,-1
         s=A(i,n+1)
         DO j=i+1,n
          s=s-A(i,j)*x(j)
         ENDDO
         x(i)=s/A(i,i)
         WRITE(19,*)"x[",i,"]=",x(i)
        ENDDO
        !closing input and output files
        CLOSE(11)
        CLOSE(19)
        END PROGRAM matrixGaussElimination
!-------------------------------------------------!
!RESULT
!-------------------------------------------------!
! The Augmented matrix is:
!  0.10  7.00 -0.30-19.30
!  3.00 -0.10 -0.20  7.85
!  0.30 -0.20 10.00 71.40
!The reduced matrix is:
!  3.00 -0.10 -0.20  7.85
! -0.00  7.00 -0.29-19.56
! -0.00  0.00 10.01 70.08
!Result:
! x[           3 ]=   6.99999952    
! x[           2 ]=  -2.50000000    
! x[           1 ]=   3.00000000    
!-------------------------------------------------!
