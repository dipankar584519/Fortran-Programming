!-------------------------------------------------------!
!Lab PH707 (LG1)
!Title: EigenvalueShiftOrigin
!Date: 27/10/2020
!Name: Dipankar Barman
!Roll No: 206121013
!Email: dipankar1998@iitg.ac.in
!-------------------------------------------------------!
!-------------------------------------------------------!
        PROGRAM EigenvaluePowerMethod
        IMPLICIT NONE
        REAL::v,ea1,stopping,Norm2,GUESS
        INTEGER::i,j,k,m,n,h
        REAL,ALLOCATABLE,DIMENSION(:)::x,p1,y1,q1,y3,w1
        REAL,ALLOCATABLE,DIMENSION(:,:)::Mm,A,A_inv
        OPEN(UNIT=11,FILE="eigenvalue.inp")
        WRITE(*,*)"GIVE THE ORDER OF THE MATRIX:"
        READ(*,*)n     !oredr of the matrix
       ! WRITE(*,*)"GIVE THE ACCURACY:"
        !READ(*,*)stopping     !oredr of the matrix
        stopping=0.01
        ALLOCATE(x(n),Mm(n,n),A(n,n),A_inv(n,n),p1(1000),y1(n),q1(1000),y3(n),w1(1000))
        WRITE(*,*)"ENTER THE GUESS VALUE:"
        READ(*,*)GUESS
        WRITE(*,*)"Mm is:"
        DO i=1,n
            READ(11,*)(Mm(i,j),j=1,n)
            WRITE(*,*)(Mm(i,j),j=1,n)
        ENDDO
        WRITE(*,*)"A=M-(Identity*GUESS) is:"
        DO i=1,n
            DO j=1,n
                IF(j==i)THEN
                    A(i,i)=Mm(i,i)-GUESS
                ELSE 
                    A(i,j)=Mm(i,j)
                ENDIF
            ENDDO
        ENDDO
        DO i=1,n
            WRITE(*,*)(A(i,j),j=1,n)
        ENDDO
        WRITE(*,*)"A_inv is:"
         CALL Gauss_Jordon_Inversion(A,n,A_inv)       
        
        DO i=1,n
            WRITE(*,*)(A_inv(i,j),j=1,n)
        ENDDO
        WRITE(*,*)"Initial guess:"
        DO i=1,n
            x(i)=1
            WRITE(*,*)x(i)
        ENDDO
 
        y1=matmul(A_inv,x)/maxval(matmul(A_inv,x),n)
        WRITE(*,*)"without Coefficient A_inv.X is:"
        DO i=1,n
         WRITE(*,*)y1(i)
        ENDDO
        
    Loop2:DO k=2,100
        y1=matmul(A_inv,y1)/maxval(matmul(A_inv,y1),n)
        y3=matmul(A_inv,y1)
        q1(k)=0
        DO h=1,n
        q1(k)=q1(k)+y3(h)*y1(h)
        ENDDO
        w1(k)=0
        DO h=1,n
        w1(k)=w1(k)+y1(h)*y1(h)
        ENDDO
        p1(k)=GUESS+w1(k)/q1(k)
        v=p1(k)
        ea1=abs((p1(k)-p1(k-1))/p1(k))*100
        WRITE(*,*)"Eigenvalue=",v,"Loop N0=",k
        WRITE(*,*)"Del(Eigenvalue)=",ea1
         Norm2=sqrt(sum(y1**2))
         WRITE(*,*)"Normalization is:",Norm2
         WRITE(*,*)"Normalized Eigenvector is:"
         DO h=1,n
        WRITE(*,*)y1(h)/Norm2
        ENDDO      
        IF(ea1<stopping)EXIT Loop2
        ENDDO Loop2
        CLOSE(11)
     
        END PROGRAM EigenvaluePowerMethod
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

