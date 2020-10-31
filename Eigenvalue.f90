!-------------------------------------------------------!
!Lab PH707 (LG1)
!Title: EigenvaluePowerMethod
!Date: 27/10/2020
!Name: Dipankar Barman
!Roll No: 206121013
!Email: dipankar1998@iitg.ac.in
!-------------------------------------------------------!
!Modified in Lab Class
!-------------------------------------------------------!
        PROGRAM EigenvaluePowerMethod
        IMPLICIT NONE
        REAL::f,ea,l,l1,u,v,ea0,ea1,stopping,Norm1,Norm2
        INTEGER::i,j,k,m,n,h
        REAL,ALLOCATABLE,DIMENSION(:)::x,y,z,p,p1,y1,x1,q1,q,y2,y3,w,w1
        REAL,ALLOCATABLE,DIMENSION(:,:)::A,A_inv,Identity
        OPEN(UNIT=11,FILE="eigenvalue.inp")
        OPEN(UNIT=19,FILE="eigenvalue.out")
        WRITE(*,*)"GIVE THE ORDER OF THE MATRIX:"
        READ(*,*)n     !oredr of the matrix
        WRITE(*,*)"GIVE THE ACCURACY:"
        READ(*,*)stopping     !oredr of the matrix
        ALLOCATE(x(n),y(n),A(n,n),A_inv(n,n),p(1000),p1(1000),x1(n),y1(n),q(1000))
        ALLOCATE(q1(1000),y2(n),y3(n),w(1000),Identity(n,n),w1(1000))
        WRITE(*,*)"A is:"
        DO i=1,n
            READ(11,*)(A(i,j),j=1,n)
            WRITE(*,*)(A(i,j),j=1,n)
        ENDDO
        WRITE(*,*)"A_inv is:"
        CALL Gauss_Jordon_Inversion(A,n,A_inv)       
        DO i=1,n
            WRITE(*,*)(A_inv(i,j),j=1,n)
        ENDDO
        Identity=matmul(A_inv,A)
        !DO i=1,n
         !   WRITE(*,*)(Identity(i,j),j=1,n)
        !ENDDO
        
        WRITE(*,*)"Initial guess:"
        DO i=1,n
            x(i)=1
            x1(i)=1
           !READ(*,*)(x(i),i=1,n)
           WRITE(*,*)x(i)
        ENDDO
        
        y=matmul(A,x)/maxval(matmul(A,x),n)
        y1=matmul(A_inv,x1)/maxval(matmul(A_inv,x1),n)
        WRITE(*,*)"without Coefficient AX is:"
        DO i=1,n
        WRITE(*,*)y(i)
        ENDDO
        WRITE(*,*)"without Coefficient A_inv.X is:"
        DO i=1,n
        WRITE(*,*)y1(i)
        ENDDO
        
     Loop1:DO i=2,1000
        y=matmul(A,y)/maxval(matmul(A,y),n)
        !WRITE(*,*)"without Coefficient AX^{(",i,")} is:"
        !DO j=1,n
        !WRITE(*,*)y(j)
        !ENDDO
        y2=matmul(A,y)
        q(i)=0
        DO j=1,n
        q(i)=q(i)+y2(j)*y(j)
        ENDDO
        w(i)=0
        DO j=1,n
        w(i)=w(i)+y(j)*y(j)
        ENDDO
        p(i)=q(i)/w(i)
        !WRITE(*,*)"Eigenvalue=",p(i)
         u=p(i)
         p(1)=0
         WRITE(*,*)p(i)
        ea0=abs((p(i)-p(i-1))/p(i))*100
        WRITE(*,*)ea0
         WRITE(*,*)"Largest Eigenvalue=",u,"Loop No=",i
        WRITE(*,*)"Del(Largest Eigenvalue)=",ea0
         Norm1=sqrt(sum(y**2))
        WRITE(*,*)"Normalization is:",Norm1
        WRITE(*,*)"Normalized Eigenvector is:"
         DO j=1,n
        WRITE(*,*)y(j)/Norm1
        ENDDO
        IF(ea0<stopping)EXIT Loop1
        ENDDO Loop1
       
        
    Loop2:DO k=2,1000
        y1=matmul(A_inv,y1)/maxval(matmul(A_inv,y1),n)
        !WRITE(*,*)"without Coefficient A_inv.X^{(",i,")} is:"
        !DO j=1,n
        !WRITE(*,*)y1(j)
        !ENDDO
        y3=matmul(A_inv,y1)
        q1(k)=0
        DO h=1,n
        q1(k)=q1(k)+y3(h)*y1(h)
        ENDDO
        w1(k)=0
        DO h=1,n
        w1(k)=w1(k)+y1(h)*y1(h)
        ENDDO
        p1(k)=w1(k)/q1(k)
        !WRITE(*,*)"Eigenvalue=",p1(i)
        v=p1(k)
        ea1=abs((p1(k)-p1(k-1))/p1(k))*100
        WRITE(*,*)"Smallest Eigenvalue=",v,"Loop N0=",k
        WRITE(*,*)"Del(Smallest Eigenvalue)=",ea1
         Norm2=sqrt(sum(y1**2))
         WRITE(*,*)"Normalization is:",Norm2
         WRITE(*,*)"Normalized Eigenvector is:"
         DO h=1,n
        WRITE(*,*)y1(h)/Norm2
        ENDDO      
        IF(ea1<stopping)EXIT Loop2
        ENDDO Loop2
       WRITE(*,*)"FINAL RESULT:"
        WRITE(*,*)"Largest Eigenvalue=",u,"Loop No=",i
        WRITE(*,*)"Del(Largest Eigenvalue)=",ea0
        Norm1=sqrt(sum(y**2))
        WRITE(*,*)"Normalization is:",Norm1
        WRITE(*,*)"Normalized Eigenvector is:"
        DO j=1,n
        WRITE(*,*)y(j)/Norm1
        ENDDO
        WRITE(*,*)"Smallest Eigenvalue=",v,"Loop N0=",k
        WRITE(*,*)"Del(Smallest Eigenvalue)=",ea1
         !Norm2=sqrt(sum(y1**2))
         WRITE(*,*)"Normalization is:",Norm2
         WRITE(*,*)"Normalized Eigenvector is:"
         DO h=1,n
        WRITE(*,*)y1(h)/Norm2
        ENDDO
        
        CLOSE(11)
        CLOSE(19)
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
       

 
