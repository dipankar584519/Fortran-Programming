!-------------------------------------------------------!
!Lab PH707 (LG1)
!Title: Non_Linear_Fitting
!Date: 19/10/2020
!Name: Dipankar Barman
!Roll No: 206121013
!Email: dipankar1998@iitg.ac.in
!-------------------------------------------------------!
        PROGRAM Non_Linear_Fitting
        IMPLICIT NONE
        REAL::f,a0,a1,fx,df_da0,df_da1,ea0,ea1,es,cr,ss,sy
        INTEGER::i,j,k,m,n
        REAL,ALLOCATABLE,DIMENSION(:)::x,y,dy,u,v,g
        REAL,ALLOCATABLE,DIMENSION(:,:)::ZTZ,ZTZ_inv,s,A,Z,ZT
!Input file given in lab will be saved as "Input_NonLinearFit.txt" to avoid any modification after lab class
        OPEN(UNIT=11,FILE="NonLinearFitting.inp")
        OPEN(UNIT=19,FILE="Output_NonLinearFit.txt")
        
        WRITE(*,*)"GIVE NUMBER OF POINTS:"
        READ(*,*)n
        WRITE(*,*)"GIVE THE ACCURACY:"
        READ(*,*)es
        WRITE(*,*)"GIVE THE INITIAL GUESS FOR a0, a1:"
        READ(*,*)a0,a1
        m=2
        ALLOCATE(x(n),y(n),dy(n),g(m),A(n,2),Z(n,m),ZT(m,n),u(m),v(m),ZTZ(m,m),s(m,m),ZTZ_inv(m,m))
        WRITE(19,*)"GIVEN DATA POINTS:"        
        DO i=1,n
            READ(11,*)(A(i,j),j=1,2)
            WRITE(19,*)(A(i,j),j=1,2)
            x(i)=A(i,1)
            y(i)=A(i,2)
        ENDDO
 
 Loop1:DO k=1,100
        WRITE(19,*)
        WRITE(19,*)"LOOP NUMBER=",k 
        !WRITE(19,*)"Z is:"
        DO i=1,n
           Z(i,1)=df_da0(a0,a1,x(i))           
           Z(i,2)=df_da1(a0,a1,x(i))
           !WRITE(19,*)(Z(i,j),j=1,m)
        ENDDO
        ZT=TRANSPOSE(Z)
        
        !WRITE(19,*)"ZTZ is:"
        ZTZ=matmul(ZT,Z)
        !DO i=1,m
        !   WRITE(19,*)(ZTZ(i,j),j=1,m)
        !ENDDO
        Call Gauss_Jordon_Inversion(ZTZ,m,ZTZ_inv)
        
        !WRITE(19,*)"ZTZ_inv IS:"
        !DO i=1,m
        ! WRITE(19,*)(ZTZ_inv(i,j),j=1,m)
        !ENDDO
        !WRITE(19,*)"ZTZ_inv.ZTZ IS:"
        s=matmul(ZTZ_inv,ZTZ)
        !DO i=1,m
        ! WRITE(19,*)(s(i,j),j=1,m)
        !ENDDO
        !WRITE(19,*)"(Y-fx0) is:"
        DO i=1,n
          dy(i)=y(i)-fx(a0,a1,x(i))
          !WRITE(19,*)dy(i)
        ENDDO
        !WRITE(19,*)"ZT(Y-fx0) is:"
        v=matmul(ZT,dy)
        !DO i=1,m
        !WRITE(19,*)v(i)
        !ENDDO
        !WRITE(19,*)"Del(a0,a1)=ZTZ_inv.ZT.(Y-fx0) IS:"
        g=matmul(ZTZ_inv,v)
        !DO i=1,m
        ! WRITE(19,*)g(i)
        !ENDDO
        a0=a0+g(1)
        a1=a1+g(2)
        WRITE(19,*)"a0=",a0,"  a1=",a1
        ea0=ABS(g(1)/a0)*100
        ea1=ABS(g(2)/a1)*100
        WRITE(19,*)"ea0=",ea0,",ea1=",ea1 
        IF(ea0*a0<es.OR.ea1*a1<es)EXIT Loop1
         ENDDO Loop1

        WRITE(*,*)"FINAL RESULT:"
        WRITE(*,*)"a0=",a0,",and a1=",a1,"Loop No:",k 
        WRITE(19,*)"FINAL RESULT: a0=",a0,",and a1=",a1 
        WRITE(*,*)"ea0=",ea0,",ea1=",ea1 
        CALL least_square_fit(a0,a1,x,y,n,cr,ss,sy)
        WRITE(*,*)"GOODNESS OF FITTING:"
        WRITE(*,*)"S=",ss,",S_y=",sy,",C_r=",cr
        WRITE(19,*)"GOODNESS OF FITTING:"
        WRITE(19,*)"S=",ss,",S_y=",sy,",C_r=",cr
        
        CLOSE(11)
        CLOSE(19)
        END PROGRAM Non_Linear_Fitting
        
!------------------FUNCTION fx------------------------!   
        REAL FUNCTION fx(a0,a1,xi)
        IMPLICIT NONE
        REAL,INTENT(IN)::a0,a1,xi
        fx=a0*(1-exp(-a1*xi))
        END FUNCTION fx
        
!------------------FUNCTION df_da0--------------------!
        REAL FUNCTION df_da0(a0,a1,xi)
        IMPLICIT NONE
        REAL,INTENT(IN)::a0,a1,xi    
        df_da0=(1-exp(-a1*xi))
        END FUNCTION df_da0
        
!-------------------FUNCTION df_da1-------------------!
        REAL FUNCTION df_da1(a0,a1,xi)
        IMPLICIT NONE
        REAL,INTENT(IN)::a0,a1,xi
        df_da1=a0*xi*exp(-a1*xi)
        END FUNCTION df_da1
        
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
         
!------------SUBROUTINE CORRELATION FUNCTION--------------!   
 SUBROUTINE least_square_fit(a0,a1,x,y,n,cr,ss,sy)
        IMPLICIT NONE
        INTEGER::i,j,k,n
        REAL::d,cr,sy,ss,a0,a1,fx
        REAL,DIMENSION(:)::x(n),y(n),s(n),l(n),q(n)
                
        d=SUM(y)/n
        DO i=1,n
!storing function values in an array (fx is a function defined earlier!
        l(i)=fx(a0,a1,x(i))     
        s(i)=(y(i)-l(i))**2         !(y_i-f(p;x_i))^2
        q(i)=(y(i)-d)**2            !(y_i-y_{avg})^2
        ENDDO
        ss=SUM(s)
        sy=SUM(q)
        cr=SQRT((sy-ss)/sy)
        END SUBROUTINE least_square_fit
 
