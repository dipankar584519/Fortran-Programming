!-------------------------------------------------------!
!Lab PH707 (LG1)
!Title: NewtonDividedInterpolation
!Date: 31/10/2020
!Name: Dipankar Barman
!Roll No: 206121013
!Email: dipankar1998@iitg.ac.in
!-------------------------------------------------------!
!Modified in Lab Class
!-------------------------------------------------------!
        PROGRAM NewtonDividedInterpolation
        IMPLICIT NONE
        REAL::f,xi
        INTEGER::i,j,n
        REAL,ALLOCATABLE,DIMENSION(:)::x,y,ea,yint
        REAL,ALLOCATABLE,DIMENSION(:,:)::Fdd
        OPEN(UNIT=11,FILE="interpolation.inp")
        READ(11,*)n     !number of Points
        WRITE(*,*)"GIVE THE NUMBER OF POINTS:",n
        ALLOCATE(x(n),y(n),Fdd(n,n),yint(n),ea(n))
        WRITE(*,*)"Give The Point:"
        READ(*,*)xi
        WRITE(*,*)"Given data set:"
        WRITE(*,*)"       x                y     "
        DO i=1,n
            READ(11,*)x(i),y(i)
            WRITE(*,*)x(i),y(i)
        ENDDO
        WRITE(*,*)"A_inv is:"
        CALL NewtonInterpolation(x,y,n,xi,yint,Fdd,ea)  
        WRITE(*,*)"fdd"
        DO i=1,n
        WRITE(*,*)(Fdd(i,j),j=1,n)
        ENDDO
        WRITE(*,*)"      Y             Del(Y)       "
        DO i=1,n
            WRITE(*,*)yint(i),ea(i)
        ENDDO
        
        CLOSE(11)
        END PROGRAM NewtonDividedInterpolation
!-------------SUBROUTINE FOR INVERSE MATRIX-----------!
        SUBROUTINE NewtonInterpolation(x,y,n,xi,yint,fdd,ea)
        IMPLICIT NONE
        INTEGER,INTENT(IN)::n
        INTEGER::i,j,k
        REAL,INTENT(IN)::xi
        REAL,INTENT(IN),DIMENSION(:)::x(n),y(n)
        REAL,INTENT(OUT),DIMENSION(:)::yint(n),ea(n)
        REAL,DIMENSION(:,:)::fdd(n,n)
        REAL::xterm,yint2
        DO i=1,n
        fdd(i,1)=y(i)
        ENDDO
        DO j=2,n
        DO i=1,n-j+1
        fdd(i,j)=(fdd(i+1,j-1)-fdd(i,j-1))/(x(i+j-1)-x(i))
        ENDDO
        ENDDO
        xterm=1
        yint(1)=fdd(1,1)
        DO i=2,n
        xterm=xterm*(xi-x(i-1))
        yint2=yint(i-1)+fdd(1,i)*xterm
        ea(i-1)=yint2-yint(i-1)
        yint(i)=yint2
        ENDDO
        END SUBROUTINE NewtonInterpolation

 
