!----------------------------------------------------!
!Lab PH707
!Title: Euler_Heuns_Midpoint Method
!Date: 11/11/2020
!Name: Dipankar Barman
!Roll No: 206121013
!Email: dipankar1998@iitg.ac.in
!----------------------------------------------------!
        PROGRAM Euler_Heuns_Midpoint
        implicit NONE
        DOUBLE PRECISION::xi,xf,y0,dx
        DOUBLE PRECISION::f,ff
        INTEGER::i,j,k,n 
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::x,y,y1,y2,y3
        OPEN(UNIT=19,FILE="Assignment8.txt")
        WRITE(*,*)"Initial point:"
        READ(*,*)xi
        WRITE(*,*)"Final point:"
        READ(*,*)xf
        WRITE(*,*)"Number of points:"
        READ(*,*)n
        WRITE(*,*)"Value of y at initial point:"
        READ(*,*)y0
!        y0=0.2; xi=0; xf=2; n=11
        ALLOCATE(x(n),y(n),y1(n),y3(n),y2(n))
        dx=(xf-xi)/(n-1)
        !dx=0.5
        x(1)=xi
        x(n)=xf
        DO i=2,n-1
        x(i)=x(i-1)+dx
        ENDDO
        WRITE(*,*)"Step Size:"
        WRITE(*,'(f10.4)')dx
        WRITE(19,*)"Initial point:",xi
        WRITE(19,*)"Final point:",xf
        WRITE(19,*)"Number of points:",n
        WRITE(19,*)"Value of y at initial point:",y0
        WRITE(19,*)"Step Size:",dx
        
 WRITE(*,*)"      x       y_Euler   y_Heuns   y_Midpoint"
WRITE(19,*)"      x       y_Euler   y_Heuns   y_Midpoint"

        CALL Euler(x,y0,y,n,dx)
        CALL HEUNS(x,y0,y1,n,dx)
        CALL MIDPOINT(x,y0,y3,y2,n,dx)
        DO i=1,n
        WRITE(*,'(f10.4,xf10.4,xf10.4,xf10.4,xf10.4,xf14.4)')x(i), y(i),y1(i),y2(i)
        WRITE(19,'(f10.4,xf10.4,xf10.4,xf10.4,xf10.4,xf14.4)')x(i), y(i),y1(i),y2(i)
        ENDDO
        CLOSE(19)
        END PROGRAM Euler_Heuns_Midpoint
   
    DOUBLE PRECISION FUNCTION f(x,y)
        DOUBLE PRECISION,INTENT(IN)::x,y
        !f=-2*x**3+12*x**2-20*x+8.5
        f=x*y-0.1*x*y*y
        END FUNCTION f
        
   DOUBLE PRECISION FUNCTION ff(x,y)
        DOUBLE PRECISION,INTENT(IN)::x,y
        ff=10/(1-Exp(-x*x/2+3.89182))
        END FUNCTION ff
    
    SUBROUTINE Euler(x,y0,y,n,dx)
    implicit NONE
       DOUBLE PRECISION,INTENT(IN)::y0,dx
       INTEGER,INTENT(IN)::n
       DOUBLE PRECISION,DIMENSION(:)::x(n),y(n)
       INTEGER::i 
       DOUBLE PRECISION::f
       y(1)=y0
        DO i=1,n-1
        y(i+1)=y(i)+f(x(i),y(i))*dx
        ENDDO
       END SUBROUTINE Euler
       
    SUBROUTINE HEUNS(x,y0,y,n,dx)
        implicit NONE
        DOUBLE PRECISION,INTENT(IN)::y0,dx
        INTEGER,INTENT(IN)::n
        DOUBLE PRECISION,DIMENSION(:)::x(n),y1(n),y(n)
        DOUBLE PRECISION::f
        INTEGER::i 
        y1(1)=y0
        y(1)=y0
        DO i=1,n-1
        y1(i+1)=y(i)+f(x(i),y(i))*dx
        y(i+1)=y(i)+(f(x(i),y(i))+f(x(i+1),y1(i+1)))*dx/2
        ENDDO
       END SUBROUTINE HEUNS

    SUBROUTINE MIDPOINT(x,y0,y,y1,n,dx)
        implicit NONE
        DOUBLE PRECISION,INTENT(IN)::y0,dx
        INTEGER,INTENT(IN)::n
        DOUBLE PRECISION,DIMENSION(:)::x(n),y1(n),y(n)
        DOUBLE PRECISION::f
        INTEGER::i 
        y1(1)=y0
        y(1)=y0
        DO i=1,n-1
        y(i+1/2)=y1(i)+f(x(i),y1(i))*dx/2
        y1(i+1)=y1(i)+f(x(i)+dx/2,y(i+1/2))*dx
        ENDDO
       END SUBROUTINE MIDPOINT
