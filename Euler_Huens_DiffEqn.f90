        PROGRAM DIFFEQN
        implicit NONE
        REAL::xi,xf,y0,dx,dy,f,ff
        INTEGER::i,j,k,n 
        REAL,ALLOCATABLE,DIMENSION(:)::x,y,y1
        !WRITE(*,*)"n, xi,xf, y0"
        !READ(*,*)n,xi,xf,y0
        y0=2; xi=0; xf=4; n=21
        ALLOCATE(x(n),y(n),y1(n))
        dx=(xf-xi)/(n-1)
        x(1)=xi
        x(n)=xf
        DO i=2,n-1
        x(i)=x(i-1)+dx
        ENDDO
WRITE(*,*)"     x(i),        y_Euler(i),      y_Heuns(i)       True y(i)"
        CALL Euler(x,y0,y,n,dx)
        CALL HUENS(x,y0,y1,n,dx)
        DO i=1,n
        WRITE(*,*)x(i), y(i), y1(i), ff(x(i))
        ENDDO
        
        END PROGRAM DIFFEQN
        
    REAL FUNCTION f(x,y)
        REAL,INTENT(IN)::x,y
        f=4*Exp(0.8*x)-0.5*y
        END FUNCTION f
        
    REAL FUNCTION ff(x)
        REAL,INTENT(IN)::x
        ff=4*(Exp(0.8*x)-exp(-0.5*x))/1.3+2*Exp(-0.5*x)
        END FUNCTION ff
    
    SUBROUTINE Euler(x,y0,y,n,dx)
    implicit NONE
       REAL,INTENT(IN)::y0,dx
       INTEGER,INTENT(IN)::n
       REAL,DIMENSION(:)::x(n),y(n)
       INTEGER::i 
       REAL::f
       y(1)=y0
        DO i=1,n-1
        y(i+1)=y(i)+f(x(i),y(i))*dx
        ENDDO
       END SUBROUTINE Euler
       
    SUBROUTINE HUENS(x,y0,y1,n,dx)
        implicit NONE
        REAL,INTENT(IN)::y0,dx
        INTEGER,INTENT(IN)::n
        REAL,DIMENSION(:)::x(n),y1(n),y(n)
        REAL::f
        INTEGER::i 
        y1(1)=y0
        y(1)=y0
        DO i=1,n-1
        y(i+1)=y(i)+f(x(i),y(i))*dx
        y1(i+1)=y1(i)+(f(x(i),y(i))+f(x(i+1),y(i+1)))*dx/2
        ENDDO
       END SUBROUTINE HUENS
