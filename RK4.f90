        PROGRAM DIFFEQN
        implicit NONE
        REAL::xi,xf,y0,dx,f,ff
        INTEGER::i,j,k,n 
        REAL,ALLOCATABLE,DIMENSION(:)::x,y,k1,k2,k3,k4
        !WRITE(*,*)"n, xi,xf, y0"
        !READ(*,*)n,xi,xf,y0
        y0=1200; xi=0; xf=480; n=13
        ALLOCATE(x(n),y(n),k1(n),k2(n),k3(n),k4(n))
        dx=(xf-xi)/(n-1)
        x(1)=xi
        x(n)=xf
        DO i=2,n-1
        x(i)=x(i-1)+dx
        ENDDO
!WRITE(*,*)"     x(i),        y_Euler(i),      y_Heuns(i)       True y(i)"
        CALL RK4(x,y0,y,k1,k2,k3,k4,n,dx)
        DO i=1,n
        WRITE(*,*)x(i),y(i),k1(i),k2(i),k3(i),k4(i)
        ENDDO
        
        END PROGRAM DIFFEQN
        
    REAL FUNCTION f(x,y)
        REAL,INTENT(IN)::x,y
        !f=4*Exp(0.8*x)-0.5*y
        f=-2.2067E-12*(y**4-81E8)
        END FUNCTION f
        
    REAL FUNCTION ff(x)
        REAL,INTENT(IN)::x
        ff=4*(Exp(0.8*x)-exp(-0.5*x))/1.3+2*Exp(-0.5*x)
        END FUNCTION ff
    
    SUBROUTINE RK4(x,y0,y,k1,k2,k3,k4,n,dx)
        implicit NONE
        REAL,INTENT(IN)::y0,dx
        INTEGER,INTENT(IN)::n
        REAL,DIMENSION(:)::x(n),y(n),k1(n),k2(n),k3(n),k4(n)
        REAL::f
        INTEGER::i 
        y(1)=y0
        
        DO i=1,n-1
        k1(i)=f(x(i),y(i))
        k2(i)=f(x(i)+dx/2,y(i)+k1(i)*dx/2)
        k3(i)=f(x(i)+dx/2,y(i)+k2(i)*dx/2)
        k4(i)=f(x(i)+dx,y(i)+k3(i)*dx)
        y(i+1)=y(i)+(k1(i)+2*k2(i)+2*k3(i)+k4(i))*dx/6
        ENDDO
       END SUBROUTINE Rk4
