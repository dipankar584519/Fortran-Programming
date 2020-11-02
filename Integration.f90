        PROGRAM integration
        REAL::f,x0,xn,Simpson,trapezoidal,I,I1
        INTEGER::j,k,n,n2
        REAL,ALLOCATABLE,DIMENSION(:)::x
        WRITE(*,*)"n,x0,xn"
        READ(*,*)n,x0,xn
        ALLOCATE(x(n))
        n2=(n+1)/2
        WRITE(*,*)"trapezoidal integration, I=",trapezoidal(x0,xn,n)
        I=(4*trapezoidal(x0,xn,n)-trapezoidal(x0,xn,n2))/3
              
        IF(MOD(n,2)==0)THEN
        WRITE(*,*)"The interval isn't appropriate for Simpson's 1/3rd rule!"
        WRITE(*,*)"The interval isn't appropriate for calculating O(h^2)-order of trapezoidal rule!"
        ELSE
        WRITE(*,*)"I=I+E=",I,"E=",abs(I-trapezoidal(x0,xn,n))
        WRITE(*,*)"Simpson's 1/3rd rule, I=", Simpson(x0,xn,n)
        !I1=(16*Simpson(x0,xn,n)-Simpson(x0,xn,n2))/15
        !WRITE(*,*)"I=I+E=",I1,"E=",abs(I1-Simpson(x0,xn,n))
        ENDIF
        END PROGRAM integration
        
    REAL FUNCTION f(x)
        REAL,INTENT(IN)::x
    
        f=5+2*x+3*x**2+4*x**3+5*x**4+6*x**5+7*x**6
        END FUNCTION f
        
    REAL FUNCTION trapezoidal(x0,xn,n)
        REAL,INTENT(IN)::x0,xn
        INTEGER,INTENT(IN)::n
        REAL,DIMENSION(:)::x(n)
        REAL::h,f,s
        h=(xn-x0)/(n-1)
        x(1)=x0
        x(n)=xn
        DO j=2,n-1
        x(j)=x(j-1)+h 
        ENDDO
        s=f(x(1))+f(x(n))
        DO j=2,n-1
        s=s+2*f(x(j))
        ENDDO
        trapezoidal=h*s/2
        END FUNCTION trapezoidal

        
    REAL FUNCTION Simpson(x0,xn,n)
        REAL,INTENT(IN)::x0,xn
        INTEGER,INTENT(IN)::n
        REAL,DIMENSION(:)::x(n)
        REAL::h,s1,f 
        h=(xn-x0)/(n-1)
        x(1)=x0
        x(n)=xn
        DO j=2,n-1
        x(j)=x(j-1)+h 
        ENDDO
        s1=f(x(1))+f(x(n))
        DO j=2,n-1,2
        s1=s1+4*f(x(j))
        ENDDO
        DO j=3,n-2,2
        s1=s1+2*f(x(j))
        ENDDO
        Simpson=h*s1/3
        END FUNCTION Simpson
