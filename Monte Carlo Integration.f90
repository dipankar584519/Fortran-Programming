    PROGRAM Exp_x     
    IMPLICIT NONE
    INTEGER::i,j,k,n,hit_miss
    DOUBLE PRECISION::r,x,y,y1,x1,f,random,A_r,importance,analytical,A_i,A_h,E_i,E_h,E_r,z,y2,g
    DOUBLE PRECISION,PARAMETER::e=2.71828182
    random=0; hit_miss=0; importance=0
    WRITE(*,*)"ENTER n:"
    READ(*,*)n
    analytical=e-1.0
    
    DO i=1,n
    CALL RANDOM_NUMBER(x)
    CALL RANDOM_NUMBER(y)
!Random x
    random=random+f(x)
!Hit&Miss Algorithm
    y1=y*e;x1=x
    IF(y1<f(x1))THEN
    hit_miss=hit_miss+1
    ENDIF
!importance Sampling
    y2=3.0*x/2.0
    z=SQRT(1+2*y2)
    g=Exp(z)/(e*z)
    importance=importance+g
    ENDDO
    
    A_i=3.0*importance/(2.0*REAL(n))
    E_i=ABS((analytical-A_i)/analytical)*100
    A_h=e*REAL(hit_miss)/REAL(n)
    E_h=ABS((analytical-A_h)/analytical)*100
    A_r=random/REAL(n)
    E_r=ABS((analytical-A_r)/analytical)*100
    WRITE(*,*)"Analytical value=",analytical
    WRITE(*,*)"Integral_r=",A_r,"E_r=",E_r
    WRITE(*,*)"Integral_h=",A_h,"E_h=",E_h
    WRITE(*,*)"Integral_i=",A_i,"E_i=",E_i

    END PROGRAM Exp_x
    
    DOUBLE PRECISION FUNCTION f(x)
    DOUBLE PRECISION,INTENT(IN)::x
    f=Exp(x)
    END FUNCTION f
    
   
