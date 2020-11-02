!####################################################!
!Lab PH707
!Title: Bisection Method
!Date: 01/11/2020
!Name: Dipankar Barman
!Roll No: 206121013
!Email: dipankar1998@iitg.ac.in
!#####################################################!

      PROGRAM Bisection
      implicit none
                         
        integer::i
        real::x1,x2,f,x,xr
!"The given potential is,U(x)=-14.3996/x+728.0Exp(-x/0.317)-0.486/x^6. We need to find the stable point of the potential near x=1"
        write(*,*)"Enter the initial two value"
        read(*,*)x1,x2             !reading the initial guess!
        IF(f(x1)*f(x2)<0)THEN
        write(*,*)"There is a root in the given range"
        Loop1:DO i=1,100
            xr=x
            x=(x1+x2)/2
            write(*,*)"f(",x,")*f(",x1,")=",f(x)*f(x1)
            IF(f(x)*f(x1)<0)THEN
             x2=x
             ELSE
             x1=x
            END IF
         write(*,*)"the equilibrium point is:",x ,"f(x)=",f(x)!Output Result!
         WRITE(*,*)"abs((x-xr)/xr)*100=",abs((x-xr)/xr)*100
            IF(abs((x-xr)/xr)*100<0.001)THEN
            EXIT Loop1
            ENDIF
        END DO Loop1
        ELSE
        WRITE(*,*)"There is no root in the given range."
        ENDIF
        END PROGRAM Bisection
        
        REAL FUNCTION f(x)
        REAL,INTENT(IN)::x
        REAL::d
        d=2*x**3-11.7*x**2+17.7*x-5
        f=d
        END FUNCTION f
!------------------------------------------------------!
       ! RESULT!
!------------------------------------------------------!
!Enter the initial value
!1
! The Program Precision is set as (x1-x0)<=0.001
!  value(x)  Potential energy Value(f/df)
!    1.3104    16.1684        -0.3104
!    1.5875     0.5799        -0.2771
!    1.8214    -4.2346        -0.2338
!    1.9836    -5.5921        -0.1622
!    2.0545    -5.8723        -0.0709
!    2.0660    -5.8999        -0.0114
!    2.0662    -5.9004        -0.0003
! the equilibrium point is:   2.06621480   
!------------------------------------------------------!
    
