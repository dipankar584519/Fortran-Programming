!####################################################!
!Lab PH707
!Title: Newton-Raphson Method
!Date: 25/09/2020
!Name: Dipankar Barman
!Roll No: 206121013
!Email: dipankar1998@iitg.ac.in
!#####################################################!

      PROGRAM newton_raphson
      implicit none
                         
        integer::i
        real::x,f,df,U
!"The given potential is,U(x)=-14.3996/x+728.0Exp(-x/0.317)-0.486/x^6. We need to find the stable point of the potential near x=1"
        write(*,*)"Enter the initial value"
        read(*,*)x              !reading the initial guess!
        write(*,*)"The Program Precision is set as (x1-x0)<=0.001"
                !write(*,*)"value(x) "," Potential energy"," Value(f/df)"
        Loop1:DO i=1,100        !Maximum No of Iterations is 100!
                U=-14.3996/x+728.0*Exp(-x/0.317)-0.486/x**6 
               ! f=14.3996/x**2-(728.0/0.317)*Exp(-x/0.317)+0.486*6/x**7
                
        f=2*x**3-11.7*x**2+17.7*x-5
        !df=-2*14.3996/x**3+(728/0.317**2)*exp(-x/0.317)-42*0.486/x**8
        df=6*x**2-23.4*x+17.7     
         x=x-(f/df)

               ! write(*,'(f10.4x,f10.4x,xf13.4)')x,U,(f/df)
            IF(abs(f/df)<=0.001)THEN
             EXIT Loop1
            END IF
            END DO Loop1
         write(*,*)"the equilibrium point is:",x !Output Result!

        END PROGRAM newton_raphson
        

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
    
