        Program prime
        implicit none
        integer::i,j,k
        real::x,y,z
        logical::yn=.TRUE.!

        
        write(*,*)"          2 -is a prime number"
        DO i=3,100
        yn=.TRUE.
                DO j=2,i-1
                !write(*,*)i,j
                    IF(MOD(i,j)==0)THEN
                       yn=.FALSE.
                      exit
                    ENDIF
                ENDDO                
            IF(yn)write(*,*)i,"-is a prime number"    
        ENDDO

        End Program prime 
    
