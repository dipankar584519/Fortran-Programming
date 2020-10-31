        Program prime
        implicit none
        integer::i,j,k
        real::x,y,z
        logical::yn=.TRUE.!

        
        write(*,*)"          2 -is a prime number"
        Loop1: DO i=3,100
            Innerdo: DO j=2,i-1
            write(*,*)i,j
                IF(MOD(i,j)==0) CYCLE Loop1!EXIT Loop1
             ENDDO Innerdo               
        write(*,*)i,"-is a prime number"    
        ENDDO Loop1

        End Program prime 
    
