        program area_perimeter
        implicit none        
        write(*,*)"This is my 1st fortran program."
        write(*,*)"Give value of a,b:"
        integer::a,b,area,perimeter
        real::c,f
        read(*,*)a,b
        area=a*b
        perimeter=2*(a+b)
        write(*,*)area,perimeter
        end program area_perimeter
