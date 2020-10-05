        program area_perimeter
        implicit none        
        integer::a,b,area,perimeter
        real::c,f
        read(*,*)a,b
        area=a*b
        perimeter=2*(a+b)
        write(*,*)area,perimeter
    
        end program area_perimeter
