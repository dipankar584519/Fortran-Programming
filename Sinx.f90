      program sine_x
        implicit none
        real::sine,x,y,d,t
        integer::i,n
        write(*,*)"value of x, and n:"
        read(*,*)x,n
        y=3.14*x/180
        sine=y
        t=y
        do i=2,n
          d=(2*i-2)*(2*i-1)
          t=t*(-y)*(y)/d
          sine=sine+t
        enddo
        write(*,*)"Sin(x)=",sine

      end program sine_x



