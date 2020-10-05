        Program matrix
        implicit none
        integer::i,j,k
        real::x,y,a,b
        real,dimension(:)::c(3)
        real,dimension(:,:)::d(3,3),e(3,3),t(3,3)

        OPEN(unit=11,File="matrix.inp")
!       read/write(*,*)    !READING/WRITING FROM TERMINAL SCREEN
!       read/write(x,*)    !READING/WRITING FROM/ON UNIT=x
        OPEN(unit=19,File="matrix.out")
!       read/write(x,*)      !JUST FOR SKIPPING A LINE IN UNIT=x
        DO i=1,3
           read(11,*)(d(i,j), j=1,3)
        ENDDO
        read(11,*)

        DO i=1,3
           read(11,*)(e(i,j), j=1,3)
        ENDDO

        DO i=1,3
           write(19,*)(d(i,j), j=1,3)
        ENDDO
        write(19,*)
        DO i=1,3
           write(19,*)(e(i,j), j=1,3)
        ENDDO
        t=MATMUL(d,e)   !A DECENT FUNTION FOR MATRIX MULTIPLICATION
        write(19,*) 
        DO i=1,3
           write(19,*)(t(i,j), j=1,3)
        ENDDO

        CLOSE(11)
        CLOSE(19)
        End Program matrix


