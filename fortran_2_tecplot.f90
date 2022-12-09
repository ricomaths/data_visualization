subroutine oned_totecplot_complex(x,y,path_data)
    implicit none

    real, intent(in) :: x(:)
    complex, intent(in) :: y(:)
    character, intent(in) :: path_data

    integer :: i,n

    n = size(x(:))
    open(unit=1,file=path_data)
    write(1,*) 'VARIABLES = "x","real(y)","imag(y)"'
    do i=1,n
        write(1,*) x(i),real(y(i)),aimag(y(i))
    end do

    close(1)
end subroutine

subroutine oned_totecplot_real(x,y,path_data)
    implicit none

    real, intent(in) :: x(:),y(:)
    character, intent(in) :: path_data

    integer :: i,n

    n = size(x(:))
    open(unit=1,file=path_data)
    write(1,*) 'VARIABLES = "x","y"'
    do i=1,n
        write(1,*) x(i),y(i)
    end do

    close(1)
end subroutine
