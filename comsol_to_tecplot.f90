! module to plot the results of the multipole method in tecplot
! read the mesh data from .mphtxt file at path.
! use that data to create a new file at path. You need to compute a vector of n_points dimension called solution

module from_comsol_to_tecplot

contains

subroutine read_comsol_mesh(path,n_points,points_cloud,n_elements,connectivity_matrix)
    implicit none

    integer, intent(out) :: n_points,n_elements
    integer :: i
    double precision,allocatable,dimension(:,:),intent(out) :: points_cloud
    integer, allocatable, dimension(:,:),intent(out) :: connectivity_matrix
    character(50), dimension(17) :: dummy_begining
    character(50), dimension(3) :: dummy_2
    character(50), dimension(9) :: dummy_3
    character(50), intent(in) :: path
    character(50) :: dummy

    open(unit=1,file=path)
    do i=1,17
        read(1,'(A)') dummy_begining(i)
    end do

    read(1,*) n_points
    allocate(points_cloud(n_points,2))

    do i=1,3
        read(1,"(A)") dummy_2(i)
    end do

    do i=1,n_points
        read(1,*) points_cloud(i,1),points_cloud(i,2)
    end do

    do i=1,9
        read(1,"(A)") dummy_3(i)
    end do

    read(1,*) n_elements
    allocate(connectivity_matrix(n_elements,3))

    read(1,"(A)") dummy

    do i=1,n_elements
        read(1,*) connectivity_matrix(i,1),connectivity_matrix(i,2),connectivity_matrix(i,3)
    end do
end subroutine

subroutine write_tecplot_mesh(path,n_points,points_cloud,n_elements,connectivity_matrix,solution)

    implicit none
    integer :: i,j
    integer, intent(in) :: n_points,n_elements
    double precision, dimension(n_points),intent(in) :: solution
    double precision,dimension(n_points,2),intent(in) :: points_cloud
    integer, dimension(n_elements,3),intent(in) :: connectivity_matrix
    character(50),intent(in) :: path
    character(15) :: str_npoints,str_elements

    write(str_npoints,*) n_points
    write(str_elements,*) n_elements

    open(unit=1,file=path)
    write(1,"(A)") "VARIABLES = x, y,<Greek>Q</Greek>"
    write(1,"(A)") "ZONE T='P_1', DATAPACKING=POINT, NODES=" // trim(adjustl(str_npoints)) //", ELEMENTS="//trim(adjustl(str_elements))//(", ZONETYPE=FETRIANGLE")


    do i=1,n_points
    write(1,*) points_cloud(i,1) ,",", points_cloud(i,2) ,"," , solution(i)
    end do

    write(1,*) " "

    do i=1,n_elements
    write(1,"(i4,$)") (connectivity_matrix(i,j)+1,j=1,3)
    write(1,*)" "
    end do

    close(1)
end subroutine

end module