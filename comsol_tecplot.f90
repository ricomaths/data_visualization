! module to plot the results of the multipole method in tecplot
! read the mesh data from .mphtxt file at path.
! use that data to create a new file at path. You need to compute a vector of n_points dimension called solution

module comsol_tecplot_functions

implicit none

contains

! ! It receives the path where the mesh from comsol is and return all info about it
subroutine read_comsol_mesh(path,n_points,points_cloud,n_elements,connectivity_matrix)
    implicit none

    integer, intent(out) :: n_points,n_elements
    integer :: i
    real,allocatable,dimension(:,:),intent(out) :: points_cloud
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

! ! ! It receives the mesh file and the solution file from comsol and creates a new file to read in tecplot
subroutine comsol_to_tecplot_solution(path_mesh,path_solution,path_tecplot)

    implicit none

    character(50), intent(in) :: path_mesh,path_solution,path_tecplot
    integer :: n_points,n_elements,i,j
    real, allocatable :: points_cloud(:,:),theta_real(:),theta_imag(:),points_dummy(:,:)
    integer, allocatable :: connectivity_matrix(:,:)
    character(50), dimension(9) :: dummy_beginning
    character(50) :: str_npoints,str_elements
    integer :: indice
    call read_comsol_mesh(path_mesh,n_points,points_cloud,n_elements,connectivity_matrix)

    open(unit=1,file=path_solution)
    do i=1,9
        read(1,'(A)') dummy_beginning(i)
    end do

    allocate(points_dummy(n_points,n_points))
    allocate(theta_real(n_points),theta_imag(n_points))

    do i=1,n_points
        read(1,*) points_dummy(i,1),points_dummy(i,2),theta_real(i),theta_imag(i)
    end do
    close(1)
    
    ! ! writting
    write(str_npoints,*) n_points
    write(str_elements,*) n_elements
    open(unit=2,file=path_tecplot)
    write(2,"(A)") 'VARIABLES = "x","y","Re(<Greek>Q</Greek>)","Im(<Greek>Q</Greek>)"'
    write(2,"(A)") "ZONE T = 'P_1', DATAPACKING = POINT, NODES=" // trim(adjustl(str_npoints)) // ", ELEMENTS=" // trim(adjustl(str_elements))&
    //(",ZONETYPE=FETRIANGLE")

    do i=1,n_points
    indice = minloc(abs(points_cloud(i,1)-points_dummy(:,1))+abs(points_cloud(i,2)-points_dummy(:,2)),dim=1)
    write(2,*) points_cloud(i,1),",",points_cloud(i,2),",",theta_real(indice),",",theta_imag(indice)
    end do

    write(2,*) " "

    do i =1, n_elements
    write(2,"(i6,$)") (connectivity_matrix(i,j)+1,j=1,3)
    write(2,*) " "
    end do

    close(2)

end subroutine

end module