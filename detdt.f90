
module dt
implicit none 
contains
subroutine detdt(box)
    use defstruct
    !$use omp_lib
    implicit none
    type(cell) :: box[cox,coy,coz,*]

    double precision, allocatable :: v2(:,:,:)
    double precision :: d, dt(cox,coy,coz)
    integer :: i,j,k
    allocate(v2(ix,iy,iz))
    
    d = min(box%con%dx,box%con%dy,box%con%dz) 

    !$omp parallel workshare
    v2 = (box%rovx**2 + box%rovy**2 + box%rovz**2)/(box%ro**2)
    v2 = v2 + (box%bx**2 + box%by**2 + box%bz**2)/box%ro
    v2 = v2 + box%con%gam * box%pr / box%ro
    !$omp end parallel workshare

    box%con%dt = box%con%a * d / sqrt(maxval(v2))
    sync all
    do i=1,cox
        do j=1,coy 
            do k=1,coz
                dt(i,j,k)=box[i,j,k,1]%con%dt
            end do
        end do
    end do
    sync all
    box%con%dt = minval(dt)

    deallocate(v2)
end subroutine

end module
