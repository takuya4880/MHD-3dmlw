
module dt
implicit none 
contains
subroutine detdt(box)
    use defstruct
    !$use omp_lib
    implicit none
    type(cell) :: box[cox,coy,coz,*]

    double precision, allocatable :: v2(:,:,:)
    double precision :: d, dtwav, dtdif, dt(cox,coy,coz)
    integer :: i,j,k
    double precision :: alp=0.01, etamax=1., vc=1000
    double precision :: jx,jy,jz,largest, eta
    allocate(v2(ix,iy,iz))
    
    d = min(box%con%dx,box%con%dy,box%con%dz) 

    !$omp parallel workshare
    v2 = (box%rovx**2 + box%rovy**2 + box%rovz**2)/(box%ro**2)
    v2 = v2 + (box%bx**2 + box%by**2 + box%bz**2)/box%ro
    v2 = v2 + box%con%gam * box%pr / box%ro
    !$omp end parallel workshare
    dtwav = box%con%a * d / sqrt(maxval(v2))

    largest = 1.e-5
    
    do k=1,iz-1
      do j=1,iy-1
        do i=1,ix-1
            jx = (box%bz(i,j+1,k)-box%bz(i,j-1,k))/(2.*box%con%dy) &
                        -(box%by(i,j,k+1)-box%by(i,j,k-1))/(2.*box%con%dz)
            jy = (box%bx(i,j,k+1)-box%bx(i,j,k-1))/(2.*box%con%dz) &
                        -(box%bz(i+1,j,k)-box%bz(i-1,j,k))/(2.*box%con%dx)
            jz = (box%by(i+1,j,k)-box%by(i-1,j,k))/(2.*box%con%dx) &
                        -(box%bx(i,j+1,k)-box%bx(i,j-1,k))/(2.*box%con%dy)

            eta = sqrt((jx**2+jy**2+jz**2)*16.*atan(1.0))/box%ro(i,j,k)  !calculate vd
            
            if (eta<vc) then
                eta = 0.
            else 
                eta = alp*(eta/vc-1.)**2
                if (eta>etamax) then
                    eta = etamax
                end if
            end if

            if (eta>largest) then
                largest = eta
            end if
        end do
      end do
    end do
    dtdif = box%con%a * 0.5 * d**2/largest

    box%con%dt = min(dtwav,dtdif)

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
