module bc
implicit none 
contains
subroutine boundary(box, uboundary)
    use defstruct
    implicit none
    type(cell) :: box[cox,coy,coz,*]
    double precision :: uboundary(9,marg)
    integer :: imz

    imz = box%con%imz

    call periodicbc(box)    

    if (imz==coz) then
        call upgradbc(box%ro, 1, uboundary)
        call upgradbc(box%rovx, 2, uboundary)
        call upgradbc(box%rovy, 3, uboundary)
        call upgradbc(box%rovz, 4, uboundary)
        call upgradbc(box%bx, 5, uboundary)
        call upgradbc(box%by, 6, uboundary)
        call upgradbc(box%bz, 7, uboundary)
        call upgradbc(box%e, 8, uboundary)
        call upgradbc(box%pr, 9, uboundary)
    end if

    if (imz==1) then
        call lowmrbc(box%ro)
        call lowmrbc(box%rovx)
        call lowmrbc(box%rovy)
        call lowmrbc2(box%rovz)
        call lowmrbc(box%bx)
        call lowmrbc(box%by)
        call lowmrbc2(box%bz)
        call lowmrbc(box%e)
        call lowmrbc(box%pr)
    end if


end subroutine 

subroutine upgradbc(arr, k, ub)   !gradient bc for upper boundary
    use defstruct
    implicit none
    double precision :: arr(ix,iy,iz)
    integer :: k
    double precision :: ub(9,marg)

    integer :: i
    
    do i=1,marg
        arr(:,:,iz-marg+i) = arr(:,:,iz-marg) + ub(k,i)
    end do    

end subroutine

subroutine lowmrbc(arr)   !mirror bc for lower boundary
    use defstruct
    implicit none
    double precision :: arr(ix,iy,iz)

    integer :: i
    
    do i=1,marg
        arr(:,:,i) = arr(:,:,2*marg+1-i)   
    end do    

end subroutine


subroutine lowmrbc2(arr)   !mirror bc for lower boundary
    use defstruct
    implicit none
    double precision :: arr(ix,iy,iz)

    integer :: i
    
    do i=1,marg
        arr(:,:,i) = - arr(:,:,2*marg+1-i)   
    end do    

end subroutine

subroutine periodicbc(box) 
    use defstruct
    implicit none
    type(cell) :: box[cox,coy,coz,*]
    
    integer :: imx, imy, imz  
    integer :: right, left, up, down, behind, before

    imx = box%con%imx
    imy = box%con%imy
    imz = box%con%imz
    
    right = imx+1
    left = imx-1
    behind = imy+1
    before = imy-1
    up = imz+1
    down = imz-1

    if (imx==1) then
        left = cox
    end if
    if (imx==cox) then
        right = 1
    end if
    if (imy==1) then
        before = coy
    end if
    if (imy==coy) then
        behind = 1
    end if
    if (imz==1) then
        down = 0
    end if 
    if (imz==coz) then
        up = 0
    end if

    !x direction
    box%ro(1:marg,:,:) = box[left,imy,imz,1]%ro(ix-2*marg+1:ix-marg,:,:) 
    box%rovx(1:marg,:,:) = box[left,imy,imz,1]%rovx(ix-2*marg+1:ix-marg,:,:) 
    box%rovy(1:marg,:,:) = box[left,imy,imz,1]%rovy(ix-2*marg+1:ix-marg,:,:) 
    box%rovz(1:marg,:,:) = box[left,imy,imz,1]%rovz(ix-2*marg+1:ix-marg,:,:) 
    box%bx(1:marg,:,:) = box[left,imy,imz,1]%bx(ix-2*marg+1:ix-marg,:,:) 
    box%by(1:marg,:,:) = box[left,imy,imz,1]%by(ix-2*marg+1:ix-marg,:,:) 
    box%bz(1:marg,:,:) = box[left,imy,imz,1]%bz(ix-2*marg+1:ix-marg,:,:) 
    box%e(1:marg,:,:) = box[left,imy,imz,1]%e(ix-2*marg+1:ix-marg,:,:) 
    box%pr(1:marg,:,:) = box[left,imy,imz,1]%pr(ix-2*marg+1:ix-marg,:,:) 

    box%ro(ix-marg+1:ix,:,:) = box[right,imy,imz,1]%ro(marg+1:2*marg,:,:)
    box%rovx(ix-marg+1:ix,:,:) = box[right,imy,imz,1]%rovx(marg+1:2*marg,:,:)
    box%rovy(ix-marg+1:ix,:,:) = box[right,imy,imz,1]%rovy(marg+1:2*marg,:,:)
    box%rovz(ix-marg+1:ix,:,:) = box[right,imy,imz,1]%rovz(marg+1:2*marg,:,:)
    box%bx(ix-marg+1:ix,:,:) = box[right,imy,imz,1]%bx(marg+1:2*marg,:,:)
    box%by(ix-marg+1:ix,:,:) = box[right,imy,imz,1]%by(marg+1:2*marg,:,:)
    box%bz(ix-marg+1:ix,:,:) = box[right,imy,imz,1]%bz(marg+1:2*marg,:,:)
    box%e(ix-marg+1:ix,:,:) = box[right,imy,imz,1]%e(marg+1:2*marg,:,:)
    box%pr(ix-marg+1:ix,:,:) = box[right,imy,imz,1]%pr(marg+1:2*marg,:,:)

    !y direction
    box%ro(:,1:marg,:) = box[imx,before,imz,1]%ro(:,iy-2*marg+1:iy-marg,:) 
    box%rovx(:,1:marg,:) = box[imx,before,imz,1]%rovx(:,iy-2*marg+1:iy-marg,:) 
    box%rovy(:,1:marg,:) = box[imx,before,imz,1]%rovy(:,iy-2*marg+1:iy-marg,:) 
    box%rovz(:,1:marg,:) = box[imx,before,imz,1]%rovz(:,iy-2*marg+1:iy-marg,:) 
    box%bx(:,1:marg,:) = box[imx,before,imz,1]%bx(:,iy-2*marg+1:iy-marg,:) 
    box%by(:,1:marg,:) = box[imx,before,imz,1]%by(:,iy-2*marg+1:iy-marg,:) 
    box%bz(:,1:marg,:) = box[imx,before,imz,1]%bz(:,iy-2*marg+1:iy-marg,:) 
    box%e(:,1:marg,:) = box[imx,before,imz,1]%e(:,iy-2*marg+1:iy-marg,:) 
    box%pr(:,1:marg,:) = box[imx,before,imz,1]%pr(:,iy-2*marg+1:iy-marg,:) 

    box%ro(:,iy-marg+1:iy,:) = box[imx,behind,imz,1]%ro(:,marg+1:2*marg,:)
    box%rovx(:,iy-marg+1:iy,:) = box[imx,behind,imz,1]%rovx(:,marg+1:2*marg,:)
    box%rovy(:,iy-marg+1:iy,:) = box[imx,behind,imz,1]%rovy(:,marg+1:2*marg,:)
    box%rovz(:,iy-marg+1:iy,:) = box[imx,behind,imz,1]%rovz(:,marg+1:2*marg,:)
    box%bx(:,iy-marg+1:iy,:) = box[imx,behind,imz,1]%bx(:,marg+1:2*marg,:)
    box%by(:,iy-marg+1:iy,:) = box[imx,behind,imz,1]%by(:,marg+1:2*marg,:)
    box%bz(:,iy-marg+1:iy,:) = box[imx,behind,imz,1]%bz(:,marg+1:2*marg,:)
    box%e(:,iy-marg+1:iy,:) = box[imx,behind,imz,1]%e(:,marg+1:2*marg,:)
    box%pr(:,iy-marg+1:iy,:) = box[imx,behind,imz,1]%pr(:,marg+1:2*marg,:)

    !z direction
    if (.not. down==0) then
        box%ro(:,:,1:marg) = box[imx,imy,down,1]%ro(:,:,iz-2*marg+1:iz-marg) 
        box%rovx(:,:,1:marg) = box[imx,imy,down,1]%rovx(:,:,iz-2*marg+1:iz-marg) 
        box%rovy(:,:,1:marg) = box[imx,imy,down,1]%rovy(:,:,iz-2*marg+1:iz-marg) 
        box%rovz(:,:,1:marg) = box[imx,imy,down,1]%rovz(:,:,iz-2*marg+1:iz-marg) 
        box%bx(:,:,1:marg) = box[imx,imy,down,1]%bx(:,:,iz-2*marg+1:iz-marg) 
        box%by(:,:,1:marg) = box[imx,imy,down,1]%by(:,:,iz-2*marg+1:iz-marg) 
        box%bz(:,:,1:marg) = box[imx,imy,down,1]%bz(:,:,iz-2*marg+1:iz-marg) 
        box%e(:,:,1:marg) = box[imx,imy,down,1]%e(:,:,iz-2*marg+1:iz-marg) 
        box%pr(:,:,1:marg) = box[imx,imy,down,1]%pr(:,:,iz-2*marg+1:iz-marg) 
    end if 

    if (.not. up==0) then
        box%ro(:,:,iz-marg+1:iz) = box[imx,imy,up,1]%ro(:,:,marg+1:2*marg)
        box%rovx(:,:,iz-marg+1:iz) = box[imx,imy,up,1]%rovx(:,:,marg+1:2*marg)
        box%rovy(:,:,iz-marg+1:iz) = box[imx,imy,up,1]%rovy(:,:,marg+1:2*marg)
        box%rovz(:,:,iz-marg+1:iz) = box[imx,imy,up,1]%rovz(:,:,marg+1:2*marg)
        box%bx(:,:,iz-marg+1:iz) = box[imx,imy,up,1]%bx(:,:,marg+1:2*marg)
        box%by(:,:,iz-marg+1:iz) = box[imx,imy,up,1]%by(:,:,marg+1:2*marg)
        box%bz(:,:,iz-marg+1:iz) = box[imx,imy,up,1]%bz(:,:,marg+1:2*marg)
        box%e(:,:,iz-marg+1:iz) = box[imx,imy,up,1]%e(:,:,marg+1:2*marg)
        box%pr(:,:,iz-marg+1:iz) = box[imx,imy,up,1]%pr(:,:,marg+1:2*marg)
    end if
    
end subroutine
    
end module
