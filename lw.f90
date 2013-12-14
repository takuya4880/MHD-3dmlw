module lw
implicit none 
contains
subroutine lw1(box, h, d, fx, fy, fz, s)
    use defstruct
    implicit none
    type(cell) :: box, h, d, fx, fy, fz, s
    
    double precision :: dx, dy, dz, dt
    dx = box%con%dx
    dy = box%con%dy
    dz = box%con%dz
    dt = box%con%dt

    call each1(box%ro,h%ro,d%ro,fx%ro,fy%ro,fz%ro,s%ro,dx,dy,dz,dt)
    call each1(box%rovx,h%rovx,d%rovx,fx%rovx,fy%rovx,fz%rovx,s%rovx,dx,dy,dz,dt)
    call each1(box%rovy,h%rovy,d%rovy,fx%rovy,fy%rovy,fz%rovy,s%rovy,dx,dy,dz,dt)
    call each1(box%rovz,h%rovz,d%rovz,fx%rovz,fy%rovz,fz%rovz,s%rovz,dx,dy,dz,dt)
    call each1(box%bx,h%bx,d%bx,fx%bx,fy%bx,fz%bx,s%bx,dx,dy,dz,dt)
    call each1(box%by,h%by,d%by,fx%by,fy%by,fz%by,s%by,dx,dy,dz,dt)
    call each1(box%bz,h%bz,d%bz,fx%bz,fy%bz,fz%bz,s%bz,dx,dy,dz,dt)
    call each1(box%e,h%e,d%e,fx%e,fy%e,fz%e,s%e,dx,dy,dz,dt)
    
end subroutine

subroutine each1(u,h,d,fx,fy,fz,s,dx,dy,dz,dt)
    use defstruct
    implicit none
    double precision :: u(ix,iy,iz),d(ix,iy,iz),h(ix,iy,iz)
    double precision :: fx(ix,iy,iz),fy(ix,iy,iz),fz(ix,iy,iz),s(ix,iy,iz)
    double precision :: dx, dy, dz, dt
    
    integer i,j,k
    double precision fffx, fffy, fffz, ss
    double precision ddx, ddy, ddz
    
    ddx = dt/dx
    ddy = dt/dy
    ddz = dt/dz
    
    do k=2,iz-1
     do j=2,iy-1
      do i=2,ix-1
          d(i,j,k) = -0.5*ddx*(0.5*(fx(i+1,j,k)-fx(i-1,j,k)))&
                     -0.5*ddy*(0.5*(fy(i,j+1,k)-fy(i,j-1,k)))&
                     -0.5*ddz*(0.5*(fz(i,j,k+1)-fz(i,j,k-1)))&
                     +0.5*dt*s(i,j,k)
      end do
     end do
    end do

    do k=1,iz-1
     do j=1,iy-1
      do i=1,ix-1
          fffx = 0.25* ( fx(i+1,j+1,k+1)+fx(i+1,j,k+1)-fx(i,j+1,k+1)-fx(i,j,k+1) &
                        +fx(i+1,j+1,k  )+fx(i+1,j,k  )-fx(i,j+1,k  )-fx(i,j,k  ) )
          fffy = 0.25* ( fy(i+1,j+1,k+1)-fy(i+1,j,k+1)+fy(i,j+1,k+1)-fy(i,j,k+1) &
                        +fy(i+1,j+1,k  )-fy(i+1,j,k  )+fy(i,j+1,k  )-fy(i,j,k  ) )
          fffz = 0.25* ( fz(i+1,j+1,k+1)+fz(i+1,j,k+1)+fz(i,j+1,k+1)+fz(i,j,k+1) &
                        -fz(i+1,j+1,k  )-fz(i+1,j,k  )-fz(i,j+1,k  )-fz(i,j,k  ) )
          ss   = 0.125* ( s(i+1,j+1,k+1)+ s(i+1,j,k+1)+ s(i,j+1,k+1)+ s(i,j,k+1) &
                        + s(i+1,j+1,k  )+ s(i+1,j,k  )+ s(i,j+1,k  )+ s(i,j,k  ) )
          h(i,j,k)=0.125*(u(i+1,j+1,k+1)+ u(i+1,j,k+1)+ u(i,j+1,k+1)+ u(i,j,k+1) &
                        + u(i+1,j+1,k  )+ u(i+1,j,k  )+ u(i,j+1,k  )+ u(i,j,k  ) )&
                        - (ddx*fffx + ddy*fffy + ddz*fffz - dt*ss)
      end do
     end do
    end do
end subroutine

subroutine lw2(box, d, fx, fy, fz, s)
    use defstruct
    implicit none
    type(cell) :: box, d, fx, fy, fz, s
    
    double precision :: dx, dy, dz, dt
    dx = box%con%dx
    dy = box%con%dy
    dz = box%con%dz
    dt = box%con%dt

    call each2(box%ro,d%ro,fx%ro,fy%ro,fz%ro,s%ro,dx,dy,dz,dt)
    call each2(box%rovx,d%rovx,fx%rovx,fy%rovx,fz%rovx,s%rovx,dx,dy,dz,dt)
    call each2(box%rovy,d%rovy,fx%rovy,fy%rovy,fz%rovy,s%rovy,dx,dy,dz,dt)
    call each2(box%rovz,d%rovz,fx%rovz,fy%rovz,fz%rovz,s%rovz,dx,dy,dz,dt)
    call each2(box%bx,d%bx,fx%bx,fy%bx,fz%bx,s%bx,dx,dy,dz,dt)
    call each2(box%by,d%by,fx%by,fy%by,fz%by,s%by,dx,dy,dz,dt)
    call each2(box%bz,d%bz,fx%bz,fy%bz,fz%bz,s%bz,dx,dy,dz,dt)
    call each2(box%e,d%e,fx%e,fy%e,fz%e,s%e,dx,dy,dz,dt)

end subroutine

subroutine each2(box,d,fx,fy,fz,s,dx,dy,dz,dt)
    use defstruct
    implicit none
    double precision :: box(ix,iy,iz),d(ix,iy,iz)
    double precision :: fx(ix,iy,iz),fy(ix,iy,iz),fz(ix,iy,iz),s(ix,iy,iz)
    double precision :: dx, dy, dz, dt
    
    integer i,j,k
    double precision fffx,fffy,fffz,ss
    double precision ddx,ddy,ddz

    ddx = dt/dx
    ddy = dt/dy
    ddz = dt/dz
    
    do k=1,iz-2
     do j=1,iy-2
      do i=1,ix-2
          fffx = 0.25* ( fx(i+1,j+1,k+1)+fx(i+1,j,k+1)-fx(i,j+1,k+1)-fx(i,j,k+1) &
                        +fx(i+1,j+1,k  )+fx(i+1,j,k  )-fx(i,j+1,k  )-fx(i,j,k  ) )
          fffy = 0.25* ( fy(i+1,j+1,k+1)-fy(i+1,j,k+1)+fy(i,j+1,k+1)-fy(i,j,k+1) &
                        +fy(i+1,j+1,k  )-fy(i+1,j,k  )+fy(i,j+1,k  )-fy(i,j,k  ) )
          fffz = 0.25* ( fz(i+1,j+1,k+1)+fz(i+1,j,k+1)+fz(i,j+1,k+1)+fz(i,j,k+1) &
                        -fz(i+1,j+1,k  )-fz(i+1,j,k  )-fz(i,j+1,k  )-fz(i,j,k  ) )
          ss   = 0.125* ( s(i+1,j+1,k+1)+ s(i+1,j,k+1)+ s(i,j+1,k+1)+ s(i,j,k+1) &
                        + s(i+1,j+1,k  )+ s(i+1,j,k  )+ s(i,j+1,k  )+ s(i,j,k  ) )
          d(i+1,j+1,k+1) = d(i+1,j+1,k+1) - 0.5*(ddx*fffx + ddy*fffy + ddz*fffz - dt*ss)
      end do
     end do
    end do

end subroutine 

subroutine artvis(box, d)
    use defstruct
    implicit none
    type(cell) :: box, d
    double precision, allocatable :: kapx(:,:,:),kapy(:,:,:),kapz(:,:,:)
    allocate(kapx(ix,iy,iz),kapy(ix,iy,iz),kapz(ix,iy,iz))
    
    kapx(2:ix,:,:) = box%con%q * abs(box%rovx(2:ix,:,:)/box%ro(2:ix,:,:) &
                                    - box%rovx(1:ix-1,:,:)/box%ro(1:ix-1,:,:))
    kapy(:,2:iy,:) = box%con%q * abs(box%rovy(:,2:iy,:)/box%ro(:,2:iy,:) &
                                    - box%rovy(:,1:iy-1,:)/box%ro(:,1:iy-1,:))
    kapz(:,:,2:iz) = box%con%q * abs(box%rovz(:,:,2:iz)/box%ro(:,:,2:iz) &
                                    - box%rovz(:,:,1:iz-1)/box%ro(:,:,1:iz-1)) 

    call eachav(box%ro, d%ro, kapx, kapy, kapz, box%con)
    call eachav(box%rovx, d%rovx, kapx, kapy, kapz, box%con)
    call eachav(box%rovy, d%rovy, kapx, kapy, kapz, box%con)
    call eachav(box%rovz, d%rovz, kapx, kapy, kapz, box%con)
    call eachav(box%bx, d%bx, kapx, kapy, kapz, box%con)
    call eachav(box%by, d%by, kapx, kapy, kapz, box%con)
    call eachav(box%bz, d%bz, kapx, kapy, kapz, box%con)
    call eachav(box%e, d%e, kapx, kapy, kapz, box%con)

    deallocate(kapx,kapy,kapz)

end subroutine

subroutine eachav(box,d,kapx,kapy,kapz,con)
    use defstruct
    double precision :: box(ix,iy,iz),d(ix,iy,iz)
    double precision :: kapx(ix,iy,iz),kapy(ix,iy,iz),kapz(ix,iy,iz)
    type(constants) con
    
    double precision :: ddx, ddy, ddz
    double precision, allocatable :: difx(:,:,:),dify(:,:,:),difz(:,:,:)
    allocate(difx(ix,iy,iz),dify(ix,iy,iz),difz(ix,iy,iz))
    ddx = con%dt / con%dx
    ddy = con%dt / con%dy
    ddz = con%dt / con%dz
    
    difx(2:ix,:,:) = box(2:ix,:,:) - box(1:ix-1,:,:)
    dify(:,2:iy,:) = box(:,2:iy,:) - box(:,1:iy-1,:)
    difz(:,:,2:iz) = box(:,:,2:iz) - box(:,:,1:iz-1)

    box(3:ix-2,3:iy-2,3:iz-2) = box(3:ix-2,3:iy-2,3:iz-2) + d(3:ix-2,3:iy-2,3:iz-2) &
             + ddx * ( kapx(4:ix-1,3:iy-2,3:iz-2)*difx(4:ix-1,3:iy-2,3:iz-2) &
                     - kapx(3:ix-2,3:iy-2,3:iz-2)*difx(3:ix-2,3:iy-2,3:iz-2) ) &
             + ddy * ( kapy(3:ix-2,4:iy-1,3:iz-2)*dify(3:ix-2,4:iy-1,3:iz-2) &
                     - kapy(3:ix-2,3:iy-2,3:iz-2)*dify(3:ix-2,3:iy-2,3:iz-2) ) &
             + ddz * ( kapz(3:ix-2,3:iy-2,4:iz-1)*difz(3:ix-2,3:iy-2,4:iz-1) &
                     - kapz(3:ix-2,3:iy-2,3:iz-2)*difz(3:ix-2,3:iy-2,3:iz-2) ) 
    deallocate(difx,dify,difz)
end subroutine 

subroutine flux(box, fx, fy, fz)
    use defstruct
    implicit none
    type(cell) :: box, fx, fy, fz
    
    integer :: i,j,k
    double precision :: alp=0.01, etamax=1., vc=1000
    double precision :: b2,roi,h,vx,vy,vz,bx,by,bz,pr
    double precision :: eta, ex, ey, ez
    double precision :: jx, jy, jz

    do k=2,iz-1
      do j=2,iy-1
        do i=2,ix-1
            roi = 1./box%ro(i,j,k)
            vx = box%rovx(i,j,k) * roi
            vy = box%rovy(i,j,k) * roi
            vz = box%rovz(i,j,k) * roi
            bx = box%bx(i,j,k)
            by = box%by(i,j,k)
            bz = box%bz(i,j,k)
            pr = box%pr(i,j,k)
            b2 = bx**2 + by**2 + bz**2
            h = 0.5*(vx**2+vy**2+vz**2)*box%ro(i,j,k) &
                        + pr*box%con%gam/(box%con%gam-1.)

            jx = (box%bz(i,j+1,k)-box%bz(i,j-1,k))/(2.*box%con%dy) &
                        -(box%by(i,j,k+1)-box%by(i,j,k-1))/(2.*box%con%dz)
            jy = (box%bx(i,j,k+1)-box%bx(i,j,k-1))/(2.*box%con%dz) &
                        -(box%bz(i+1,j,k)-box%bz(i-1,j,k))/(2.*box%con%dx)
            jz = (box%by(i+1,j,k)-box%by(i-1,j,k))/(2.*box%con%dx) &
                        -(box%bx(i,j+1,k)-box%bx(i,j-1,k))/(2.*box%con%dy)

            eta = sqrt((jx**2+jy**2+jz**2)*16.*atan(1.0))*roi  !calculate vd
            
            if (eta<vc) then
                eta = 0.
            else 
                eta = alp*(eta/vc-1.)**2
                if (eta>etamax) then
                    eta = etamax
                end if
            end if
            ex = eta*jx + (-vy*bz+vz*by)
            ey = eta*jy + (-vz*bx+vx*bz)
            ez = eta*jz + (-vx*by+vy*bx)

            fx%ro(i,j,k) = box%rovx(i,j,k)
            fx%rovx(i,j,k) = vx*box%rovx(i,j,k) - bx*bx + pr + 0.5*b2
            fx%rovy(i,j,k) = vx*box%rovy(i,j,k) - bx*by
            fx%rovz(i,j,k) = vx*box%rovz(i,j,k) - bx*bz
            fx%bx(i,j,k) = 0.
            fx%by(i,j,k) = -ez
            fx%bz(i,j,k) = ey
            fx%e(i,j,k) = h*vx + (ey*bz - ez*by) 

            fy%ro(i,j,k) = box%rovy(i,j,k)
            fy%rovx(i,j,k) = vy*box%rovx(i,j,k) - by*bx
            fy%rovy(i,j,k) = vy*box%rovy(i,j,k) - by*by + pr + 0.5*b2
            fy%rovz(i,j,k) = vy*box%rovz(i,j,k) - by*bz
            fy%bx(i,j,k) = ez
            fy%by(i,j,k) = 0.
            fy%bz(i,j,k) = -ex
            fy%e(i,j,k) = h*vy + (ez*bx - ex*bz) 

            fz%ro(i,j,k) = box%rovz(i,j,k)
            fz%rovx(i,j,k) = vz*box%rovx(i,j,k) - bz*bx
            fz%rovy(i,j,k) = vz*box%rovy(i,j,k) - bz*by
            fz%rovz(i,j,k) = vz*box%rovz(i,j,k) - bz*bz + pr + 0.5*b2
            fz%bx(i,j,k) = -ey
            fz%by(i,j,k) = ex
            fz%bz(i,j,k) = 0.
            fz%e(i,j,k) = h*vz + (ex*by - ey*bx)
        end do
      end do
    end do
    

end subroutine

subroutine source(box, s)
    use defstruct
    implicit none
    type(cell) :: box, s
    integer :: i
    double precision :: fugou(iz)
   
    fugou = 1.
    if (box%con%imz==1) then
        if (box%con%a==-1) then 
            fugou(1:marg-1) = -1.
            fugou(marg) = 0.
        else
            fugou(1:marg) = -1.
        end if
    end if

    s%ro = 0.
    s%bx = 0.
    s%by = 0.
    s%bz = 0.
    s%rovx = box%ro*box%con%gx
    s%rovy = box%ro*box%con%gy
    forall(i=1:iz) s%rovz(:,:,i) = box%ro(:,:,i)*box%con%gz*fugou(i)
    s%e = box%rovx*box%con%gx + box%rovy*box%con%gy + box%rovz*box%con%gz

end subroutine

end module 
