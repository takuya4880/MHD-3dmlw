module ic
implicit none 
contains
subroutine initial(box, uboundary)
    use defstruct
    implicit none
    type(cell) :: box[cox,coy,coz,*]
    double precision :: uboundary(9,marg)

    integer :: i,j,m,origin
    integer :: offset
    double precision :: gami               !inberse of gamma
    double precision :: wid
    double precision :: amp, tpt, tpho, tcor, x, z, a, ad, lp, phicor
    double precision :: w, ztr, zfsl, zfsu, zmgc, zinv, betafs, betacor
    double precision :: den(iiz), pre(iiz),temp(iiz),phi(iiz)
    double precision :: beta(iiz),beta1i(iiz),beta2i(iiz),b(iiz),ee(iiz) 
    double precision :: zz(iiz)
    m = marg
    wid = box%con%wid
    gami = 1./box%con%gam

    amp = 0.05
    lp = 20. 
    phicor = 4.*atan(1.d0) 
    tpt = 25.
    tpho = 1.
    tcor = tpho * tpt
    w = 0.5
    ztr = 8.
    zfsl = -4.
    zfsu = -2.
    zmgc = 9.5
    zinv = 5.
    betafs=4.
    betacor=0.2
    a = 2.
    ad = a*(box%con%gam-1.)/box%con%gam

    origin = int(5./box%con%hig*nnz)+1+m
    
    forall(i=1:ix) box%x(i)=box%con%dx*(nx*(box%con%imx-1)+i-m)
    forall(i=1:iy) box%y(i)=box%con%dy*(ny*(box%con%imy-1)+i-m)
    forall(i=1:iz) box%z(i)=box%con%dz*(nz*(box%con%imz-1)+i-origin)
    forall(i=1:iiz) zz(i)=box%con%dz*(i-origin)
    
    box%con%gx = 0.
    box%con%gy = 0.
    box%con%gz = -gami 

    do i=origin+1,iiz
        temp(i) = tpho + 0.5 * (tcor-tpho)*(tanh((zz(i)-ztr)/w) + 1.)
    end do
    do i=origin,1,-1
        temp(i) = tpho - ad*zz(i)
    end do

    den(origin) = 1.
    pre(origin) = gami*temp(origin)
    
    beta1i = 0.25/betafs*(tanh((zz-zfsl)/w) + 1.) * (-tanh((zz-zfsu)/w) + 1.)
    beta2i = 0.5/betacor*(tanh((zz-zmgc)/w) + 1.)
    beta = 1./(beta1i+beta2i)
    do i=origin+1,iiz
        den(i) = den(i-1) * ((1.+1./beta(i-1))*temp(i-1) + 0.5*box%con%gam*box%con%dz*box%con%gz)&
                          / ((1.+1./beta(i))*temp(i) - 0.5*box%con%gam*box%con%dz*box%con%gz)
    end do
    do i=origin-1,1,-1
        den(i) = den(i+1) * ((1.+1./beta(i+1))*temp(i+1) - 0.5*box%con%gam*box%con%dz*box%con%gz)&
                          / ((1.+1./beta(i))*temp(i) + 0.5*box%con%gam*box%con%dz*box%con%gz)
    end do
    
    pre = pre(origin) * (den/den(origin)) * (temp/temp(origin))
    b = sqrt(2.*pre/beta)
    
    do i=1,iiz
        if (zz(i)<zinv) then
            phi(i)=0.
        else
            phi(i)=phicor
        end if
    end do

    !open(24,file="initial.dat",status="replace")
    !do i=1,iz
    !   write (24,*) box%z(i), den(i), pre(i), temp(i), 1./beta(i), b(i), phi(i)
    !end do
    !close(24)

    offset = nz*(box%con%imz-1) 

    forall(i=1:iz) box%ro(:,:,i) = den(i+offset)
    box%rovx = 0.
    box%rovy = 0.
    box%rovz = 0.
    do i=1,ix
        do j=1,iz
            x = box%x(i)
            z = box%z(j)
            if (z>zfsl .and. z<zfsu .and. x>0.5*wid-0.25*lp .and. x<0.5*wid+0.25*lp) then
                box%rovz(i,:,j) = box%ro(i,:,j)*amp*cos(8.*atan(1.d0)*(x-0.5*wid)/lp)
            end if
        end do
    end do

    forall(i=1:iz) box%bx(:,:,i) = b(i+offset)*cos(phi(i+offset))
    box%by = 0.
    forall(i=1:iz) box%bz(:,:,i) = b(i+offset)*sin(phi(i+offset))
    forall(i=1:iz) box%pr(:,:,i) = pre(i+offset) 
    box%e = 0.5*(box%rovx**2 + box%rovy**2 + box%rovz**2)/box%ro &
            + box%pr/(box%con%gam-1.) &
            + 0.5*(box%bx**2 + box%by**2 + box%bz**2)

    uboundary(1,1:marg) = den(iiz-marg+1:iiz) - den(iiz-marg)
    uboundary(2,1:marg) = box%rovx(10,10,iz-marg+1:iz) - box%rovx(10,10,iz-marg)
    uboundary(3,1:marg) = box%rovy(10,10,iz-marg+1:iz) - box%rovy(10,10,iz-marg)
    uboundary(4,1:marg) = box%rovz(10,10,iz-marg+1:iz) - box%rovz(10,10,iz-marg)
    uboundary(5,1:marg) = box%bx(10,10,iz-marg+1:iz) - box%bx(10,10,iz-marg)
    uboundary(6,1:marg) = box%by(10,10,iz-marg+1:iz) - box%by(10,10,iz-marg)
    uboundary(7,1:marg) = box%bz(10,10,iz-marg+1:iz) - box%bz(10,10,iz-marg)
    uboundary(8,1:marg) = box%e(10,10,iz-marg+1:iz) - box%e(10,10,iz-marg)
    uboundary(9,1:marg) = pre(iiz-marg+1:iiz) - pre(iiz-marg)

end subroutine
end module

