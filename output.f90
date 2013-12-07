module op
implicit none
contains
subroutine outp(box,t)
    use defstruct
    use cansio
    implicit none
    type(cell) :: box
    double precision :: t

    write(box%op%mf_t) t
    write(box%op%mf_ro) box%ro
    write(box%op%mf_pr) box%pr
    write(box%op%mf_vx) box%rovx/box%ro
    write(box%op%mf_vy) box%rovz/box%ro
    write(box%op%mf_bx) box%bx
    write(box%op%mf_by) box%bz
    write(box%op%mf_az) box%bpot

end subroutine 

subroutine outpinit(box)
    use defstruct
    use cansio
    implicit none
    type (cell) :: box
    
    integer :: mpe
    character :: cno*4
    mpe = (box%con%imz-1)*cox + box%con%imx -1
    write(cno,'(i4.4)') mpe

    box%op%mf_params=9
    box%op%mf_t=10
    box%op%mf_x=11
    box%op%mf_y=12
    box%op%mf_ro=20
    box%op%mf_pr=21
    box%op%mf_vx=22
    box%op%mf_vy=23
    box%op%mf_bx=25
    box%op%mf_by=26
    box%op%mf_az=28
    
    call dacdefparam(box%op%mf_params,'params.txt.'//cno)
    call dacdef0s(box%op%mf_t,'t.dac.'//cno,6)
    call dacdef2s(box%op%mf_ro,'ro.dac.'//cno,6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_pr,'pr.dac.'//cno,6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_vx,'vx.dac.'//cno,6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_vy,'vy.dac.'//cno,6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_bx,'bx.dac.'//cno,6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_by,'by.dac.'//cno,6,box%con%ix,box%con%iz)
    call dacdef2s(box%op%mf_az,'az.dac.'//cno,6,box%con%ix,box%con%iz)

    call dacputparami(box%op%mf_params,'ix',box%con%ix)
    call dacputparami(box%op%mf_params,'jx',box%con%iz)
    call dacputparami(box%op%mf_params,'margin',box%con%marg)
    call dacputparami(box%op%mf_params,'mpi',1)
    call dacputparami(box%op%mf_params,'mpex',cox*coz)
    call dacputparami(box%op%mf_params,'mpe',mpe)
    call dacputparami(box%op%mf_params,'ipex',cox)
    call dacputparami(box%op%mf_params,'ipe',box%con%imx-1)
    call dacputparami(box%op%mf_params,'jpex',coz)
    call dacputparami(box%op%mf_params,'jpe',box%con%imz-1)
    

    call dacdef1d(box%op%mf_x,'x.dac.'//cno,6,box%con%ix)
    write(box%op%mf_x) box%x
    call dacdef1d(box%op%mf_y,'y.dac.'//cno,6,box%con%iz)
    write(box%op%mf_y) box%z
    call dacputparamd(box%op%mf_params,'gm',box%con%gam)

    close(box%op%mf_x)
    close(box%op%mf_y)
    close(box%op%mf_params)
     
end subroutine

subroutine readdata(box,t)
    use defstruct
    use cansio
    implicit none
    type (cell) :: box
    double precision :: t

    integer :: ndi,n,nx0,ix0,jx0,mtype
    integer :: mpe
    character :: cno*4
    mpe = (box%con%imz-1)*cox + box%con%imx -1
    write(cno,'(i4.4)') mpe

    ndi=1000

    box%op%mfi_t=60
    box%op%mfi_ro=70
    box%op%mfi_pr=71
    box%op%mfi_vx=72
    box%op%mfi_vy=73
    box%op%mfi_bx=74
    box%op%mfi_by=75
    box%op%mfi_az=76
    
    call dacopnr0s(box%op%mfi_t,'in/t.dac.'//cno,mtype,nx0)
    call dacopnr2s(box%op%mfi_ro,'in/ro.dac.'//cno,mtype,ix0,jx0,nx0)
    call dacopnr2s(box%op%mfi_pr,'in/pr.dac.'//cno,mtype,ix0,jx0,nx0)
    call dacopnr2s(box%op%mfi_vx,'in/vx.dac.'//cno,mtype,ix0,jx0,nx0)
    call dacopnr2s(box%op%mfi_vy,'in/vy.dac.'//cno,mtype,ix0,jx0,nx0)
    call dacopnr2s(box%op%mfi_bx,'in/bx.dac.'//cno,mtype,ix0,jx0,nx0)
    call dacopnr2s(box%op%mfi_by,'in/by.dac.'//cno,mtype,ix0,jx0,nx0)
    call dacopnr2s(box%op%mfi_az,'in/az.dac.'//cno,mtype,ix0,jx0,nx0)

    do n=1,ndi
        read(box%op%mfi_t,end=9900) t
        read(box%op%mfi_ro) box%ro
        read(box%op%mfi_pr) box%pr
        read(box%op%mfi_vx) box%rovx
        read(box%op%mfi_vy) box%rovz
        read(box%op%mfi_bx) box%bx
        read(box%op%mfi_by) box%bz
        read(box%op%mfi_az) box%bpot
    end do
9900  continue
    
    box%rovx = box%rovx*box%ro
    box%rovz = box%rovz*box%ro
    box%e = 0.5*(box%rovx**2 + box%rovy**2 + box%rovz**2)/box%ro &
            + box%pr/(box%con%gam-1.) &
            + 0.5*(box%bx**2 + box%by**2 + box%bz**2)

end subroutine

end module 
