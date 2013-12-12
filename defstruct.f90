module defstruct
    implicit none
    integer,parameter :: nx=100
    integer,parameter :: ny=70
    integer,parameter :: nz=100
    integer,parameter :: cox=7
    integer,parameter :: coy=3
    integer,parameter :: coz=6
    integer,parameter :: nnx=nx*cox
    integer,parameter :: nny=ny*coy
    integer,parameter :: nnz=nz*coz
    integer,parameter :: marg=4
    integer,parameter :: ix=nx+2*marg
    integer,parameter :: iy=ny+2*marg
    integer,parameter :: iz=nz+2*marg
    integer,parameter :: iix=nnx+2*marg
    integer,parameter :: iiy=nny+2*marg
    integer,parameter :: iiz=nnz+2*marg
    

    type constants 
        integer imx,imy,imz
        double precision dx, dy, dz, dt, wid, dep, hig
        double precision gam, q, a
        double precision gx, gy, gz 
    end type

    type output
        integer mf_params, mf_t, mf_ro, mf_pr, mf_vx, mf_vy
        integer mf_vz, mf_bx, mf_by, mf_bz, mf_x, mf_y, mf_z
        integer mfi_t, mfi_ro, mfi_pr, mfi_vx, mfi_vy
        integer mfi_vz, mfi_bx, mfi_by, mfi_bz
    end type
    
    type cell
        type(constants) con
        type(output) op
        double precision x(ix), y(iy), z(iz)
        double precision ro(ix,iy,iz)
        double precision rovx(ix,iy,iz), rovy(ix,iy,iz), rovz(ix,iy,iz)
        double precision bx(ix,iy,iz), by(ix,iy,iz), bz(ix,iy,iz)
        double precision pr(ix,iy,iz), e(ix,iy,iz), eta(ix,iy,iz) 
    end type    


    
end module

