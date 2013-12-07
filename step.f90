module st
implicit none 
contains
subroutine step(box)
    use defstruct
    use lw
    use pr

    implicit none
    type(cell) :: box
    type(cell), pointer :: fx, fy, fz, s, h, d
    integer :: i
    allocate(fx, fy, fz, s, h, d)
    fx%con = box%con
    fy%con = box%con
    fz%con = box%con
    s%con = box%con
    h%con = box%con
    d%con = box%con
    h%con%a = -1.   !use a to knowtify this is half step value for source term
        
    h%x = box%x + 0.5*box%con%dx
    h%y = box%y + 0.5*box%con%dy
    h%z = box%z + 0.5*box%con%dz

    call flux(box, fx, fy, fz)
    call source(box, s)
    call lw1(box, h, d, fx, fy, fz, s)
    call pressure(h)
    
    call flux(h, fx, fy, fz)
    call source(h, s)
    call lw2(box, d, fx, fy, fz, s)

    call artvis(box, d)
    call pressure(box)

    deallocate(fx, fy, fz, s, h, d)

end subroutine
end module
