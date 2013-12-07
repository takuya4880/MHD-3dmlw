program main
    use defstruct   
    use ic
    use bc
    use op 
    use pr
    use dt
    use st
    !$use omp_lib
    implicit none
    
    type(cell) :: box[cox,coz,*]
    integer :: i,j
    double precision :: uboundary(9,marg)
    double precision :: t, tint, tend, tnxt
    integer :: start(8), time(8), minits, timelimit
    integer :: flag(cox,coz), timeup[cox,coz,*] 
    integer :: ns, nsout 
    character*10 :: tmp
    integer :: mcont

    call omp_set_num_threads(1)
    !allocate(box)
    !open(23,file="result.dat",status="replace")

    call date_and_time(tmp,tmp,tmp,start)

    mcont = 0
    timelimit = 210 !: 210:3.5hours
    box%con%nx = nx
    box%con%nz = nz
    box%con%ix = ix
    box%con%iz = iz
    box%con%imx = this_image(box,1)
    box%con%imz = this_image(box,2)
    box%con%marg = marg
    box%con%wid = 150.
    box%con%hig = 60.
    box%con%dx = box%con%wid/dble(nnx-1)
    box%con%dz = box%con%hig/dble(nnz-1)
    box%con%a = 0.4
    box%con%q = 3.
    box%con%gam = 5./3.

    ns = 0
    nsout = 1e5
    t = 0.
    tint = 1.
    tnxt = tint
    tend = 80.

    call initial(box, uboundary)
    sync all
    call boundary(box, uboundary)
    call outpinit(box)
    if (mcont==1) then
        call readdata(box,t)
        tnxt = dint(t) + tint
    end if
    call outp(box,t)
    call pressure(box)

    do
        call detdt(box)    
        call step(box)
        sync all
        call boundary(box, uboundary)
        t = t + box%con%dt
        ns = ns + 1
        if (box%con%imx*box%con%imz==1) print *,t,box%con%dt 
        if (t>=tnxt .or. ns>=nsout) then
            call outp(box,t)
            tnxt = tnxt + tint
            ns = 0
        endif
        if (t>tend) exit
        if (box%con%dt<1.e-10) exit

        call date_and_time(tmp,tmp,tmp,time)
        time = time - start
        minits = time(3)*24*60+time(5)*60+time(6)
        timeup = minits/timelimit   
        sync all
        do i=1,cox
            do j=1,coz
                flag(i,j) = timeup[i,j,1]
            end do
        end do 
        if (product(flag)==1) exit
        
    end do

end program main
        
