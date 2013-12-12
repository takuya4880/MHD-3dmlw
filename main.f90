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
    
    type(cell) :: box[cox,coy,coz,*]
    integer :: i,j,k
    double precision :: uboundary(9,marg)
    double precision :: t, tint, tend, tnxt
    integer :: start(8), time(8), minits, timelimit
    integer :: flag(cox,coy,coz), timeup[cox,coy,coz,*] 
    integer :: ns, nsout 
    character*10 :: tmp
    integer :: mcont

    call omp_set_num_threads(1)
    !allocate(box)

    call date_and_time(tmp,tmp,tmp,start)

    mcont = 0
    timelimit = 210 !: 210:3.5hours
    box%con%imx = this_image(box,1)
    box%con%imy = this_image(box,2)
    box%con%imz = this_image(box,3)
    box%con%wid = 160.
    box%con%dep = 50.
    box%con%hig = 85.
    box%con%dx = box%con%wid/dble(nnx-1)
    box%con%dy = box%con%dep/dble(nny-1)
    box%con%dz = box%con%hig/dble(nnz-1)
    box%con%a = 0.4
    box%con%q = 3.
    box%con%gam = 1.4

    t = 0.
    tint = 1.
    tnxt = tint
    tend = 80.
    ns = 0
    nsout = 1e5

    call initial(box, uboundary)
    sync all
    call boundary(box, uboundary)
    sync all
    call outpinit(box)
    if (mcont==1) then
        call readdata(box,t)
        tnxt = t + tint
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
        if (box%con%imx*box%con%imy*box%con%imz==1) print *,t,box%con%dt 
        if (t>=tnxt .or. ns>=nsout) then
            call outp(box,t)
            tnxt = tnxt + tint
            ns = 0
        endif

        call date_and_time(tmp,tmp,tmp,time)
        time = time - start
        minits = time(3)*24*60+time(5)*60+time(6)
        timeup = minits/timelimit   
        sync all
        do i=1,cox
            do j=1,coy
                do k=1,coz
                    flag(i,j,k) = timeup[i,j,k,1]
                end do
            end do
        end do 
        if (product(flag)==1 .or. t>tend .or. box%con%dt<1.e-10) exit
        
    end do

end program main
        
