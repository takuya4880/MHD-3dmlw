program main
    use defstruct   
    use ic
    use bc
    use op 
    use pr
    use dt
    use st
    implicit none
    
    type(cell) :: box[cox,coy,coz,*]
    integer :: i,j,k
    double precision :: uboundary(9,marg)
    double precision :: t, tint, tend, tnxt
    integer :: start(8), time(8), minits, timelimit
    integer :: flag(cox,coy,coz), timeup[cox,coy,coz,*] 
    character*10 :: tmp
    integer :: mcont

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

    t = 0
    tint = 5.
    tnxt = tint
    tend = 120.

    call initial(box, uboundary)
    sync all
    call boundary(box, uboundary)
    sync all
    call outpinit_fl(box)
    if (mcont==1) then
        call readdata(box,t)
        if (t>79) tint=1.
        tnxt = t - mod(t,tint) + tint
    else
        call outp_fl(box,t)
    end if
    call pressure(box)

    do
        if (t>79) tint=1.
        call detdt(box)    
        sync all
        call step(box)
        sync all
        call boundary(box, uboundary)
        sync all
        t = t + box%con%dt
        if (t>=tnxt) then
            if (this_image()==1) print *,t,box%con%dt 
            call outp_fl(box,t)
            tnxt = tnxt + tint
        endif
        if (t>tend) exit

        sync all
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
        sync all
        if (product(flag)==1 .or. box%con%dt<1.e-10) then
            if (this_image()==1) print *,t,box%con%dt 
            call outpinit_db_in(box)
            call outp_db(box,t)
            exit
        end if
    end do

end program main
        
