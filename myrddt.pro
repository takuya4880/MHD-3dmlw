narg=-1
;narg=[15]
;narg=indgen(7)*10
dir='./in1/'
sq4pi=sqrt(16.*atan(1))

mpest='.'+string(0,form='(i4.4)')
dacgetparam,dir+'params.txt'+mpest   ,'mpex',mpex
dacgetparam,dir+'params.txt'+mpest   ,'margin',margin
dacgetparam,dir+'params.txt'+mpest   ,'gm',gm
dacget0s,dir+'t.dac'+mpest,t,narg=narg
t=float(t)
nx=n_elements(t)

mpestar='.'+string(indgen(mpex),form='(i4.4)')

ipear=intarr(mpex)
jpear=intarr(mpex)
kpear=intarr(mpex)
for mpe=0,mpex-1 do begin
  mpest=mpestar[mpe]
  dacgetparam,dir+'params.txt'+mpest   ,'ipe',ipe
  dacgetparam,dir+'params.txt'+mpest   ,'jpe',jpe
  dacgetparam,dir+'params.txt'+mpest   ,'kpe',kpe
  ipear[mpe]=ipe
  jpear[mpe]=jpe
  kpear[mpe]=kpe
endfor

files=dir+'x.dac'+mpestar[uniq(ipear,sort(ipear))]
dacget1d,files,x,margin=margin
files=dir+'y.dac'+mpestar[uniq(jpear,sort(jpear))]
dacget1d,files,y,margin=margin
files=dir+'z.dac'+mpestar[uniq(kpear,sort(kpear))]
dacget1d,files,z,margin=margin
ix=n_elements(x)
jx=n_elements(y)
kx=n_elements(z)

if 1 then begin
files=dir+'ro.dac'+mpestar
dacget4s,files,ro,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
ro=float(ro)
files=dir+'pr.dac'+mpestar
dacget4s,files,pr,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
pr=float(pr)
files=dir+'vx.dac'+mpestar
dacget4s,files,vx,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
vx=float(vx)
files=dir+'vy.dac'+mpestar
dacget4s,files,vy,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
vy=float(vy)
files=dir+'vz.dac'+mpestar
dacget4s,files,vz,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
vz=float(vz)
files=dir+'bx.dac'+mpestar
dacget4s,files,bx,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
bx=float(bx)*sq4pi
files=dir+'by.dac'+mpestar
dacget4s,files,by,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
by=float(by)*sq4pi
files=dir+'bz.dac'+mpestar
dacget4s,files,bz,narg=narg,margin=margin,ipear=ipear,jpear=jpear,kpear=kpear
bz=float(bz)*sq4pi
endif

te=pr/ro*float(gm)
x=float(x)
y=float(y)
z=float(z)
dx=x[1]-x[0]
dy=y[1]-y[0]
dz=z[1]-z[0]

cden=fltarr(ix,jx,kx,nx)
for n=0,nx-1 do begin
    for k=1,kx-2 do begin
        for j=1,jx-2 do begin
            for i=1,ix-2 do begin
                cx = (bz(i,j+1,k,n)-bz(i,j-1,k,n))/(2.*dy) - (by(i,j,k+1,n)-by(i,j,k-1,n))/(2.*dz)
                cy = (bx(i,j,k+1,n)-bx(i,j,k-1,n))/(2.*dz) - (bz(i+1,j,k,n)-bz(i-1,j,k,n))/(2.*dx)
                cz = (by(i+1,j,k,n)-by(i-1,j,k,n))/(2.*dx) - (bx(i,j+1,k,n)-bx(i,j-1,k,n))/(2.*dy)
                cden[i,j,k,n] = sqrt(cx^2+cy^2+cz^2)
            endfor
        endfor
    endfor
endfor

divb=fltarr(ix,jx,kx,nx)
for n=0,nx-1 do begin
    for k=1,kx-2 do begin
        for j=1,jx-2 do begin
            for i=1,ix-2 do begin
                dbx = (bx(i+1,j,k,n)-bx(i-1,j,k,n))/(2.*dx) 
                dby = (by(i,j+1,k,n)-by(i,j-1,k,n))/(2.*dy)
                dbz = (bz(i,j,k+1,n)-bz(i,j,k-1,n))/(2.*dz)
                divb[i,j,k,n] = dbx+dby+dbz
            endfor
        endfor
    endfor
endfor



end
