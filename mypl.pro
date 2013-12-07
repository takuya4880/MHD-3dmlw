set_plot,'z'

xsize=750
ysize=300
timeoffset=0
device,set_resolution=[xsize,ysize]
device,set_character_size=[6,9]

!p.charsize=2.0
!p.background=255
!p.color=0
!p.position=[0.1,0.1,0.9,0.9]
;!x.margin=[5,2]
;!y.margin=[3,2]
ct = 39
loadct, ct

mout = 0
gam = 5./3.

mx=7
namear=strarr(mx)
dminar=fltarr(mx)
dmaxar=fltarr(mx)
dlogar=intarr(mx)

m=0 & namear[m]='log(ro)' & dminar[m]=-9.00d0 & dmaxar[m]=1.00d0 & dlogar[m]=1
m=1 & namear[m]='vx'      & dminar[m]=-1.00d1 & dmaxar[m]=1.00d1
m=2 & namear[m]='vy'      & dminar[m]=-3.00d0 & dmaxar[m]=3.00d0
m=3 & namear[m]='log(pr)' & dminar[m]=-7.00d0 & dmaxar[m]=1.00d0 & dlogar[m]=1
m=4 & namear[m]='bx'      & dminar[m]=-4.40d0 & dmaxar[m]=4.40d0
m=5 & namear[m]='by'      & dminar[m]=-1.80d0 & dmaxar[m]=1.80d0
m=6 & namear[m]='T'       & dminar[m]=-2.00d2 & dmaxar[m]=5.00d2

m=mout
nlevels=128
dmin=dminar(m) & dmax=dmaxar(m)
dlevel=(dmax-dmin)/nlevels
levels1=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels-1)

for time=0,n_elements(t)-1 do begin

nlevels2=64
dmax=max(az[*,*,time]) & dmin=min(az[*,*,time])
dlevel=(dmax-dmin)/nlevels2
levels2=(dmin+dlevel/2.d0)+dlevel*findgen(nlevels2-1)

time_st=strmid(strcompress(string(t[time]),/remove_all),0,4)

xrange=[min(x),max(x)] & yrange=[min(y),max(y)]
cbrange=[dminar[m],dmaxar[m]]

case m of
    0: contour, alog10(ro[*,*,time]),x,y,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,/fill,levels=levels1
    1: contour, vx[*,*,time],x,y,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,/fill,levels=levels1
    2: contour, vy[*,*,time],x,y,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,/fill,levels=levels1
    3: contour, alog10(pr[*,*,time]),x,y,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,/fill,levels=levels1
    4: contour, bx[*,*,time],x,y,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,/fill,levels=levels1
    5: contour, by[*,*,time],x,y,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,/fill,levels=levels1
    6: contour, (gam*pr[*,*,time]/ro[*,*,time]),x,y,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,/fill,levels=levels1
endcase

contour, az[*,*,time],x,y,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,/noerase,levels=levels2
xyouts,0.45,0.93,namear[m],/normal,charsize=2
xyouts,0.8,0.93,'t='+time_st,/normal,charsize=2
color_bar, cbrange,0.92,0.1,0.95,0.9, ct=ct,charsize=1,/vertical

;filepng='tmp.png'
filepng='png'+'/'+string(time+timeoffset,format='(i3.3)')+'.png'
img=tvrd()
tvlct,red,green,blue,/get
write_png,filepng,img,red,green,blue

endfor

end
