dims=size(ro)
ix=dims[1]
iy=dims[2]
iz=dims[3]

id = ncdf_create('ro.nc', /clobber)
ncdf_control, id, /fill

xid = ncdf_dimdef(id, 'x', ix)
yid = ncdf_dimdef(id, 'y', iy)
zid = ncdf_dimdef(id, 'z', iz)
tid = ncdf_dimdef(id, 't', /unlimited)
vid = ncdf_vardef(id, 'ro', [xid,yid,zid,tid], /float)
vxid = ncdf_vardef(id, 'x', [xid], /float)
vyid = ncdf_vardef(id, 'y', [yid], /float)
vzid = ncdf_vardef(id, 'z', [zid], /float)

ncdf_control, id, /endef

ncdf_varput, id, vid, alog10(float(ro[*,*,*,0:-2]))
ncdf_varput, id, vxid, x
ncdf_varput, id, vyid, y
ncdf_varput, id, vzid, z

ncdf_close, id

end
