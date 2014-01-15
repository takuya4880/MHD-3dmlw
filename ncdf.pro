dims=size(ro)
ix=dims[1]
iy=dims[2]
iz=dims[3]

id = ncdf_create('ro.nc', /clobber)
ncdf_control, id, /fill

xid = ncdf_dimdef(id, 'x', ix)
yid = ncdf_dimdef(id, 'y', iy)
zid = ncdf_dimdef(id, 'z', iz)
roid = ncdf_vardef(id, 'ro', [xid,yid,zid], /float)
prid = ncdf_vardef(id, 'pr', [xid,yid,zid], /float)
bxid = ncdf_vardef(id, 'bx', [xid,yid,zid], /float)
byid = ncdf_vardef(id, 'by', [xid,yid,zid], /float)
bzid = ncdf_vardef(id, 'bz', [xid,yid,zid], /float)
vxid = ncdf_vardef(id, 'vx', [xid,yid,zid], /float)
vyid = ncdf_vardef(id, 'vy', [xid,yid,zid], /float)
vzid = ncdf_vardef(id, 'vz', [xid,yid,zid], /float)
teid = ncdf_vardef(id, 'te', [xid,yid,zid], /float)
cxid = ncdf_vardef(id, 'x', [xid], /float)
cyid = ncdf_vardef(id, 'y', [yid], /float)
czid = ncdf_vardef(id, 'z', [zid], /float)

ncdf_control, id, /endef

ncdf_varput, id, roid, ro[*,*,*,25]
ncdf_varput, id, prid, pr[*,*,*,25]
ncdf_varput, id, bxid, bx[*,*,*,25]
ncdf_varput, id, byid, by[*,*,*,25]
ncdf_varput, id, bzid, bz[*,*,*,25]
ncdf_varput, id, vxid, vx[*,*,*,25]
ncdf_varput, id, vyid, vy[*,*,*,25]
ncdf_varput, id, vzid, vz[*,*,*,25]
ncdf_varput, id, teid, te[*,*,*,25]
ncdf_varput, id, cxid, x
ncdf_varput, id, cyid, y
ncdf_varput, id, czid, z

ncdf_close, id

end
