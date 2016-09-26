pro allphysj

; find ang mom of ALL gas particles at the time they enter the halo.

tunit=1.223128407d18/3.1556926d7/1d9 ; Gyr
munit=1.84793d16
vunit = 1260.
lunit=50000.


allgas = mrdfits('grp1.allgas.iord.fits',0)
early = mrdfits('early.iord.fits',0)
accrz = mrdfits('grp1.accrz.fits',0)
zs = accrz[uniq(accrz,sort(accrz))]
zs = reverse(zs)
hubble = wmap3_hubble(zs)
a = 1./(1.+zs)


match,early,allgas,e1,e2
lremove,e2,allgas
nunique = n_elements(allgas)

; now need to read in all the files and align them.
haloidfile = file_search("*.haloid.dat")
readcol,haloidfile[0],filename,halos,format='a,i',/silent
filename=reverse(filename)
halos = reverse(halos)
nfiles = n_elements(filename)
filediff = nfiles-n_elements(zs)

jhalo = {iords:lonarr(nunique),jx:fltarr(nunique),jy:fltarr(nunique),jz:fltarr(nunique)}

for i=0,nfiles-filediff-1 do begin
    accretednow = where(accrz eq zs[i],nnow)
    print,nnow,' gas particles accreted '
    rtipsy,filename[i+filediff],h,g,d,s
    grp = read_lon_array(filename[i+filediff]+'.amiga.grp')
    iord = read_lon_array(filename[i+filediff]+'.iord')
    gas = g[where(grp[0:h.ngas-1] eq halos[i+filediff])]
    dark = d[where(grp[h.ngas:h.ngas+h.ndark-1] eq halos[i+filediff])]
    star = s[where(grp[h.ngas+h.ngas:h.n-1] eq halos[i+filediff])]
    iordhalo = iord[where(grp eq halos[i+filediff])]
    readcol,filename[i+filediff]+'.shrinkcenters',haloid,cx,cy,cz,/silent
    location = where(haloid eq halos[i+filediff])
    cx = cx[location]
    cy = cy[location]
    cz = cz[location]
    gas.x = gas.x-cx[0]
    dark.x = dark.x-cx[0]
    star.x = star.x-cx[0]
    gas.y = gas.y-cy[0]
    dark.y = dark.y-cy[0]
    star.y = star.y-cy[0]
    gas.z = gas.z-cz[0]
    dark.z = dark.z-cz[0]
    star.z = star.z-cz[0]

    ; align the particles to the ang mom of the gas in the disk
    align,gas,dark,star,3.0/lunit  ; within 3 kpc

   ; now calculate J for the gas particles I want to trace.
    ; match returns the first instance when a particle got accreted.
    match,iordhalo,allgas[accretednow],i1,i2
    accretedgas = gas[i1]
    if total(i1) gt 0 then print,n_elements(i1),' are unique' else print,'youve done it wrong again'

    ex = accretedgas.x*a[i]
    why = accretedgas.y*a[i]
    zee = accretedgas.z*a[i]
    ; alyson says:  v_phys = a*(H*x + xdot)
    hubble = hubble*2.894405/73.
    vex = a[i] * (hubble[i]*accretedgas.x + accretedgas.vx)
    vwhy = a[i] * (hubble[i]*accretedgas.y + accretedgas.vy)
    vzee = a[i] * (hubble[i]*accretedgas.z + accretedgas.vz)

    Jx = why*vzee - zee*vwhy
    Jy = zee*vex - ex*vzee
    Jz = ex*vwhy - why*vex

    ; put in array to save
    jhalo.iords[accretednow] = allgas[accretednow]
    jhalo.jx[accretednow] = jx
    jhalo.jy[accretednow] = jy
    jhalo.jz[accretednow] = jz


endfor
save,jhalo,filename='jhalophys.dat'

stop

end
