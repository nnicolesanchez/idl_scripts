pro tracej

tunit=1.223128407d18/3.1556926d7/1d9 ; Gyr
munit=1.84793d16
vunit = 1260.
lunit=50000.

haloidfile = file_search("../*.haloid.dat")
readcol,'../files.list',file,format='a',/silent
readcol,'../times.list',filetime,filez,format='f',/silent
readcol,haloidfile[0],filename,halos,format='a,i',/silent
filename=reverse(filename)
halos = reverse(halos)
nfiles=n_elements(file)
filediff = nfiles-n_elements(halos)
restore,'allgas.sav'
ngas = n_elements(allgas.zaccr)
uniquezs = allgas.zaccr[uniq(allgas.zaccr,sort(allgas.zaccr))]
uniquezs = reverse(uniquezs)
hubble = wmap3_hubble(uniquezs)
a = 1./(1.+uniquezs)
diff = n_elements(halos)-n_elements(uniquezs)
; need to make sure files and z's are matching up.
;deltafile = n_elements(filez) - n_elements(halos)
;goodz = filez[deltafile:n_elements(filetime)-1]  ; now zfile has the same indices as filename
; looking for a direct match in zs hasn't been working
; perhaps find the beginning and end of the indices instead?
;match,goodz,uniquezs,zm1,zm2
; but "early" is designated at z=50, which won't match up.  add it in.
;beginind = min(zm1)-1
;endind = max(zm1)
;filestouse = filename[beginind:endind]

; the last files aren't counted because nothing enters that the BH accretes.
angmom = {iords:allgas.iords,jx:fltarr(ngas),jy:fltarr(ngas),jz:fltarr(ngas),z:fltarr(ngas)}
stop
; still off by one.  z = 13 goes with step 12.  z = 15 is step 10.
for i=0,nfiles-diff-filediff-1 do begin
    ; which particles to trace?
    ; the ones who accreted at this timestep.
    thisgas = where(allgas.zaccr eq uniquezs[i],nthisgas)
    if nthisgas eq 0 then print,'no accreted gas at redshift',uniquezs[i]
    if nthisgas eq 0 then continue
    gasiords = allgas.iords[thisgas]
    ; which halo is it entering?  (not halo 1 initially!)
    thishalo = halos[i+diff-1]
    print,'halo file '+filename[i+diff-1]
    ; select particles from that halo.
    print,'reading file '+file[i+diff+filediff-1]
    print,'z = ',uniquezs[i]
    rtipsy,'../'+file[i+diff+filediff-1],h,g,d,s
    grp = read_lon_array('../'+file[i+diff+filediff-1]+'.amiga.grp')
    iord = read_lon_array('../'+file[i+diff+filediff-1]+'.iord')
    gas = g[where(grp[0:h.ngas-1] eq thishalo)]
    dark = d[where(grp[h.ngas:h.ngas+h.ndark-1] eq thishalo)]
    star = s[where(grp[h.ngas+h.ngas:h.n-1] eq thishalo)]
    ; need to center it before aligning.
    readcol,'../'+file[i+diff+filediff-1]+'.shrinkcenters',haloid,cx,cy,cz,/silent
    location = where(haloid eq thishalo)
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
    align,gas,dark,star,3.0/lunit  ; within 5 kpc

    ; now calculate J for the gas particles I want to trace.
    thishalogasiords = iord[where(grp[0:h.ngas-1] eq thishalo)]
    match,thishalogasiords,gasiords,g1,g2

    if n_elements(g1) ne nthisgas then begin
        ;stop
        print,'bad! not all accreted particles are in halo'
        print,'recovering ',n_elements(g1),' out of',nthisgas
        ;stop
        ; include particles whose halo inds are 0
        zerogasiords = iord[where(grp[0:h.ngas-1] eq 0)]
        match,zerogasiords,gasiords,z1,z2
        if total(z1) gt -0.5 then begin
            print,'include',n_elements(z1),' from halo 0'
            if total(g1) ge -0.5 then g1=[g1,z1] else g1=z1
            g1=g1[sort(g1)]
        ;stop
        ; now fix g2, which won't match g1 anymore
            match,iord[where(grp[0:h.ngas-1] eq thishalo OR grp[0:h.ngas-1] eq 0) ],gasiords,g11,g22
            g2=g22
        endif
    endif

    ; do all in physical coordinates.
    accretedgas = gas[g1] ; calculate J for this gas.
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
    angmom.jx[thisgas[g2]] = jx
    angmom.jy[thisgas[g2]] = jy
    angmom.jz[thisgas[g2]] = jz
    angmom.z[thisgas[g2]] = uniquezs[i]
;stop

endfor

haloentryj = angmom
save,haloentryj,filename='haloentryj.dat'



end
