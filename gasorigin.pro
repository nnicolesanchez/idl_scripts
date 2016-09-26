pro gasorigin,plotter=plotter

; in spherical shells around a central (most massive) BH
; find the mass of gas as it's broken into
; clumpy, cold, and shocked origins


; units
lunit = 50000.
vunit = 1260.
munit=1.84793e16
lunit=50000.
tunit=1.223e18/3.15d7/1d9 ; Gyr

file = file_search('*.00???')

; files
rtipsy,file[0],h,g,d,s
grp = read_lon_array(file[0]+'.amiga.grp')
iord = read_lon_array(file[0]+'.iord')
clumpy = mrdfits('../clumpy.accr.iord.fits',0)
shocked = mrdfits('../shocked.iord.fits',0)
cold = mrdfits('../unshock.iord.fits',0)
;early = mrdfits('../early.iord.fits',0)
;smooth = mrdfits('../smooth.accr.iord.fits',0)
;tracestep = mrdfits('../tracedtostep.fits',0)
;accrz = mrdfits('../grp1.accrz.fits',0)
;allgas = mrdfits('../grp1.allgas.iord.fits',0)


; cut out halo 1
g = g[where(grp[0:h.ngas-1] eq 1)]
d = d[where(grp[h.ngas:h.ngas+h.ndark-1] eq 1)]
s = s[where(grp[h.ngas+h.ndark:h.n-1] eq 1)]
afac = h.time

bh=where(s.tform lt 0.0,nbh)
biggestbh = max(s[bh].mass,wbig); central BH
x = s[bh[wbig]].x
y = s[bh[wbig]].y
z = s[bh[wbig]].z

; center galaxy on the BH
g.x = g.x-x
g.y = g.y-y
g.z = g.z-z
radius = sqrt(g.x^2. + g.y^2. + g.z^2.)*lunit*afac ; physical radius

; annuli out to 100 kpc
radii = findgen(78)/2. ; kpc
nrad = n_elements(radii)

coldmass=fltarr(nrad)
shockedmass= fltarr(nrad)
clumpymass = fltarr(nrad)

for i=0,nrad-2 do begin
    thisshell = where(radius ge radii[i] AND radius lt radii[i+1])
    thisiord = iord[thisshell]
    thisgas = g[thisshell]
    ; which is clumpy
    match,thisiord,clumpy,cl1,cl2
    ; which is cold
    match,thisiord,cold,c1,c2
    ; which is shocked
    match,thisiord,shocked,s1,s2
    ; determine masses of gas in each annulus
    clumpymass[i] = total(g[thisshell[cl1]].mass)*munit
    coldmass[i]    = total(g[thisshell[c1]].mass)*munit
    shockedmass[i] = total(g[thisshell[s1]].mass)*munit
endfor

; which is biggest?  use to set scale
maxscale = max(clumpymass) > max(coldmass) > max(shockedmass)

if keyword_set(plotter) then paperplot,filename=file[0]+'.gasorigin.ps'

plot,radii,clumpymass,xtit="Radius (physical kpc)",ytit='Gas Mass (M'+sunsymbol()+')',title="z = "+trim(1./afac-1.),yra=[0,maxscale];,/ylog,yra=[1e6,1e9]
  oplot,radii,coldmass,color=80
  oplot,radii,shockedmass,color=240
  oplot,radii,clumpymass,color=125

legend,['clumpy','cold','shocked'],lines=[0,0,0],colors=[125,80,240],/right,charsize=1

if keyword_set(plotter) then paperplot,/close

stop

; gas mass fraction vs radius

totalmass=clumpymass+coldmass+shockedmass

plot,radii,coldmass/totalmass,yra=[0,1],xtit="Radius (physical kpc)",ytit='Gas Mass Fraction',title="z = "+trim(1./afac-1.)
  oplot,radii,coldmass/totalmass,color=80
 oplot,radii,shockedmass/totalmass,color=240
  oplot,radii,clumpymass/totalmass,color=125

stop

end
