;# just input your favourite galaxy: everything else should be aok!
;# "MW1.1024g1bwK"," ", " "
;#
pro accrmode, haloidoutfile, lunit, wmap1=wmap1, skip=skip

;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; This code creates a structure that finds the maximum temperature of a 
; gas particle in its history.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if keyword_set(wmap1) then densunit = 135.98 else densunit = 147.8344

readcol, haloidoutfile, files, halo, format='a,l'
;Output of haloid runs from low z to high z, but find_shocked requires high to low.  Fix here.
files = reverse(files)
halo = reverse(halo)
list = files+'.allgas.history.fits'
gtp_file = files+'.amiga.gtp'
statfile = files+'.amiga.stat'
iord_files = files+'.iord'

if keyword_set(skip) then goto, jump1
; ****************************************************
; define some stuff
;*******************************************************

alliords = mrdfits('grp1.allgas.iord.fits',0)
nall = n_elements(alliords)
nsteps = n_elements(halo)
print,nsteps
readcol, statfile[0], grp, mvir, rvir0, format='l,x,x,x,x,f,f', /silent
grpind = where(grp eq halo[0])
print,'Memory I am using',memory(/high)
allstruct = replicate({iord:0L, mass:LONARR(nsteps), grp:LONARR(nsteps), temp:FLTARR(nsteps), rho:FLTARR(nsteps), entropy:FLTARR(nsteps), radius:FLTARR(nsteps), vx:FLTARR(nsteps), vy:FLTARR(nsteps), vz:FLTARR(nsteps), x:FLTARR(nsteps), y:FLTARR(nsteps), z:FLTARR(nsteps)}, nall)
print,'Memory I am using',memory(/high)
allstruct.iord = alliords
history = mrdfits(list[0],1)
rtipsy, gtp_file[0], h,g,d,s
allstruct.temp[0] = history.temp
allstruct.mass[0] = history.mass
;allstruct.rho[0] = history.rho*densunit/(h.time^3.)
allstruct.rho[0] = history.rho*densunit
;allstruct.entropy[0] = alog10(history.temp^1.5/(history.rho*densunit/h.time^3.))
allstruct.entropy[0] = alog10(history.temp^1.5/(history.rho*densunit))
allstruct.grp[0] = history.haloid
;allstruct.radius[0] = ((h.time*(history.x-lunit*s[0].x))^2.+(h.time*(history.y-lunit*s[0].y))^2.+(h.time*(history.z-lunit*s[0].z))^2.)^0.5 
;allstruct.radius[0] = ((h.time*(history.x-lunit*s[grpind].x))^2.+(h.time*(history.y-lunit*s[grpind].y))^2.+(h.time*(history.z-lunit*s[grpind].z))^2.)^0.5 
allstruct.radius[0] = ((history.x-lunit*s[grpind].x)^2.+(history.y-lunit*s[grpind].y)^2.+(history.z-lunit*s[grpind].z)^2.)^0.5 
allstruct.vx[0] = history.vx
allstruct.vy[0] = history.vy
allstruct.vz[0] = history.vz
allstruct.x[0] = history.x
allstruct.y[0] = history.y
allstruct.z[0] = history.z
print,'Memory I am using',memory(/high)
FOR j=1L,nsteps-1 do begin
   help,h
   history = mrdfits(list[j],1)
   rtipsy, gtp_file[j], h,g,d,s
   readcol, statfile[j], grp, mvir, rvir0, format='l,x,x,x,x,f,f', /silent
   grpind = where(grp eq halo[j])
   iord = read_lon_array(iord_files[j])
   inds = binfind(iord, alliords)
   exist = where(inds NE -1, comp=del)
   ninds = n_elements(exist)
   allstruct[exist].temp[j] = history[exist].temp
   allstruct[exist].mass[j] = history[exist].mass
   allstruct[exist].vx[j] = history[exist].vx
   allstruct[exist].vy[j] = history[exist].vy
   allstruct[exist].vz[j] = history[exist].vz
   allstruct[exist].x[j] = history[exist].x
   allstruct[exist].y[j] = history[exist].y
   allstruct[exist].z[j] = history[exist].z
   ;allstruct[exist].rho[j] = history[exist].rho*densunit/(h.time^3.)
   allstruct[exist].rho[j] = history[exist].rho*densunit
   ;allstruct[exist].entropy[j] = alog10(history[exist].temp^1.5/(history[exist].rho*densunit/h.time^3.))
   allstruct[exist].entropy[j] = alog10(history[exist].temp^1.5/(history[exist].rho*densunit))
   allstruct[exist].grp[j] = history[exist].haloid
   ;allstruct[exist].radius[j] = ((h.time*(history[exist].x-lunit*s[grpind].x))^2.+(h.time*(history[exist].y-lunit*s[grpind].y))^2.+(h.time*(history[exist].z-lunit*s[grpind].z))^2.)^0.5 
   allstruct[exist].radius[j] = ((history[exist].x-lunit*s[grpind].x)^2.+ (history[exist].y-lunit*s[grpind].y)^2.+(history[exist].z-lunit*s[grpind].z)^2.)^0.5 
   if del[0] ne -1 then begin
   allstruct[del].temp[j] = -1
   allstruct[del].mass[j] = -1
   allstruct[del].rho[j] = -1
   allstruct[del].entropy[j] = -1
   allstruct[del].grp[j] = -1
   allstruct[del].radius[j] = -1
   allstruct[del].vx[j] = -1
   allstruct[del].vy[j] = -1
   allstruct[del].vz[j] = -1
   allstruct[del].x[j] = -1
   allstruct[del].y[j] = -1
   allstruct[del].z[j] = -1
   endif
ENDFOR
stop

;Write for testing purposes.  But this file is huge.  Delete if possible.
mwrfits, allstruct, 'grp1.allgas.entropy.test.fits', /create

;Pass structure to find_smooth, to find early and smoothly accreted
;gas particles
;stop
jump1: allstruct = mrdfits('grp1.allgas.entropy.test.fits',1)
find_smooth, allstruct, halo
stop

;Now pass the structure to find_shocked, to find the shocked smooth particles 
if keyword_set(wmap1) then find_shocked, allstruct, haloidoutfile, /wmap1 else $
find_shocked, allstruct, haloidoutfile
stop

;Now sort the results to get clumpy and unshocked
all = mrdfits('grp1.allgas.iord.test.fits',0)
early = mrdfits('early.iord.test.fits',0)
test = binfind(early,all)
late = all(where(test eq -1))
sm = mrdfits('smooth.accr.iord.test.fits',0)
test = binfind(sm,late)
clumpy = late(where(test eq -1))
mwrfits, clumpy, 'clumpy.accr.iord.test.fits', /create
shock = mrdfits('shocked.iord.test.fits',0)
if total(shock) gt 0.0 then begin
   test = binfind(shock,sm)
   unshock = sm(where(test eq -1))
endif else unshock = sm
mwrfits, unshock, 'unshock.iord.test.fits', /create

;Get info on accretion time, mass at accretion 
accretion, haloidoutfile, allstruct

end

