;# just input your favourite galaxy: everything else should be aok!
;# "MW1.1024g1bwK"," ", " "
;#
pro accrmode_changa, haloidoutfile, lunit, wmap1=wmap1, skip=skip, vpac38=vpac38, laststep=laststep

;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; This code creates a structure that finds the maximum temperature of a 
; gas particle in its history.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if keyword_set(wmap1) then densunit = 135.98 else densunit = 147.8344

readcol, haloidoutfile, files, halo, format='a,l'
if keyword_set(vpac38) then files='/net/vpac40/astro1/h258/3072/changa/df/'+files else print,'nope'
print,files[0]
;stop
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
print,'These are all the particles',nall
;stop
nsteps = n_elements(halo)
print,'This is nsteps',nsteps
stop
firststep=1L 

if keyword_set(laststep) then begin 
   nsteps = 170 
   firststep = 169
   print,'doing the last timestep only'
endif

print,'first step =',firststep
print,'ending step =',nsteps[n_elements(nsteps)-1]
print,'This needs to be an interger = ',nsteps
;stop

readcol, statfile[0], grp, mvir, rvir0, format='l,x,x,x,x,f,f', /silent
grpind = where(grp eq halo[0])

print,'Memory I am using',memory(/high)
;stop
allstruct = replicate({iord:0L, grp:LONARR(nsteps)}, nall);, temp:FLTARR(nsteps), rho:FLTARR(nsteps), entropy:FLTARR(nsteps), radius:FLTARR(nsteps)}, nall)mass:LONARR(nsteps), vx:FLTARR(nsteps), vy:FLTARR(nsteps), vz:FLTARR(nsteps), x:FLTARR(nsteps), y:FLTARR(nsteps), z:FLTARR(nsteps)}, nall)
stop


;;iords = replicate({iord:0L}, nall)
;;grps  = replicate({grp:LONARR(nsteps)}, nall)



print,'Memory I am using',memory(/high)
;stop
;print,allstruct[0:100].iord
allstruct.iord = alliords

;;iords = alliords
;help,allstruct.iord,alliords
;print,allstruct[0:100].iord
;print,alliords[0:100]


history = mrdfits(list[0],1)
rtipsy, gtp_file[0], h,g,d,s
;allstruct.temp[0] = history.temp
;allstruct.mass[0] = history.mass
;allstruct.rho[0] = history.rho*densunit/(h.time^3.)
;allstruct.rho[0] = history.rho*densunit
;allstruct.entropy[0] = alog10(history.temp^1.5/(history.rho*densunit/h.time^3.))
;allstruct.entropy[0] = alog10(history.temp^1.5/(history.rho*densunit))

;;grps[0] = history.haloid
allstruct.grp[0] = history.haloid
;allstruct.radius[0] = ((h.time*(history.x-lunit*s[0].x))^2.+(h.time*(history.y-lunit*s[0].y))^2.+(h.time*(history.z-lunit*s[0].z))^2.)^0.5 
;allstruct.radius[0] = ((h.time*(history.x-lunit*s[grpind].x))^2.+(h.time*(history.y-lunit*s[grpind].y))^2.+(h.time*(history.z-lunit*s[grpind].z))^2.)^0.5 
;allstruct.radius[0] = ((history.x-lunit*s[grpind].x)^2.+(history.y-lunit*s[grpind].y)^2.+(history.z-lunit*s[grpind].z)^2.)^0.5 
;allstruct.vx[0] = history.vx
;allstruct.vy[0] = history.vy
;allstruct.vz[0] = history.vz
;allstruct.x[0] = history.x
;allstruct.y[0] = history.y
;allstruct.z[0] = history.z

;stop
print,'Memory I am using',memory(/high)
;stop
;firststep = 1L
FOR j=firststep,nsteps-1 do begin
   print,j
   help,h

   history = mrdfits(list[j],1)
   rtipsy, gtp_file[j], h,g,d,s
   readcol, statfile[j], grp, mvir, rvir0, format='l,x,x,x,x,f,f', /silent
   grpind = where(grp eq halo[j])
   iord = read_lon_array(iord_files[j])
   inds = binfind(iord, alliords)
   exist = where(inds NE -1, comp=del)
   ninds = n_elements(exist)
;   allstruct[exist].temp[j] = history[exist].temp
;   allstruct[exist].mass[j] = history[exist].mass
;   allstruct[exist].vx[j] = history[exist].vx
;   allstruct[exist].vy[j] = history[exist].vy
;   allstruct[exist].vz[j] = history[exist].vz
;   allstruct[exist].x[j] = history[exist].x
;   allstruct[exist].y[j] = history[exist].y
;   allstruct[exist].z[j] = history[exist].z
   ;allstruct[exist].rho[j] = history[exist].rho*densunit/(h.time^3.)
;   allstruct[exist].rho[j] = history[exist].rho*densunit
   ;allstruct[exist].entropy[j] = alog10(history[exist].temp^1.5/(history[exist].rho*densunit/h.time^3.))
;   allstruct[exist].entropy[j] =
;   alog10(history[exist].temp^1.5/(history[exist].rho*densunit))
   
   allstruct[exist].grp[j] = history[exist].haloid
   ;allstruct[exist].radius[j] = ((h.time*(history[exist].x-lunit*s[grpind].x))^2.+(h.time*(history[exist].y-lunit*s[grpind].y))^2.+(h.time*(history[exist].z-lunit*s[grpind].z))^2.)^0.5 
;   allstruct[exist].radius[j] = ((history[exist].x-lunit*s[grpind].x)^2.+ (history[exist].y-lunit*s[grpind].y)^2.+(history[exist].z-lunit*s[grpind].z)^2.)^0.5 
   if del[0] ne -1 then begin
;   allstruct[del].temp[j] = -1
;   allstruct[del].mass[j] = -1
;   allstruct[del].rho[j] = -1
;   allstruct[del].entropy[j] = -1
   allstruct[del].grp[j] = -1
;   allstruct[del].radius[j] = -1
;   allstruct[del].vx[j] = -1
;   allstruct[del].vy[j] = -1
;   allstruct[del].vz[j] = -1
;   allstruct[del].x[j] = -1
;   allstruct[del].y[j] = -1
;   allstruct[del].z[j] = -1
   print,n_elements(allstruct)
   close,1
   endif
ENDFOR
;print,'you need to type ulimit'
;stop


;Write for testing purposes.  But this file is huge.  Delete if possible.
mwrfits, allstruct, 'grp1.allgas.entropy_iordgrps.fits';, /create
help,allstruct
print,allstruct[0:10]

;Pass structure to find_smooth, to find early and smoothly accreted
;gas particles
;stop
jump1: allstruct = mrdfits('grp1.allgas.entropy_iordgrps.fits',1)
print,memory(/high)
stop
find_smooth_changa, allstruct, halo
;stop

;Now pass the structure to find_shocked, to find the shocked smooth particles 
if keyword_set(wmap1) then find_shocked, allstruct, haloidoutfile, /wmap1 else $
find_shocked, allstruct, haloidoutfile
;stop

;Now sort the results to get clumpy and unshocked
all = mrdfits('grp1.allgas.iord.fits',0)
early = mrdfits('early.iord.fits',0)
test = binfind(early,all)
late = all(where(test eq -1))
sm = mrdfits('smooth.accr.iord.fits',0)
test = binfind(sm,late)
clumpy = late(where(test eq -1))
mwrfits, clumpy, 'clumpy.accr.iord.fits', /create
shock = mrdfits('shocked.iord.fits',0)
if total(shock) gt 0.0 then begin
   test = binfind(shock,sm)
   unshock = sm(where(test eq -1))
endif else unshock = sm
mwrfits, unshock, 'unshock.iord.fits', /create

;Get info on accretion time, mass at accretion 
; Nicole note: This uses mass,velocities, and x,y,z
accretion, haloidoutfile, allstruct

end
