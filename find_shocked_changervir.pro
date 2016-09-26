pro find_shocked_changervir, phase, haloidoutput, WMAP1=wmap1
;# just input your favourite galaxy: everything else shoule be aok!
;# "MW1.1024g1bwK"," ", " "
;#

readcol, haloidoutput, files, halo, format='a,l'
files = reverse(files)
halo = reverse(halo)
nsteps = n_elements(files)
gtpfile = files+'.amiga.gtp'
statfile = files+'.amiga.stat'

;Get tvir
tvir = fltarr(nsteps)
rvir = fltarr(nsteps)
mu=0.7229612 ; effective atomic weight
for j=0,nsteps-1 do begin
  readcol, statfile[j], mvir, rvir0, format='x,x,x,x,x,f,f', skipline=halo[j], numline=1, /silent
  rtipsy, gtpfile[j], h,g,d,s
  tvir[j] = 1.7306277e-07*mvir*mu/(rvir0/1000.)
  rvir[j] = rvir0
endfor

;phase = mrdfits('grp1.allgas.entropy.fits',1)
smaccr = mrdfits('smooth.accr.iord.fits',0)

;Concentrate only on smoothly accreted gas, to avoid complications with SF in 
;accreted halos
ind = findex(phase.iord, smaccr)
smphase = phase[ind]

; Get a at each timestep to calculate min SF density

test = n_elements(smphase[0].grp)
if test ne nsteps then stop

; Assume there's an entropy floor set by the UV background temp and the 
; omega_baryon.  When gas shocks, it jumps in density by at least a factor 
; of 4, and to T_vir.  
;omega_m = rho_m(t)/rho_c(t)

;IF keyword_set(wmap1) then begin
; om_m = 0.3  ;for this cosmology
; z1 = [5.89,4,3,2,1,0.5,0]
; rho_c = rho_crit(z1)
; rho_m = om_m*rho_c
; t_uv = [18781.3,24207.5,25697.6,24083.4,19720.4,16283.6,11554.3]
; smean = alog10(t_uv^1.5/rho_m)
; a = fltarr(nsteps)
;  FOR j=0,nsteps-1 do begin
;    rtipsy, gtpfile[j], h,g,d,s
;    a[j] = h.time
;  ENDFOR
; z = (1./a)-1.
; linterp, z1, t_uv, z, tuv
; linterp,z1,smean,z,ent  ;This interpolates mean entropy
; sshock = alog10((3.*tvir/8.)^1.5/(4.*om_m*rho_crit(z)))
; sshock2 = alog10(tvir^1.5/(4.*om_m*rho_crit(z)))
;; Note that I used om_m here rather than om_b, but om_b factors into both sshock and 
; ent in the same way, so that deltas comes out the same.
; deltas = sshock-ent
; deltas2 = sshock2-ent
;ENDIF ELSE BEGIN
 om_m = 0.24  ;for this cosmology
 ;z1 = [5.89,4,3,2,1,0.5,0]
 z1 = [20,15,10,9,8,7,6,5,4,3,2,1.5,1,0.75,0.5,0.25,0.125,0.1,0.05,0.]
 rho_c = wmap3_rho_crit(z1)
 rho_m = om_m*rho_c
 ;t_uv = [18781.3,24207.5,25697.6,24083.4,19720.4,16283.6,11554.3]  from original uniform box in grad school -- think it had a bug
 ;t_uv = [19687.7,23646.5,30831.8,33438.5,31486.3,24760.4,19516.2,21699.4,24856.0,25089.7,23668.0,21868.2,19205.8,17608.6,15748.6,13548.1,12297.7,12031.6,11507.2,10983.2] ;uniform box from 5 years ago, now with more z.  Still think it has a bug.
 ;t_uv = [12261.2,13805.8,15046.8,15173.2,14982.1,14307.4,13047.0,12830.3,12558.5,11425.5,9598.2,8501.2,7220.7,6513.1,5720.4,4837.1,4353.8,4251.3,4052.3,3843.9]  ;uniform cosmo25 box, primordial cooling + .UV file.  Note much cooler than the box I built ages ago.  Whoops.
 ;t_uv = [89.3,49.7,24.2,16917.3,11959.4,10710.4,11093.1,12313.0,12259.1,11040.0,9298.0,8245.4,6990.2,6305.8,5541.4,4692.3,4229.0,4130.8,3940.1,3739.9];uniform cosmo25 box using OLD cooltable_xdr (12Mb version).  Very different at high z than using the .UV file (see email from Sijing).
 t_uv = [89.3,6964.4,8819.9,8695.6,8511.5,8298.4,7910.5,7859.6,12068.7,13292.9,10419.8,8967.6,7402.0,6589.2,5701.5,4743.1,4233.4,4126.8,3921.6,3709.0];uniform cosmo25 box using NEW cooltable_xdr (25Mb version).  Very different at high z than using the .UV file (see email from Sijing).
 smean = alog10(t_uv^1.5/rho_m)
 a = fltarr(nsteps)
 FOR j=0,nsteps-1 do begin
   rtipsy, gtpfile[j], h,g,d,s
   a[j] = h.time
 ENDFOR
 z = (1./a)-1.
 linterp,z1,smean,z,ent  ;This interpolates mean entropy
 sshock = alog10((3.*tvir/8.)^1.5/(4.*om_m*wmap3_rho_crit(z)))
 deltas = sshock-ent
 sshock2 = alog10(tvir^1.5/(4.*om_m*wmap3_rho_crit(z)))
 deltas2 = sshock2-ent
;ENDELSE


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Have decided to follow particles until they reach 30 kpc 
; from the galaxy center, or rvir (whichever is smaller).
; First, find the step where the particles reach 30 kpc. If rvir 
; at that step is smaller than 30 kpc, then for those particles 
; find where particle enters rvir instead.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ngas = n_elements(smphase.iord)
accr=lonarr(ngas)+999
lowr=fltarr(ngas)
FOR j=0L,ngas-1 do lowr[j] = min(smphase[j].radius(where(smphase[j].radius gt 0.0)))
ind = where(lowr le max(0.2*rvir), nind, comp=ind2, ncomp=nind2)
for k=0L,nind-1 do accr[ind(k)] = min(where(smphase[ind(k)].radius lt max(0.2*rvir) and smphase[ind(k)].radius ge 0.))
FOR j=0L,nind2-1 do accr[ind2(j)] = min(where(smphase[ind2(j)].radius eq lowr[ind2(j)]))  
change = where(rvir lt max(0.2*rvir))
first = max(change)+1
ind = where(accr le first,nind)
for k=0L,nind-1 do accr[ind(k)] = min(where(smphase[ind(k)].grp eq halo))

mwrfits, accr, 'tracedtostep_change.fits', /create


testmin = min(where(deltas gt 0))
testmin2 = min(where(deltas2 gt 0))
if testmin lt 0 then begin
  print, 'Shocked entropy is always less than the ambient entropy.  Stopping'
  iords=0.
  shinfo=0.
  goto, noshock
  ;stop 
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
shocked=-1
timestep=-1
radius=-1
allshock=-1
allradii=-1
alltime=-1


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find shocks 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
late = where(accr gt testmin, nlate)
IF nlate gt 0 then begin
FOR j=0L,n_elements(late)-1 do begin
  ent_shock = alog10((3.*tvir[testmin+1:accr(late[j])]/8.)^1.5/(4.*smphase[late(j)].rho[testmin:accr(late[j])-1]))
  delta_ent = smphase[late(j)].entropy[testmin+1:accr(late[j])]
  ;What if require here that temp must be 3/8 tvir as well?
  test = where(delta_ent gt ent_shock and smphase[late(j)].temp[testmin+1:accr(late[j])] gt 3.*tvir[testmin+1:accr(late[j])]/8.,ntest)
  if ntest gt 0 then begin
    prior = test-1
    nprior = lonarr(ntest)
    for k=0,ntest-1 do nprior[k] = where(test[0:ntest-1] eq prior[k])
    allinds = where(nprior eq -1,nall)
    if nall ne 0 then begin
     ind0 = where(test[allinds] ne 0,n0)
     if n0 ne 0 then begin
     allinds = allinds(ind0)
     allshock = [allshock, late(j)]
     allradii = [allradii,smphase[late(j)].radius[testmin+1+test(allinds)]]
     alltime = [alltime,test(allinds)+testmin+1]
     endif
    endif
    rtest = smphase[late(j)].radius[testmin+1:accr(late[j])]
    in = where(rtest-rvir[testmin+1:accr(late[j])] lt 0, nin)
    if nin gt 0 then begin
      keep = intersect(in,test)
      prior = min(keep)-1
      sudden = where(test eq prior)
      if sudden eq -1 and prior ne -1 then begin
        shocked = [shocked, late(j)]
        timestep = [timestep, min(keep)+testmin+1]
        radius = [radius,smphase[late(j)].radius[min(keep)+testmin+1]]
      endif
    endif
  endif
ENDFOR
ENDIF
    
;Eliminate any double entries that may occur in the high z and low z searches
shocked = shocked[1:n_elements(shocked)-1]
timestep = timestep[1:n_elements(timestep)-1]
radius = radius[1:n_elements(radius)-1]
allshock = allshock[1:n_elements(allshock)-1]
allradii = allradii[1:n_elements(allradii)-1]
alltime = alltime[1:n_elements(alltime)-1]
accrvir = accr[shocked]
accrall = accr[allshock]
print, 'Number of gas particles traced:  ', ngas
print, 'Number of particles shocked within Rvir:  ', n_elements(shocked)
print, 'Number of particles shocked anywhere:  ', n_elements(allradii)

;Sort them (just in case - they need to be sorted for other programs)
sorted = sort(shocked)
siord = shocked[sorted]
saccr = accrvir[sorted]
sstep = timestep[sorted]
srad = radius[sorted]
iords = smphase[siord].iord
shinfo = fltarr(4,n_elements(shocked))
shinfo[0,*] = sstep
shinfo[1,*] = z[sstep]
shinfo[2,*] = srad
shinfo[3,*] = saccr
noshock:
mwrfits, iords, 'shockedchange.iord.fits', /create
mwrfits, shinfo, 'shockedchange.info.fits', /create

;sorted = sort(allshock)
;siord = allshock[sorted]
;saccr = accrall[sorted]
;sstep = alltime[sorted]
;srad = allradii[sorted]
;iords = smphase[siord].iord
;mwrfits, iords, 'smooth.shocked.iord.all.fits', /create
;mwrfits, z[sstep], 'smooth.shocked.zshock.all.fits', /create
;mwrfits, saccr, 'smooth.shocked.tracestep.all.fits', /create
;mwrfits, srad, 'smooth.shocked.radius.all.fits', /create

end

