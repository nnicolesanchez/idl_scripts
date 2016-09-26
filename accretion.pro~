pro accretion, haloidoutfile, phase

readcol, haloidoutfile, files, halo, format='a,l'
files = reverse(files)
halo = reverse(halo)
nsteps = n_elements(files)
gtp_file = files+'.amiga.gtp'

;phase has been shortened to exclude particles already in the galaxy at the first step (early)
ngas = n_elements(phase.iord)

accr = lonarr(ngas)
FOR j=0L,ngas-1 do begin
   ingal = where(phase[j].grp eq halo)
   ;The following ensures that the particle has been in the halo for 
   ;two consecutive steps to count as accreted
   ;test = lonarr(nsteps-1)
   ;for i=0,n_elements(ingal)-2 do test[i] = ingal[i+1]-ingal[i]
   ;ind = min(where(test eq 1))
   ;if ind ne -1 then accr[j] = ingal[ind]
   ;if min(ingal) eq nsteps-1 then accr[j] = nsteps-1  ;for those particles that enter at the final step 
   ;I should test how this result differs from just 1 step: min(ingal) = accr
   accr[j] = min(ingal)
ENDFOR
z = fltarr(nsteps)
for j=0L,nsteps-1 do begin
  rheader, gtp_file[j], h
  z[j] = (1./h.time)-1.
endfor
accrtime = z[accr]
mwrfits, accrtime, 'grp1.accrz.fits', /create

massataccr = fltarr(n_elements(accr))
vxataccr = fltarr(n_elements(accr))
vyataccr = fltarr(n_elements(accr))
vzataccr = fltarr(n_elements(accr))
xataccr = fltarr(n_elements(accr))
yataccr = fltarr(n_elements(accr))
zataccr = fltarr(n_elements(accr))
FOR j=1,nsteps-1 do begin
  ind1 = where(accr eq j, nind)
  if nind ne 0 then begin
    massataccr[ind1] = phase[ind1].mass[j]
    vxataccr[ind1] = phase[ind1].vx[j]
    vyataccr[ind1] = phase[ind1].vy[j]
    vzataccr[ind1] = phase[ind1].vz[j]
    xataccr[ind1] = phase[ind1].x[j]
    yataccr[ind1] = phase[ind1].y[j]
    zataccr[ind1] = phase[ind1].z[j]
   endif
ENDFOR

mwrfits, massataccr, 'grp1.mass_at_accr.fits', /create
mwrfits, vxataccr, 'grp1.vx_at_accr.fits', /create
mwrfits, vyataccr, 'grp1.vy_at_accr.fits', /create
mwrfits, vzataccr, 'grp1.vz_at_accr.fits', /create
mwrfits, xataccr, 'grp1.x_at_accr.fits', /create
mwrfits, yataccr, 'grp1.y_at_accr.fits', /create
mwrfits, zataccr, 'grp1.z_at_accr.fits', /create

end
