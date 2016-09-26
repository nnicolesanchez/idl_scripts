pro find_smooth_changa,phase,haloid
; phase is a structure containing all relevant info, passed to find_smooth by accrmode
; haloid starts with high z and ends with z=0 
; 
; Finds the gas that ends up in the galaxy at z=0 that was smoothly
; accreted - i.e., never in another halo


; ****************************************************
; define some stuff
;*******************************************************

ind = -1
test0 = where(haloid ne 1,ntest,comp=grp1,ncomp=ngrp1)
nsteps = n_elements(haloid)
print,nsteps
;The next two lines eliminate particles that already belong to 
;the galaxy at the first timestep (so we don't know if they were 
;smoothly accreted or clump), and has the added bonus of removing
;some clumpy particles from the rest of the search

early = where(phase.grp[0] eq haloid[1], comp=free) ;Originally at haloid[0]; editted for h277 (which begins at timestep 120)
print,'number of early particles',n_elements(early)

if early[0] ne -1 then mwrfits, phase[early].iord, 'early.iord.fits', /create
phase = phase[free]
print,'number of gas particles',n_elements(phase)
ngal = n_elements(phase.iord)
print,'number of galaxies',ngal

;FOR j=0L,ngal-1 do begin
;  h1 = histogram(phase[j].grp, min=0, max=1)
;  if total(h1) eq nsteps then begin
;    ind = [ind,j] 
;    continue
;  endif 
;  if total(h1) ne ngrp1 then continue else begin
;    test2 = where(phase[j].grp[test0] eq haloid[test0],ntest2)
;    if ntest2 ne ntest then continue else ind = [ind,j]
;  endelse
;ENDFOR
; start at 0L
FOR j=ngal-100,ngal-1 do begin
  test1 = where(phase[j].grp eq haloid, n1)
  print,'length of test1:',n_elements(test1)
  print,'n1:',n1
  if n1 eq 0 then continue
; if n1 eq 1 and test1[0] eq nsteps-1 then begin
;    ind = [ind,j]
;    continue
;  endif
;  if n1 eq 1 then begin
;    prior = indgen(test1)
;    test2 = where(phase[j].grp[prior] ne 0, n2)
;    if n2 gt 0 then continue else ind = [ind,j]
;    continue
;  endif
; Check it's been there two consecutive steps
;  ingal = lonarr(n1)
;  for k=0,n1-2 do ingal[k] = test1[k+1]-test1[k]
;  first = min(where(ingal eq 1))+test1[0]
;
; One timestep only is ok
;  first = min(test1)
  
  print,'I think this needs to not be zero:',n_elements(test1)
  prior = indgen(min(test1))
  test2 = where(phase[j].grp[prior] ne 0, n2)
  if n2 gt 0 then continue else ind = [ind,j]
ENDFOR


ind = ind[1:n_elements(ind)-1]
iord_smooth = phase[ind].iord
mwrfits, iord_smooth, 'smooth.accr.iord.fits', /create

end

