;******************************************************************************
;This function analyzes the history of gas and star particles by
;taking a list of input particles and outputing a .fits file containing
;their properties at a given timestep.  
;
;REQUIREMENTS:
;
;     This procedure requires a minimum of one array if only gas is 
;     being traced (orig_iord), but three arrays if stars or gas+stars
;     are being traced (orig_iord, ind_stars, igasord_stars)  
;
;     orig_iord: the iord values of the objects you'd like to track, as 
;     given in the *iord files.  If you are tracing both gas and stars, 
;     this array will include the iord values for all particles 
;
;     igasord_stars: if you are tracing star particles, this is the 
;     igasorder value of those stars.  If you are tracing both gas and 
;     stars simultaneously, this array is ONLY for the stars.  If you 
;     are tracing only gas, do not include this input
;
;     ind_stars: If you are tracing stars (or gas+stars), this array 
;     is the index (the order in the list) of the STARS ONLY.  That is,
;     this array has a minimum value of 0 and maximum value of h.nstar.
;     Ignore h.ngas+h.ndark.  This is not the index to the entire 
;     (gas+dm+star) list, but ONLY the star list.  This is true even 
;     if you are tracing gas and stars simultaneously.  If you are 
;     tracing gas only, do not include this input
;
;     The program also deals with AMIGA .grp outputs (.stat could be 
;     added).  The output structure will contain the ID of the halo  
;     (history.haloid) that a star/gas particle belonged to at that 
;     timestep.  (Thus, AMIGA should be run on the timestep of choice 
;     before running this script.)
;
;     The program writes a fits file for each timestep using the 
;     name given in the simoutputlist file plus a unique string 
;     'tonamefiles' to identify a specific trace from others (e.g., gas 
;     vs stars of a MW type galaxy, traced separately).
;
; OPTIONAL KEYWORDS
;     The above assumes that you are tracing stars or gas particles that
;     still exist at the final timestep.  In some cases, you may wish to 
;     follow the history of gas particles that get deleted by the final 
;     step.  In that case, use keyword DEL.  This will write the properties 
;     of particles as in the standard trace, but when a particle no longer
;     exists it is left as zeros.
;
;     A number of our AMIGA outputs have found halos down to a minimum of 
;     16 particles.  However, our mass function may only be resolved down 
;     to 64 particles.  If you'd like to set a minimum number of particles 
;     for AMIGA to use, set MINPART to the minimum number of particles. 
;     If keyword TYPE isn't also set, then MINPART will require that there be 
;     a minimum number of MINPART dm particles.  Other options are 'gas', 'star', 
;     'baryon', and 'tot.'
;
;******************************************************************************
pro find_history, simoutputlist, tonamefiles, tipsyunitsfile, ORIG_IORD=orig_iord, IND_STARS=ind_stars, IGASORD_STARS=igasord_stars, MINPART=minpart, DEL=del, TYPE=type
; tonamefiles is a string to append to the output files, so multiple searches can 
; be done on multiple criteria and kept separate.
;
; tipsyunitsfile is tipsy.units.idl: this is the same as tipsy.units, but with the 
; ()'s removed so that idl can read it 
; 
; ind_stars should be an array at starts at h.ngas+h.ndark = 0 (that is, the maximum 
; length of ind_stars is the length of h.nstar for the timesteps they are pulled
; from). 
;
; orig_iord includes the iord values of both gas and star particles, if both are being 
; traced
;
; igasord_stars is the igasorder values for stars being traced 
;
; The number of elements in igasord_star and ind_stars should be equal, but 
; orig_iord will be longer if you are tracing both gas and stars  


readcol, tipsyunitsfile, lengthunit, massunit, velunit, format='d,d,f',/silent
root = simoutputlist 

;get number of timesteps and desired particles
ntimesteps=N_ELEMENTS(root)
norig=N_ELEMENTS(orig_iord)

;create structure to hold information of particles over time:
;history structure is subscripted by number of desired particles, many
;elements within it are subscripted by number of timesteps
history=REPLICATE({iord:0l, igasorder:0l, mark:0l, mass:DBLARR(ntimesteps), x:DBLARR(ntimesteps), y:DBLARR(ntimesteps), z:DBLARR(ntimesteps), vx:DBLARR(ntimesteps), vy:DBLARR(ntimesteps), vz:DBLARR(ntimesteps), rho:DBLARR(ntimesteps), temp:DBLARR(ntimesteps), metallicity:DBLARR(ntimesteps), haloid:LONARR(ntimesteps)}, norig)

history.iord=orig_iord
IF n_elements(igasord_stars) ne 0 THEN history.igasorder=igasord_stars 
;For each output,
FOR i = 0l, ntimesteps - 1 DO BEGIN
  ;Read in the .iord file.
  iord=read_lon_array(root[i]+'.iord')
  ;Read in valuable information from the tipsy file
  rtipsy, root[i], h,g,d,s
  s = s(where(s.tform gt 0))
  haloind=read_lon_array(root[i]+'.amiga.grp')

IF keyword_set(MINPART) then begin
  readcol, root[i]+'.amiga.stat', grp, ntot, ngas, nstar, ndark, format='l,l,l,l,l', /silent
  IF keyword_set(TYPE) then begin
    if type eq 'tot' then mingrp = grp(where(ntot lt minpart, ngrp)) 
    if type eq 'gas' then mingrp = grp(where(ngas lt minpart, ngrp)) 
    if type eq 'star' then mingrp = grp(where(ngas lt minpart, ngrp)) 
    if type eq 'baryon' then mingrp = grp(where(ngas+nstar lt minpart, ngrp)) 
    if type eq 'dm' then mingrp = grp(where(ndark lt minpart, ngrp)) 
    if ngrp ne 0 then begin
      FOR j=0,n_elements(mingrp)-1 do begin
      ind = where(haloind eq mingrp[j])
      haloind[ind] = 0
      ENDFOR
    endif
  ENDIF ELSE begin
    test = where(ndark lt minpart, ngrp)
    if ngrp ne 0 then begin
      mingrp = grp(where(ndark lt minpart, ngrp)) 
      FOR j=0,n_elements(mingrp)-1 do begin
      ind = where(haloind eq mingrp[j])
      haloind[ind] = 0
      ENDFOR
    endif
  ENDELSE
ENDIF
  
  nstars = n_elements(ind_stars)
  ngas = norig-nstars
  IF ngas NE 0 then begin
    orig_iord_gas = orig_iord[0:ngas-1]
    gas = indgen(ngas, /long)
  ENDIF
  IF nstars NE 0 then begin
    orig_iord_star = orig_iord[ngas:norig-1]
    stars = indgen(nstars, /long)+ngas
  ENDIF
  ;For stars,
  IF nstars GT 0 then begin
    ; Find which are still stars and which are gas at this step
    gasprog = where(orig_iord_star GT max(iord), nprog, comp=oldstars, ncomplement=nold)
    ind2stars = ind_stars[oldstars]	;Found the stars, now find gas indices
    if nprog ne 0 then begin
    ind2gas = ind_stars[gasprog]		;For indexing into the history array
    igasorder = igasord_stars[gasprog]
    iord=iord[0:h.ngas]	;shorten for faster reading
    gasind = FINDEX(iord, igasorder) 
    endif 
    ;Now fill up the history structure
    ;first for stars
    history[stars(0:nold-1)].mark[i]=long(h.ngas+h.ndark+ind2stars)+1
    FOR j=0L,nold-1 do begin
        history[stars(j)].mass[i]=s[ind2stars(j)].mass*massunit 
        history[stars(j)].x[i]=s[ind2stars(j)].x*lengthunit
        history[stars(j)].y[i]=s[ind2stars(j)].y*lengthunit
        history[stars(j)].z[i]=s[ind2stars(j)].z*lengthunit
        history[stars(j)].vx[i]=s[ind2stars(j)].vx*velunit
        history[stars(j)].vy[i]=s[ind2stars(j)].vy*velunit
        history[stars(j)].vz[i]=s[ind2stars(j)].vz*velunit
        history[stars(j)].rho[i]=0.
        history[stars(j)].temp[i]=0.
        history[stars(j)].metallicity[i]=s[ind2stars(j)].metals
        history[stars(j)].haloid[i]=haloind[ind2stars(j)+h.ngas+h.ndark]
    ENDFOR
    ;now for stars that are still gas at this step
    if nprog ne 0 then begin
    history[stars(nold:nstars-1)].mark[i]=long(gasind)+1
    FOR j=0l,nprog-1 do begin
    history[stars(j+nold)].mass[i]=g[gasind(j)].mass*massunit
    history[stars(j+nold)].x[i]=g[gasind(j)].x*lengthunit 
    history[stars(j+nold)].y[i]=g[gasind(j)].y*lengthunit
    history[stars(j+nold)].z[i]=g[gasind(j)].z*lengthunit
    history[stars(j+nold)].vx[i]=g[gasind(j)].vx*velunit
    history[stars(j+nold)].vy[i]=g[gasind(j)].vy*velunit
    history[stars(j+nold)].vz[i]=g[gasind(j)].vz*velunit
    history[stars(j+nold)].rho[i]=g[gasind(j)].dens;*densityunit
    history[stars(j+nold)].temp[i]=g[gasind(j)].tempg
    history[stars(j+nold)].metallicity[i]=g[gasind(j)].zmetal
    history[stars(j+nold)].haloid[i]=haloind[gasind(j)] 
    ENDFOR     
    endif
  ENDIF 
  IF ngas GT 0 then begin
    iord=iord[0:h.ngas]	;shorten for faster reading
    gasind = binfind(iord, orig_iord_gas)
    del = where(gasind eq -1)
    exist = where(gasind ne -1)
    gasind = gasind(where(gasind ne -1))
    history[exist].mark[i]=long(gasind)+1
    IF keyword_set(del) then begin
    FOR j=0L,n_elements(gasind)-1 do begin
    history[exist(j)].mass[i]=g[gasind(j)].mass*massunit
    history[exist(j)].x[i]=g[gasind(j)].x*lengthunit
    history[exist(j)].y[i]=g[gasind(j)].y*lengthunit
    history[exist(j)].z[i]=g[gasind(j)].z*lengthunit
    history[exist(j)].vx[i]=g[gasind(j)].vx*velunit
    history[exist(j)].vy[i]=g[gasind(j)].vy*velunit
    history[exist(j)].vz[i]=g[gasind(j)].vz*velunit
    history[exist(j)].rho[i]=g[gasind(j)].dens;*densityunit
    history[exist(j)].temp[i]=g[gasind(j)].tempg
    history[exist(j)].metallicity[i]=g[gasind(j)].zmetal
    history[exist(j)].haloid[i]=haloind[gasind(j)]
    ENDFOR
    ENDIF ELSE BEGIN
    FOR j=0L,ngas-1 do begin
    history[gas(j)].mass[i]=g[gasind(j)].mass*massunit
    history[gas(j)].x[i]=g[gasind(j)].x*lengthunit 
    history[gas(j)].y[i]=g[gasind(j)].y*lengthunit
    history[gas(j)].z[i]=g[gasind(j)].z*lengthunit
    history[gas(j)].vx[i]=g[gasind(j)].vx*velunit
    history[gas(j)].vy[i]=g[gasind(j)].vy*velunit
    history[gas(j)].vz[i]=g[gasind(j)].vz*velunit
    history[gas(j)].rho[i]=g[gasind(j)].dens;*densityunit
    history[gas(j)].temp[i]=g[gasind(j)].tempg
    history[gas(j)].metallicity[i]=g[gasind(j)].zmetal
    history[gas(j)].haloid[i]=haloind[gasind(j)]       
    ENDFOR
    ENDELSE
  ENDIF
  outfile = root[i]+'.'+tonamefiles+'.history.fits'
  mwrfits, history, outfile, /create

ENDFOR
;stop
END

