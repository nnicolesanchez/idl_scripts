;# This script is a redux of Jillian's gasiords.pro (find in
;# Vanderbilt VPAC40 /home/sanchenn/IDL/gasiords.pro which
;# works for low res, Gasoline simulation outputs) 
;# This script works for high res, CHANGA files instead. 

;# The steps we take:
;#    Pull out the information from BHAccLog.fits for gasiords.pro
;#       BHIORDS          (ids for each bh)
;#       BHMBHF           (mass of black hole at end of timestep;
;#                        'final') Not used in this script. Just for out.bh
;#       GASIORD          (id of eaten mass)
;#       DMQ/EATEN MASS   (eaten mass of gas only)
;#       TIME/Z           (time/redshift)
;#    Specifically, pull out info for bhs of interest (h258: B52&B62)
;#    Put in chronological order
;#    With merger time from highres_bhtrace_latez_better.pro,
;#       Concatenate in time order
;#    Run through Jillian's gasiords.pro prescription

;# NOTE: You will concatenate the info for B62 and B52 for 
;#       these arrays of data. See highres_bhtrace_latez_better.pro 
;#       for details on how central/massive BHs were determined

;# IMPORTANT NOTE: We've determined that we will ONLY
;# concentrate on the main bh and not consider the secondary as
;# it is not comparatively interesting enough to trace

;# Vanderbilt Univ.  -- VPAC40:  /home/sanchenn/IDL/highres_gasiords.pro
;# N. Nicole Sanchez -- October 18, 2015
pro highres_gasiords,mainbh,secondarybh
  
;# Parameter and Constants 
if keyword_set(mainbh) then mainbh=mainbh  else mainbh = 43553282 ;For h258
;if keyword_set(secondarybh) then secondarybh=secondarybh else secondary=43553362
munit      = 1.84793d16                    ; M_sun
loadct,4

;# Read in fits files
file       = mrdfits('../BHAccLog.fits',1)
help,file

;# Calculate z from a factor to sort by
z          = (1. /file.afac) - 1.

;# Pull out info we want
mainbh_ind = where(file.bhiord eq mainbh)
mainbh_mgi = file.mgi[mainbh_ind]          ;Only necessary to remove badlines
mainbh_id  = file.bhiord[mainbh_ind]
mainbh_gid = file.gasiord[mainbh_ind]
mainbh_mbhf = file.mbhf[mainbh_ind] * munit ;M_sun
mainbh_dmq = file.dmacc[mainbh_ind] * munit ;M_sun
mainbh_z   = z[mainbh_ind]
;# There are some bad lines where mgi=0, remove these & redefine arrays
goodlines  = where(mainbh_mgi ne 0)
mainbh_id  = mainbh_id[goodlines]
mainbh_gid = mainbh_gid[goodlines]
mainbh_mbhf = mainbh_mbhf[goodlines]
mainbh_dmq = mainbh_dmq[goodlines]
mainbh_z   = mainbh_z[goodlines]
;# Sort in chronological order and redefine arrays
z_ind      = sort(mainbh_z)  ;returns indeces
mainbh_id  = mainbh_id[z_ind]
mainbh_gid = mainbh_gid[z_ind]
mainbh_mbhf = mainbh_mbhf[z_ind] 
mainbh_dmq = mainbh_dmq[z_ind]
mainbh_z   = mainbh_z[z_ind]
print,'Starting redshift',mainbh_z[n_elements(mainbh_z)-1]
print,'Starting mass',mainbh_mbhf[n_elements(mainbh_mbhf)-1]
print,'Ending redshift',mainbh_z[0]
print,'Ending mass',mainbh_mbhf[0]


;# If there is a secondary black hole to trace
if keyword_set(secondarybh) then begin
   secondarybh = secondarybh
   secbh_ind = where(file.bhiord eq secondarybh)
   secbh_mgi = file.mgi[secbh_ind] ;Only necessary to remove badlines
   secbh_id  = file.bhiord[secbh_ind]
   secbh_gid = file.gasiord[secbh_ind]
   secbh_mbhf = file.mbhf[secbh_ind] * munit ;M_sun
   secbh_dmq = file.dmacc[secbh_ind] * munit ;M_sun
   secbh_z   = z[secbh_ind]
   ;# There are some bad lines where mgi=0, remove these & redefine arrays
   goodlines  = where(secbh_mgi ne 0)
   secbh_id  = secbh_id[goodlines]
   secbh_gid = secbh_gid[goodlines]
   secbh_mbhf = secbh_mbhf[goodlines]
   secbh_dmq = secbh_dmq[goodlines]
   secbh_z   = secbh_z[goodlines]
   ;# Sort in chronological order and redefine arrays
   z_ind     = sort(secbh_z)
   secbh_id  = secbh_id[z_ind]
   secbh_gid = secbh_gid[z_ind]
   secbh_mbhf = secbh_mbhf[z_ind]
   secbh_dmq = secbh_dmq[z_ind]
   secbh_z   = secbh_z[z_ind]
   print,'Secondary starting redshift',secbh_z[n_elements(secbh_z)-1]
   print,'Secondary ending redshift',secbh_z[0]
   
   allbh_id  = [0]
   allbh_gid = [0]
   allbh_mbhf = [0]
   allbh_dmq = [0]
   allbh_z   = [0]
   for i=0,n_elements(mainbh_id)-1 do begin
      if (mainbh_z[i] lt secbh_z[0]) then begin
         new_id  = mainbh_id[i]
         new_gid = mainbh_gid[i]
         new_mbhf = mainbh_mbhf[i]
         new_dmq = mainbh_dmq[i]
         new_z   = mainbh_z[i]
         j=0
         ;print,new_z
      endif else begin
         new_id  = secbh_id[j]
         new_gid = secbh_gid[j]
         new_mbhf = secbh_mbhf[j]
         new_dmq = secbh_dmq[j]
         new_z   = secbh_z[j]
         ;print,new_z
         j = j+1
      endelse
      
      allbh_id  = [allbh_id,new_id]
      allbh_gid = [allbh_gid,new_gid]
      allbh_mbhf = [allbh_mbhf, new_mbhf]
      allbh_dmq = [allbh_dmq,new_dmq]
      allbh_z   = [allbh_z,new_z]
   endfor

   allbh_id  = allbh_id[1:n_elements(allbh_id)-1]
   allbh_gid = allbh_gid[1:n_elements(allbh_gid)-1]
   allbh_mbhf = allbh_mbhf[1:n_elements(allbh_mbhf)-1]
   allbh_dmq = allbh_dmq[1:n_elements(allbh_dmq)-1]
   allbh_z   = allbh_z[1:n_elements(allbh_z)-1]

endif else begin
   
   allbh_id  = mainbh_id
   allbh_gid = mainbh_gid
   allbh_mbhf = mainbh_mbhf
   allbh_dmq = mainbh_dmq
   allbh_z   = mainbh_z

endelse

;print,'Staring mass',allbh_mbhf[n_elements(allbh_mbhf)-1]
;print,'Starting redshift',mainbh_z[n_elements(mainbh_z)-1]
;print,'Ending mass',allbh_mbhf[0]
;print,'Ending Redshift',allbh_z[0]


;# The following is pulled directly from Jillian's gasiords.pro
;# for low res h258; edited for above variables 
unique_iords = allbh_id[uniq(allbh_id,sort(allbh_id))]

n=n_elements(unique_iords)
; number of bhs
print,n,' black holes'
close,1

for i=0,n-1 do begin
    thisbh = where(allbh_id eq unique_iords[i])
    eatengas = allbh_gid[thisbh]
    eatenmass = allbh_dmq[thisbh]
    ; only unique iords here.
    uniquegas=uniq(eatengas,sort(eatengas)) ; indices
    eacheatengas=eatengas[uniquegas]
    filename='gasiords.'+trim(unique_iords[i])
    openw,1,filename
    for j=0L,n_elements(eacheatengas)-1 do begin
        eacheatenmass=eatenmass[where(eatengas eq eacheatengas[j])]
        printf,1,eacheatengas[j],total(eacheatenmass)
    endfor
    close,1
endfor


spawn,'ls gasiords.* > gasiords.list'
;# End of Jillian's code


;# Another attempt to just pull out gasiords and mass;
;# Use above method, but figure out why that works better
;filename='gasiords.good.all_fromNicolesway'
;openw,1,filename
;for i=0,n_elements(allbh_gid)-1 do begin
;   printf,1,allbh_gid[i],allbh_dmq[i]
;endfor
;close,1


;# Create a pseudo out.bh with bh ids, gas ids, time, dmq
filename='out.bh'
;openw,1,filename
;for i=0,n_elements(allbh_z)-1 do begin
;   printf,1,'BHiord',allbh_id[i],'gasiord',allbh_gid[i],'eaten mass/dmq',allbh_dmq[i],'redshift',al;lbh_z[i]
;endfor
;close,1


;# Make .fits files! Too big otherwise
mwrfits, allbh_id, 'centralbh_bhiords.fits', /create
mwrfits, allbh_gid, 'centralbh_gasiords.fits', /create
mwrfits, allbh_mbhf, 'centralbh_bhmass.fits', /create
mwrfits, allbh_dmq, 'centralbh_eatenmass.fits', /create
mwrfits, allbh_z, 'centralbh_redshift.fits', /create




end
