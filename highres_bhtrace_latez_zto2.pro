;# This program determines which bhs merged into our main halo prior to
;# timestep 144. We assume we already know which black hole is in our main halo
;# after timestep 144. For the case of h258, this bhiord is ; BHIORD 43553352
;# (B52) so we search for which BHs merged in early z (before timestep
;# 144) to form B52.

;# The steps we take:
;# Read in fits file with information about galaxy sim
;#     Currently examining high rez, Changa simulation of h258
;# Parse out array with B52s info
;#     file.bhiords    (bh id)
;#     file.afac       (a factor which calculates z)
;#     mbhf            (final black hole mass; we could also use mbhi)
;#     mbhi            (same as mbhf except shifted forward one step in array)
;#                     (ie mbhf[4] = mbhi[5])
;# Calculate difference in mass between each time step
;# Create array with these values that are greater than the bimodal break in a
;#     distribution of mbhf; Also get times
;# Determine how bhs and their info look once the bh has disappeared/merged
;# Parse out info for each of the bhs that have disappeared
;# Find the mass at which these bhs have disappeared
;# Compare to the array of times and masses
;# Determine which galaxies merged into the main halo
;# Boom.

;# Vanderbilt Univ.  -- VPAC40: /home/sanchenn/IDL/highres_bhtrace_latez_better.pro
;# N. Nicole Sanchez -- October 16, 2016
pro highres_bhtrace_latez_zto2,bhofinterest,bimodalcutoff

;# Constants
if keyword_set(bhofinterest) then bhofinterest=bhofinterest  else bhofinterest = 43553352 ;For h258
if keyword_set(bimodalcutoff) then bimodalcutoff = bimodalcutoff else bimodalcutoff = 1.0d4  ;For h258, & should work usually
munit      = 1.84793d16                        ; M_sun
loadct,4

;# Read in fits files
file       = mrdfits('BHAccLog.fits',1)
;help,file

;# Calculate z from a factor to sort by
z          = (1. /file.afac) - 1.
print,'Starting redshift',z[0]

;# Set main BH
B52_ind    = where(file.bhiord eq bhofinterest) 
B52_bhiord = file.bhiord[B52_ind]
B52_mbhi   = file.mbhi[B52_ind] * munit   ; M_sun
B52_mbhf   = file.mbhf[B52_ind] * munit   ; M_sun
B52_dmacc  = file.dmacc[B52_ind] * munit  ; M_sun
B52_mgi    = file.mgi[B52_ind] * munit    ; M_sun
B52_z      = z[B52_ind]
;# NOTE: dmacc is change in gas mass from accr, does not include 
;# gas accr from BH; will not use in rest of this script
;# (Included above only to remind you, Nicole)

;# THERE IS WEIRDNESS: Sometimes mgi = 0, these are bad lines
;# Remove these then redefine all lines
goodlines  = where(B52_mgi ne 0)
B52_bhiord = B52_bhiord[goodlines]
B52_mbhf   = B52_mbhf[goodlines]
B52_mbhi   = B52_mbhi[goodlines]
B52_z      = B52_z[goodlines]

;# Sort in time
z_ind     = sort(B52_z)
B52_mbhf  = B52_mbhf[z_ind]
B52_mbhi  = B52_mbhi[z_ind]
B52_z     = B52_z[z_ind]
nmbhf = n_elements(B52_mbhf)
print,'Number of BH Masses',nmbhf

;# Calculate deltambhf and find bimodal cutoff
B52_deltambhf = fltarr(n_elements(B52_mbhf))
for i=1,n_elements(B52_mbhf)-1 do begin
   B52_deltambhf[i] = B52_mbhf[i-1] - B52_mbhf[i]
endfor

;# Plot distribution of ∆mbhf to determine where your mass cutoff is
plothist,alog10(B52_deltambhf[where(B52_deltambhf gt 0.0)]),/ylog,xtitle='Log of Changes in Final BH Masses in M_sun',yrange=[0.1,100000]
stop

;# Plot mbhf across time to see all mergers (big steps)
plot,B52_z,B52_mbhf,xtitle='Redshift',ytitle='Log of B52 Mass in Msun',/ylog,xrange=[25,0]
stop

;# Create array with deltambhf values greater than cutoff value
;# Cut off value for h258 = 10^4
B52_merger_masses = [0]
B52_merger_times  = [0]
B52_mass_at_merger = [0]
;print,1.0d4
for i=0,n_elements(B52_deltambhf)-1 do begin
   if (B52_deltambhf[i] gt 1.0e4) then begin
      new_dmacc = B52_dmacc[i]
      new_deltambhf = B52_deltambhf[i]
      new_z = B52_z[i]
      new_mbhf = B52_mbhf[i]
      B52_merger_masses = [B52_merger_masses,new_deltambhf] 
      B52_merger_times   = [B52_merger_times,new_z]
      B52_mass_at_merger = [B52_mass_at_merger,new_mbhf]
      print,new_deltambhf
      print,new_z
      print,new_mbhf
   endif
endfor

;print,B52_mergers[0]
;print,n_elements(B52_mergers)
;# NOTE: I'm doing a minus 2 at the end of these array redefinitions
;# because there is some sort of bug in the highest redshift (large
;# z) values.
B52_merger_masses = B52_merger_masses[1:n_elements(B52_merger_masses)-1]
B52_merger_times = B52_merger_times[1:n_elements(B52_merger_times)-1]
B52_mass_at_merger = B52_mass_at_merger[1:n_elements(B52_mass_at_merger)-1]
print,B52_merger_masses
print,'Number of main mergers',n_elements(B52_merger_masses)
;stop

;# Figure out a way to determine when a black hole disappears
uniq_ind = uniq(file.bhiord)
uniq_ids = file.bhiord[UNIQ(file.bhiord, SORT(file.bhiord))]
;print,uniq_ids[0:100]
print,'Number of unique BH ids',n_elements(uniq_ids)

bh_maybe_id   = [0]
bh_maybe_mass = [0]
bh_maybe_time = [0]
print,'Number of Timesteps in Main BH',nmbhf
for i=0,n_elements(uniq_ids)-1 do begin
   bh_ind  = where(file.bhiord eq uniq_ids[i])
   bh_ids  = file.bhiord[bh_ind]
   ;print,'These should all be the same',bh_ids
   bh_mbhf = file.mbhf[bh_ind] * 1.84793d16   ; M_sun
   bh_z    = z[bh_ind]
   z_ind   = sort(bh_z)
   bh_mbhf = bh_mbhf[z_ind]
   bh_z   = bh_z[z_ind]
   nbh     = n_elements(bh_ind)
   ;print,'Number of times BH appears in',nbh

   if (bh_z[0] ne B52_z[0]) then begin
      new_id   = bh_ids[0]
      new_mass = bh_mbhf[0]
      new_time = bh_z[0]
      
      bh_maybe_id = [bh_maybe_id, new_id]
      bh_maybe_mass = [bh_maybe_mass, new_mass]
      bh_maybe_time = [bh_maybe_time, new_time]
   endif
endfor

bh_maybe_id = bh_maybe_id[1:n_elements(bh_maybe_id)-1]
bh_maybe_mass = bh_maybe_mass[1:n_elements(bh_maybe_mass)-1]
bh_maybe_time = bh_maybe_time[1:n_elements(bh_maybe_time)-1]
;print,bh_maybe_id
;stop

for i=0,n_elements(bh_maybe_mass)-1 do begin
   for j=0,n_elements(B52_merger_masses)-1 do begin
;      if (bh_maybe_mass[i] ge B52_merger_masses[j]-3) and (bh_maybe_mass[i] le B52_merger_masses[j]+3) then begin
;         print,bh_maybe_id[i]
;         print,bh_maybe_mass[i]
;         print,bh_maybe_time[i]
;      endif
      if (bh_maybe_time[i] ge B52_merger_times[j]-0.005) and (bh_maybe_time[i] le B52_merger_times[j]+0.005) then begin
         print,'Disappearing BH id',bh_maybe_id[i]
         print,'Disappearing BH mass',bh_maybe_mass[i]
         print,'Disappearing BH redshift',bh_maybe_time[i]
         print,'B52s mass at this time',B52_mass_at_merger[j]

      endif
   endfor
endfor

;####### Random notes for Nicole 
;# Remember: for h258, we know the main halo before 144 (z=10.10)
;# SO, we only have to look at mergers before this to determine
;# who the progenitors are
;IDL> rtipsy,'h258.cosmo50cmb.3072gst10bwepsK1BHC52.000144',h,g,d,s
;% Compiled module: RTIPSY.
;% Compiled module: SWAP_ENDIAN.
;IDL> help,h
;** Structure <1ff7d88>, 6 tags, length=32, data length=28, refs=1:
;   TIME            DOUBLE         0.090039161
;   N               LONG          43560070
;   NDIM            LONG                 3
;   NGAS            LONG           9429025
;   NDARK           LONG          34123819
;   NSTAR           LONG              7226
;IDL> print,h.time
;     0.090039161
;IDL> a=h.time
;IDL> z          = (1. /a) - 1.
;IDL> print,z
;       10.106278

end
