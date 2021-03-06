;# This program determines which bhs merged into our main halo prior to
;# timestep 336. We assume we already know which black hole is in our main halo
;# after timestep 336 (360 and above). For the case of CHANGA h258,
;# this bhiord is  BHIORD 43553282 (B82) so we search for which 
;# BHs merged in early z (before timestep 336) to form B82.

;# The steps we take:
;# Read in fits file with information about galaxy sim
;#     Currently examining high rez, Changa simulation of h258
;# Parse out array with B82s info
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

;# Vanderbilt Univ.  -- VPAC40: /home/sanchenn/IDL/highres_bhtrace_latez.pro
;# N. Nicole Sanchez -- Last edit: January 23, 2016
pro highres_bhtrace_latez,bhofinterest,bimodalcutoff

;# Input Parameters
if keyword_set(bhofinterest) then begin
   bhofinterest=bhofinterest  
endif else begin
   print,'What BH is most central/massive?' 
   stop
endelse
if keyword_set(bimodalcutoff) then begin
   bimodalcutoff = bimodalcutoff 
endif else begin
   bimodalcutoff = 1.0d4        ;For h258, & should work usually
endelse
munit      = 1.84793d16         ; M_sun
loadct,4

;# Read in fits files
file       = mrdfits('../BHAccLog.fits',1)
;help,file

;# Calculate z from a factor to sort by
z          = (1. /file.afac) - 1.
print,'Starting redshift',z[0]

;# Set Main BH (for high res h258: B82)
MBH_ind    = where(file.bhiord eq bhofinterest) 
MBH_bhiord = file.bhiord[MBH_ind]
MBH_mbhi   = file.mbhi[MBH_ind] * munit   ; M_sun
MBH_mbhf   = file.mbhf[MBH_ind] * munit   ; M_sun
MBH_dmacc  = file.dmacc[MBH_ind] * munit  ; M_sun
MBH_mgi    = file.mgi[MBH_ind] * munit    ; M_sun
MBH_z      = z[MBH_ind]
;# NOTE: dmacc is change in gas mass from accr, does not include 
;# gas accr from BH; will not use in rest of this script
;# (Included above only to remind you, Nicole)

;# THERE IS WEIRDNESS: Sometimes mgi = 0, these are bad lines
;# Remove these then redefine all lines
goodlines  = where(MBH_mgi ne 0)
MBH_bhiord = MBH_bhiord[goodlines]
MBH_mbhf   = MBH_mbhf[goodlines]
MBH_mbhi   = MBH_mbhi[goodlines]
MBH_z      = MBH_z[goodlines]

;# Sort in time
z_ind     = sort(MBH_z)
MBH_mbhf  = MBH_mbhf[z_ind]
MBH_mbhi  = MBH_mbhi[z_ind]
MBH_z     = MBH_z[z_ind]
nmbhf = n_elements(MBH_mbhf)
print,'Number of BH Masses',nmbhf

;# Calculate deltambhf and find bimodal cutoff
MBH_deltambhf = fltarr(n_elements(MBH_mbhf))
for i=1,n_elements(MBH_mbhf)-1 do begin
   MBH_deltambhf[i] = MBH_mbhf[i-1] - MBH_mbhf[i]
endfor

;# Plot distribution of ∆mbhf to determine where your mass cutoff is
plothist,alog10(MBH_deltambhf[where(MBH_deltambhf gt 0.0)]),/ylog,xtitle='Log of Changes in Final BH Masses in M_sun',yrange=[0.1,100000]
stop

;# Plot mbhf across time to see all mergers (big steps)
plot,MBH_z,MBH_mbhf,xtitle='Redshift',ytitle='Log of Central BH Mass in Msun',/ylog,xrange=[25,0]
stop

;# Create array with deltambhf values greater than cutoff value
;# Cut off value for h258 = 10^4
MBH_merger_masses = [0]
MBH_merger_times  = [0]
MBH_mass_at_merger = [0]
;print,1.0d4
for i=0,n_elements(MBH_deltambhf)-1 do begin
   if (MBH_deltambhf[i] gt bimodalcutoff) then begin
      new_dmacc = MBH_dmacc[i]
      new_deltambhf = MBH_deltambhf[i]
      new_z = MBH_z[i]
      new_mbhf = MBH_mbhf[i]
      MBH_merger_masses = [MBH_merger_masses,new_deltambhf] 
      MBH_merger_times   = [MBH_merger_times,new_z]
      MBH_mass_at_merger = [MBH_mass_at_merger,new_mbhf]
      ;print,new_deltambhf
      ;print,new_z
      ;print,new_mbhf
   endif
endfor

;print,MBH_mergers[0]
;print,n_elements(MBH_mergers)
;# NOTE: I'm doing a minus 2 at the end of these array redefinitions
;# because there is some sort of bug in the highest redshift (large
;# z) values.
MBH_merger_masses  = MBH_merger_masses[1:n_elements(B52_merger_masses)-1]
MBH_merger_times   = MBH_merger_times[1:n_elements(B52_merger_times)-1]
MBH_mass_at_merger = MBH_mass_at_merger[1:n_elements(B52_mass_at_merger)-1]
print,MBH_merger_masses
print,'Number of main mergers',n_elements(MBH_merger_masses)
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
   bh_mbhf = file.mbhf[bh_ind] * munit   ; M_sun
   bh_z    = z[bh_ind]
   z_ind   = sort(bh_z)
   bh_mbhf = bh_mbhf[z_ind]
   bh_z   = bh_z[z_ind]
   nbh     = n_elements(bh_ind)
   ;print,'Number of times BH appears in',nbh

   if (bh_z[0] ne MBH_z[0]) then begin
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
;#print,bh_maybe_id
;#stop


for i=0,n_elements(bh_maybe_mass)-1 do begin
   for j=0,n_elements(MBH_merger_masses)-1 do begin
;      if (bh_maybe_mass[i] ge B52_merger_masses[j]-3) and (bh_maybe_mass[i] le B52_merger_masses[j]+3) then begin
;         print,bh_maybe_id[i]
;         print,bh_maybe_mass[i]
;         print,bh_maybe_time[i]
;      endif
      if (bh_maybe_time[i] ge MBH_merger_times[j]-0.02) and (bh_maybe_time[i] le MBH_merger_times[j]+0.02) then begin
         print,'Disappearing BH id',bh_maybe_id[i]
         print,'Disappearing BH mass',bh_maybe_mass[i]
         print,'Disappearing BH redshift',bh_maybe_time[i]

         print,'Main BH mass at this time',MBH_mass_at_merger[j]

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
