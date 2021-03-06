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
pro highres_printredshifts

readcol,'timestep_directories.list',timesteps,format='a'
readcol,'timesteps.list',timestep_numbers,format='l'
print,timestep_numbers
print,timesteps[0]

for i=151,n_elements(timesteps)-1 do begin
   rtipsy,timesteps[i],h,g,d,s
   a = h.time
   z = (1./a) - 1
;   print,timestep_numbers[i],'is at redshift',z
   print,z
endfor

end
