pro find_all_gas, haloidoutput
;haloidoutput includes z=0 file
;order doesn't matter 
;generally will run from high z to low z

readcol, haloidoutput, basefile, halo, format='a,l'
iordfile = basefile+'.iord'
grpfile = basefile+'.amiga.grp'
ngas = lonarr(n_elements(basefile))
for j=0L,n_elements(basefile)-1 do begin
  rheader, basefile[j], h
  ngas[j] = h.ngas
endfor

z0iord = read_lon_array(iordfile[0])
z0grp = read_lon_array(grpfile[0])
z0grp_gas = z0grp[0:ngas[0]-1]
test = where(z0grp_gas eq halo[0], ntest)
if ntest eq 0 then cum_iord = -1 else begin
  z0gas_iord = z0iord(where(z0grp_gas eq halo[0]))
  cum_iord = z0gas_iord
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; The ntest = 0 was added for tracing dSphs 06/2011
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for i=1,n_elements(basefile)-1 do begin
  grp = read_lon_array(grpfile[i])
  grpgas = grp[0:(ngas[i]-1)]
  hziord = read_lon_array(iordfile[i])
  test = where(grpgas eq halo[i], ntest)
  if ntest eq 0 then cum_iord = [cum_iord,-1] else begin
    hziord_gas = hziord(where(grpgas eq halo[i]))
    cum_iord = [cum_iord,hziord_gas]
  endelse
endfor

cum_iord = cum_iord(where(cum_iord ne -1))
sorted = cum_iord(sort(cum_iord))
unique = sorted(uniq(sorted))

outfile = 'grp1.allgas.iord.fits'
mwrfits, unique, outfile, /create


end



