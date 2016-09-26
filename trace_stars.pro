; This code uses the .amiga.grp files to find the progenitor halo that 
; contained most of the z=0 dm of a given halo

pro trace_stars, simname, grp, stepsize=stepsize

;files = file_search(simname+'.0*/*.amiga.stat')
readcol,'files.list',files,format='a',/silent
;command = "ls "+simname+".0*/*.amiga.stat | grep amiga.stat | sed 's/.amiga.stat//g'"
;spawn, command, filebase
filebase = files
ind2 = fltarr(2.d7,n_elements(grp))
majarr = fltarr(n_elements(files),n_elements(grp))
grpfiles = filebase+'.amiga.grp'

openw, lun, simname+'.trace_stars.dat', /get_lun
for i=n_elements(files)-1,0,-stepsize do begin
  test = read_lon_array(grpfiles[i])
  rheader, filebase[i], h
  test = test[h.ngas+h.ndark:h.n-1]
  FOR j=0,n_elements(grp)-1 do begin
    if i eq n_elements(files)-1 then begin
     ind1 = where(test eq grp[j])
     res = histogram(test[ind1], locations=x)
    endif else begin
     indmax = where(ind2[*,j] ge h.nstar, nindmax)
     if nindmax ne 0 then ind2[indmax,j] = 0
     ind1 = where(test[ind2[*,j]] ne 0 and ind2[*,j] ne 0, nind1)
     if nind1 ne 0 then res = histogram(test[ind2[ind1,j]], locations=x) else continue
    endelse
    maj = x(where(res eq max(res)))
    majarr[i,j] = maj[0]
    ;if maj[0] eq 1 then continue else begin
    nmaj = n_elements(where(test eq maj[0]))
    ind2[0:nmaj-1,j] = where(test eq maj[0])
    ind2[nmaj:2.d7-1,j] = 0
    ;endelse
  ENDFOR
;print, maj
endfor
for j=0,n_elements(grp)-1 do begin
 for i=n_elements(files)-1,0,-stepsize do printf, lun, format='(A105, 2x, I4, 2x, I4, 2x, I6)', filebase[i], grp[j], majarr[i,j]
endfor
close, lun
free_lun, lun


end

