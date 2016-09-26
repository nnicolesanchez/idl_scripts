pro uniquegasiords

;If gasiords.all.nonunique concatinate all gasiords.#### into 1 file
readcol,'gasiords.all.nonunique',gasiord,meaten,format='l,d',/silent

uniquegas = gasiord[uniq(gasiord,sort(gasiord))]
nunique = n_elements(uniquegas)
gasiords = lonarr(nunique)
meatens = dblarr(nunique)

for i=0,nunique-1 do begin
    howmany = where(gasiord eq uniquegas[i],n)
    gasiords[i] = uniquegas[i]
    if n gt 1 then meatens[i] = total(meaten[howmany]) else meatens[i] = meaten[i]
        
endfor

openw,lun,'gasiords.all',/get_lun
for j=0,nunique-1 do printf,lun,gasiords[j],meatens[j]
close,lun
free_lun,lun






end
