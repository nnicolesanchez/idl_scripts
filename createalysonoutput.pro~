pro createalysonoutput

; creates a file for alyson to read in to halo tracing code.
; format is "filename", "BH halo id" for each BH


readcol,'files.list',file,format='a',/silent
nfiles = n_elements(file)
restore,'bhhalo.sav'
nbh = n_elements(bhhalo)

for i=0,nbh-1 do begin
	notzero = where(bhhalo[i].bhiord ne 0,n0)
	if n0 eq 0 then continue
	bhiord = bhhalo[i].bhiord[notzero[0]]
	filename = 'haloid.'+trim(bhiord)+'.dat'
	haloid = bhhalo[i].haloid
	openw,1,filename,width=500
        ;stop
	for j=0,nfiles-1 do printf,1,file[j],haloid[j]
	close,1
endfor


end
