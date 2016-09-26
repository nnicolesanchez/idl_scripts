;Finds the BH closest to the center of halo 1 for /astro1/nicole/otherh277/bhhalo.sav

pro h277centralbh

restore,'bhhalo.sav'
halo1 = where(bhhalo.haloid[133] eq 1, nin1) 

;Figure out how many BH's are in halo1 (second # in array;
;first # is # of timesteps)
;help,bhhalo[halo1].bhiord

printbhclosest=fltarr(134)
printbhdistiord=fltarr(134)

for i=0,133 do begin

	bhiords   = where(bhhalo[halo1].bhiord[i] ne 0)
	bhclosest = min(bhhalo[halo1[bhiords]].halodist[i]) 
	bhdist=bhhalo[halo1[bhiords]].bhiord[i]
	bhclosestiord = where(bhhalo[halo1].halodist[i] eq bhclosest)
	bhdistiord=bhhalo[halo1[bhclosestiord]].bhiord[i]
	printbhclosest[i] = bhclosest
	printbhdistiord[i] = bhdistiord	
endfor
;print,printbhdistiord
;help,bhdist
;help,bhdistiord

write_csv,'closebhiords.sav',printbhdistiord,printbhclosest

end
