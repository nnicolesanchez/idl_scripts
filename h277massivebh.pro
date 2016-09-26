;Finds the BH closest to the center of halo 1 for /astro1/nicole/otherh277/bhhalo.sav

pro h277massivebh

restore,'bhhalo.sav'
halo1 = where(bhhalo.haloid[133] eq 1, nin1) 

;Figure out how many BH's are in halo1 (second # in array;
;first # is # of timesteps)
;help,bhhalo[halo1].bhiord

printbhmassive=fltarr(134)
printbhmassiord=fltarr(134)

for i=0,133 do begin

	bhiords   	   = where(bhhalo[halo1].bhiord[i] ne 0)
	bhmassive 	   = max(bhhalo[halo1[bhiords]].halomass[i]) 
	bhdist		   = bhhalo[halo1[bhiords]].bhiord[i]
	bhmassiord 	   = where(bhhalo[halo1].halomass[i] eq bhmassive)
	bhdistiord         = bhhalo[halo1[bhmassiord]].bhiord[i]
	printbhmassive[i]  = bhmassive
	printbhmassiord[i] = bhmassiord
endfor
;print,printbhdistiord
;help,bhdist
;help,bhdistiord

write_csv,'massbhiords.sav',printbhmassiord,printbhmassive

end
