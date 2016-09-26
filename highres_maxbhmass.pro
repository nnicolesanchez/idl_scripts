pro highres_maxbhmass

munit  = 1.84793d16 ; msun

;readcol,'../bhfitsfiles.dat',files,format='a'
;nfiles = n_elements(files)

file = mrdfits('BHAccLog.fits',1)
uniq_afac_ind = uniq(file.afac)
uniq_afac = file.afac[uniq_afac_ind]
uniq_z = 1. /(uniq_afac - 1.)
;print,'Unique redshifts :',uniq_z[0:100]

for i=0,10 do begin
   afac = where(file.afac eq uniq_afac[i])
   bhs_at_afac = file.afac[afac]
   print,'Redshift :',(1. /(afac - 1.))
   print,'Numnber of BHs at z :',n_elements(bhs_at_afac)
   bhs_iord_at_afac = file.bhiord[afac]
   print,'BH IDs at z :',bhs_iord_at_afac


;	index         = where(max(bhfile.mass))
;	maxbhmass_id  = bhfile.bhiord[index]
;	maxbhmass     = bhfile.mass[index]
;
;	print,'Print BH Index',index
;	print,'Print BH Max Iord',maxbhmass_id
;	print,'Print BH Max Mass',maxbhmass

endfor

end
