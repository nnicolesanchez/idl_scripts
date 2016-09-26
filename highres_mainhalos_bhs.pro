pro highres_mainhalos_bhs

restore,'bhhalo.sav'
munit=1.84793d16
readcol,'timesteps.list',timestep,format='l'
ntimesteps = n_elements(timestep)

print,'Stuff in bhhalo.sav'
help,bhhalo

print,'Number of BHs'
help,bhhalo.haloid[0]
print,'Number of timesteps'
help,bhhalo[0].haloid

;bhhalo[# of BHs].haloid[timesteps]
;Note for haloids: 0 means no halo and 1 should be most masive, then
;2, then 3

max_bh_in_halo1 = strarr(ntimesteps)
max_bh_in_halo2 = strarr(ntimesteps)
max_bh_in_halo3 = strarr(ntimesteps)
max_bh_in_halo4 = strarr(ntimesteps)
max_bh_in_halo5 = strarr(ntimesteps)

max_bhmass_halo1 = strarr(ntimesteps)
max_bhmass_halo2 = strarr(ntimesteps)
max_bhmass_halo3 = strarr(ntimesteps)
max_bhmass_halo4 = strarr(ntimesteps)
max_bhmass_halo5 = strarr(ntimesteps)

for i = 0,ntimesteps-1 do begin      
   print,'For Timestep',timestep[i]
   haloid_timestep    = bhhalo.haloid[i]
   bh_in_timestep     = bhhalo.bhiord[i]
   bhmass_in_timestep = bhhalo.mass[i] * munit

   print,'Most massive BH in Halo 1'
   BH_in_halo1    = where(haloid_timestep eq 1)
   haloid1 = haloid_timestep[BH_in_halo1]
   bhiord1 = bh_in_timestep[BH_in_halo1]
   bhmass1 = bhmass_in_timestep[BH_in_halo1]
   max_bhmass1 = max(bhmass1)
   if max_bhmass1 eq 0 then print, 'No Black Hole'
   max_bh_ind1 = where(bhmass1 eq max_bhmass1)
   max_bhiord1 = bhiord1[max_bh_ind1]
   print,max_bhiord1
   print,max_bhmass1

   print,'Most massive BH in Halo 2'
   BH_in_halo2    = where(haloid_timestep eq 2)
   haloid2 = haloid_timestep[BH_in_halo2]
   bhiord2 = bh_in_timestep[BH_in_halo2]
   bhmass2 = bhmass_in_timestep[BH_in_halo2]
   max_bhmass2 = max(bhmass2)
   if max_bhmass2 eq 0 then print, 'No Black Hole'
   max_bh_ind2 = where(bhmass2 eq max_bhmass2)
   max_bhiord2 = bhiord2[max_bh_ind2]
   print,max_bhiord2
   print,max_bhmass2

   print,'Most massive BH in Halo 3'
   BH_in_halo3    = where(haloid_timestep eq 3)
   haloid3 = haloid_timestep[BH_in_halo3]
   bhiord3 = bh_in_timestep[BH_in_halo3]
   bhmass3 = bhmass_in_timestep[BH_in_halo3]
   max_bhmass3 = max(bhmass3)
   if max_bhmass3 eq 0 then print, 'No Black Hole'
   max_bh_ind3 = where(bhmass3 eq max_bhmass3)
   max_bhiord3 = bhiord3[max_bh_ind3]
   print,max_bhiord3
   print,max_bhmass3

   print,'Most massive BH in Halo 4'
   BH_in_halo4    = where(haloid_timestep eq 4)
   haloid4 = haloid_timestep[BH_in_halo4]
   bhiord4 = bh_in_timestep[BH_in_halo4]
   bhmass4 = bhmass_in_timestep[BH_in_halo4]
   max_bhmass4 = max(bhmass4)
   if max_bhmass4 eq 0 then print, 'No Black Hole'
   max_bh_ind4 = where(bhmass4 eq max_bhmass4)
   max_bhiord4 = bhiord4[max_bh_ind4]
   print,max_bhiord4
   print,max_bhmass4

   print,'Most massive BH in Halo 5'
   BH_in_halo5    = where(haloid_timestep eq 5)
   haloid5 = haloid_timestep[BH_in_halo5]
   bhiord5 = bh_in_timestep[BH_in_halo5]
   bhmass5 = bhmass_in_timestep[BH_in_halo5]
   max_bhmass5 = max(bhmass5)
   if max_bhmass5 eq 0 then print, 'No Black Hole'
   max_bh_ind5 = where(bhmass5 eq max_bhmass5)
   max_bhiord5 = bhiord5[max_bh_ind5]
   print,max_bhiord5   
   print,max_bhmass5

   max_bh_in_halo1[i]  = max_bhiord1
   max_bh_in_halo2[i]  = max_bhiord2
   max_bh_in_halo3[i]  = max_bhiord3
   max_bh_in_halo4[i]  = max_bhiord4
   max_bh_in_halo5[i]  = max_bhiord5

   max_bhmass_halo1[i] = max_bhmass1
   max_bhmass_halo2[i] = max_bhmass2
   max_bhmass_halo3[i] = max_bhmass3
   max_bhmass_halo4[i] = max_bhmass4
   max_bhmass_halo5[i] = max_bhmass5

endfor

openw,lun,'maxbh_in_mainhalos.dat',/get_lun,width=500
printf,lun,'    Timestep  ','MaxBH in H1    ','MaxBHmass H1 ','MaxBH in H2    ','MaxBHmass H2 ','MaxBH in H3    ','MaxBHmass H3 ','MaxBH in H4    ','MaxBHmass H4 ','MaxBH in H5    ','MaxBHmass H5'
for i=0,ntimesteps-1 do printf,lun,timestep[i],max_bh_in_halo1[i],max_bhmass_halo1[i],max_bh_in_halo2[i],max_bhmass_halo2[i],max_bh_in_halo3[i],max_bhmass_halo3[i],max_bh_in_halo4[i],max_bhmass_halo4[i],max_bh_in_halo5[i],max_bhmass_halo5[i]
close,lun
free_lun,lun

min_bh_in_halo1 = strarr(ntimesteps)
min_bh_in_halo2 = strarr(ntimesteps)
min_bh_in_halo3 = strarr(ntimesteps)
min_bh_in_halo4 = strarr(ntimesteps)
min_bh_in_halo5 = strarr(ntimesteps)

min_bhdist_halo1 = strarr(ntimesteps)
min_bhdist_halo2 = strarr(ntimesteps)
min_bhdist_halo3 = strarr(ntimesteps)
min_bhdist_halo4 = strarr(ntimesteps)
min_bhdist_halo5 = strarr(ntimesteps)

for i = 0,ntimesteps-1 do begin
   print,'For Timestep',timestep[i]
   haloid_timestep    = bhhalo.haloid[i]
   bh_in_timestep     = bhhalo.bhiord[i]
   bhdist_in_timestep = bhhalo.halodist[i]

   print,'Most central BH in Halo 1'
   BH_in_halo1    = where(haloid_timestep eq 1)
   haloid1 = haloid_timestep[BH_in_halo1]
   bhiord1 = bh_in_timestep[BH_in_halo1]
   bhdist1 = bhdist_in_timestep[BH_in_halo1]
   min_bhdist1 = min(bhdist1)
   if min_bhdist1 eq 0 then print, 'No Black Hole'
   min_bh_ind1 = where(bhdist1 eq min_bhdist1)
   min_bhiord1 = bhiord1[min_bh_ind1]
   print,min_bhiord1
   print,min_bhdist1

   print,'Most central BH in Halo 2'
   BH_in_halo2    = where(haloid_timestep eq 2)
   haloid2 = haloid_timestep[BH_in_halo2]
   bhiord2 = bh_in_timestep[BH_in_halo2]
   bhdist2 = bhdist_in_timestep[BH_in_halo2]
   min_bhdist2 = min(bhdist2)
   if min_bhdist2 eq 0 then print, 'No Black Hole'
   min_bh_ind2 = where(bhdist2 eq min_bhdist2)
   min_bhiord2 = bhiord2[min_bh_ind2]
   print,min_bhiord2
   print,min_bhdist2

   print,'Most central BH in Halo 3'
   BH_in_halo3    = where(haloid_timestep eq 3)
   haloid3 = haloid_timestep[BH_in_halo3]
   bhiord3 = bh_in_timestep[BH_in_halo3]
   bhdist3 = bhdist_in_timestep[BH_in_halo3]
   min_bhdist3 = min(bhdist3)
   if min_bhdist3 eq 0 then print, 'No Black Hole'
   min_bh_ind3 = where(bhdist3 eq min_bhdist3)
   min_bhiord3 = bhiord3[min_bh_ind3]
   print,min_bhiord3
   print,min_bhdist3

   print,'Most central BH in Halo 4'
   BH_in_halo4    = where(haloid_timestep eq 4)
   haloid4 = haloid_timestep[BH_in_halo4]
   bhiord4 = bh_in_timestep[BH_in_halo4]
   bhdist4 = bhdist_in_timestep[BH_in_halo4]
   min_bhdist4 = min(bhdist4)
   if min_bhdist4 eq 0 then print, 'No Black Hole'
   min_bh_ind4 = where(bhdist4 eq min_bhdist4)
   min_bhiord4 = bhiord4[min_bh_ind4]
   print,min_bhiord4
   print,min_bhdist4

   print,'Most central BH in Halo 5'
   BH_in_halo5    = where(haloid_timestep eq 5)
   haloid5 = haloid_timestep[BH_in_halo5]
   bhiord5 = bh_in_timestep[BH_in_halo5]
   bhdist5 = bhdist_in_timestep[BH_in_halo5]
   min_bhdist5 = min(bhdist5)
   if min_bhdist5 eq 0 then print, 'No Black Hole'
   min_bh_ind5 = where(bhdist5 eq min_bhdist5)
   min_bhiord5 = bhiord5[min_bh_ind5]
   print,min_bhiord5
   print,min_bhdist5

   min_bh_in_halo1[i]  = min_bhiord1
   min_bh_in_halo2[i]  = min_bhiord2
   min_bh_in_halo3[i]  = min_bhiord3
   min_bh_in_halo4[i]  = min_bhiord4
   min_bh_in_halo5[i]  = min_bhiord5

   min_bhdist_halo1[i] = min_bhdist1
   min_bhdist_halo2[i] = min_bhdist2
   min_bhdist_halo3[i] = min_bhdist3
   min_bhdist_halo4[i] = min_bhdist4
   min_bhdist_halo5[i] = min_bhdist5

endfor

openw,lun,'minbhdist_mainhalos.dat',/get_lun,width=500
printf,lun,'    Timestep   ','CBH in H1  ','MinBHdist H1  ','CBH in H2  ','MinBHdist H2  ','CBH in H3  ','MinBHdist H3  ','CBH in H4  ','MinBHdist H4  ','CBH in H5  ','MinBHdist H5'
for i=0,ntimesteps-1 do printf,lun,timestep[i],min_bh_in_halo1[i],min_bhdist_halo1[i],min_bh_in_halo2[i],min_bhdist_halo2[i],min_bh_in_halo3[i],min_bhdist_halo3[i],min_bh_in_halo4[i],min_bhdist_halo4[i],min_bh_in_halo5[i],min_bhdist_halo5[i]
close,lun
free_lun,lun

end


