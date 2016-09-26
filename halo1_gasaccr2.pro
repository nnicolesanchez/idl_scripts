pro halo1_gasaccr2

orig_iord = mrdfits('grp1.allgas.iord.fits',0)
readcol,'files.list',files,format='a',/silent
nfiles = n_elements(files)

for i=0,nfiles-1 do begin 
   find_history, files[i], 'allgas', 'tipsy.units.idl', ORIG_IORD=orig_iord, minpart=64, /del
   print,i
endfor


end
