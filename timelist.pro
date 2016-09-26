pro timelist

haloidfile = file_search('*haloid.dat')
readcol,haloidfile[0],file,haloid,format='a,i',/silent
file       = reverse(file)
nfiles     = n_elements(file)
redshift   = fltarr(nfiles)
numbers    = fltarr(nfiles)

for i=0,nfiles-1 do begin
   rtipsy,file[i],h,g,d,s
   print,nfiles
   a           = h.time
   z           = (1./a)-1. ; redshift.
   tunit       = 1.223128407d18/3.1556926d7/1d9 ; Gyr
   redshift[i] = z
   numbers[i]  = i
endfor
   openw,lun,'timelist.dat',/get_lun,width=500
   printf,lun,'File Name                      ','Redshift ','Numbers '
   for i=0,nfiles-1 do printf,lun,file[i],redshift[i],numbers[i]
   close,lun
   free_lun,lun
end
