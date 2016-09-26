pro zto2_halo1_gasaccr1

; trying to copy Alyson's script
; halo1_gasaccr1.sh
; in IDL
; this creates *history.fits files for each output

readcol,'../files.list',files,format='a',/silent
nfiles = n_elements(files)
grpfiles = files+'.amiga.grp'
iords = files+'.iord'
igasords = files+'.igasorder'
print,nfiles

for i=0,nfiles-1 do begin
    rtipsy,'../'+files[i],h,g,d,s
    sind = where(s.tform gt 0)
    grp = read_lon_array(grpfiles[i])
    grp = grp[(h.ngas+h.ndark):(h.n-1)]
    grp = grp[sind]
    dwfstar = where(grp EQ 1)
    dwf = sind[dwfstar]+(h.ngas+h.ndark)
    iord = read_lon_array(iords[i])
    igasord = read_lon_array(igasords[i])
    orig_iord = iord(dwf)
    orig_igasord = igasord(dwf)
    ;print,i
    ;stop
    find_history, files[i], 'star', 'tipsy.units.idl', IND_STARS=dwfstar, IGASORD_STARS=orig_igasord, ORIG_IORD=orig_iord, minpart=64
    print,i
endfor

end
