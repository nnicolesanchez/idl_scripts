; like read_ascii_array only reads a long integer rather than float
FUNCTION read_lon_array,filename

  openr,lun,filename,/get_lun
  arrsize = 0L
  readf,lun,arrsize
;  readcol,filename[0],arrsize,format='L',numline=1
  array = lonarr(arrsize)
  readf,lun,array
  close,lun
  free_lun,lun

  return, array

END
