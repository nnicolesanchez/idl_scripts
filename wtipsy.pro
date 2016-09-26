PRO wtipsy,outfile,header,catg,catd,cats,STANDARD=standard
;This program takes structures like those read in via rtipsy, and
;writes them to native binary output or standard binary output if
;/STANDARD is set.  This assumes you're working on a linux box

if (N_PARAMS() eq 0) then begin
  print, "wtipsy.pro -- Writes tipsy files with structures as read in"
  print, "using rtipsy.pro."
  print, ""
  print, "Usage: "
  print, "        wtipsy, outfilename ,header ,g ,d ,s, [/STANDARD]"
  print, ""
  print, "Input parameters: "
  print, "  outfilename  string containing name of output file"
  print, "  g,d,s     gas, dark and star structures"  
  print, "            if no gas or stars just include a dummy variable"    
  print, "Please read rtipsy.pro for the structure definitions"
  print, "  /STANDARD will write the output in standard format if you're using a i386 box"
  return
endif

OPENW,1,outfile


IF (keyword_set(standard) EQ 0) THEN BEGIN
;NATIVE
    WRITEU,1,header[0]

    IF (header.ngas) GT 0 THEN BEGIN
        FOR i=0l,header.ngas-1 DO BEGIN
            WRITEU,1,catg[i]
        ENDFOR
    ENDIF
    
    
    IF (header.ndark) GT 0 THEN BEGIN
        FOR i=0l,header.ndark-1 DO BEGIN
            WRITEU,1,catd[i]
        ENDFOR    
    ENDIF
    
    IF (header.nstar) GT 0 THEN BEGIN
        FOR i=0l,header.nstar-1 DO BEGIN
            WRITEU,1,cats[i]
        ENDFOR
    ENDIF    

ENDIF ELSE BEGIN

    WRITEU,1,swap_endian(header[0])

    dummy=1L
    WRITEU,1,swap_endian(dummy)

    IF (header.ngas) GT 0 THEN BEGIN
        FOR i=0l,header.ngas-1 DO BEGIN
            WRITEU,1,swap_endian(catg[i])
        ENDFOR
    ENDIF
        
    IF (header.ndark) GT 0 THEN BEGIN
        FOR i=0l,header.ndark-1 DO BEGIN
            WRITEU,1,swap_endian(catd[i])
        ENDFOR    
    ENDIF
    
    IF (header.nstar) GT 0 THEN BEGIN
        FOR i=0l,header.nstar-1 DO BEGIN
            WRITEU,1,swap_endian(cats[i])
        ENDFOR
    ENDIF
ENDELSE

CLOSE,1
FREE_LUN, 1

END
