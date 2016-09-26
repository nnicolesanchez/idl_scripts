pro paperplotlandscape,filename=filename,close=close
;plots postscripts without having to set all of that crap.
; does it in landscape, which  mucks up all other postscript plotting

;if opening, type psplot,filename='filename.ps'
;if closing, type psplot,/close

!P.THICK=4
!P.CHARTHICK=4
!X.THICK=4
!Y.THICK=4
!P.CHARSIZE=1.25


if keyword_set(close) then begin
  device,/close
  set_plot,'x'
!P.THICK=1
!P.CHARTHICK=1
!X.THICK=1
!Y.THICK=1


  return
endif else begin

  set_plot,'ps'
  device,/color,bits_per_pixel=8,filename=filename,/times,/landscape
  return

endelse



end
