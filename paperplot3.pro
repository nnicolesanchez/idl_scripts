pro paperplot3,filename=filename,close=close
;plots postscripts without having to set all of that crap.

;this one is specifically for 3x3 plots.

;if opening, type psplot3,filename='filename.ps'
;if closing, type psplot3,/close

!P.THICK=4
!P.CHARTHICK=4
!X.THICK=4
!Y.THICK=4
!X.MARGIN=[4,4]
!p.charsize=2
!p.font=0

if keyword_set(close) then begin
  device,/close
  set_plot,'x'
  !P.THICK=1
  !P.CHARTHICK=1
  !X.THICK=1
  !Y.THICK=1
  !p.charsize=1
  !p.font=-1

  return
endif else begin

  set_plot,'ps'
  device,/color,bits_per_pixel=8,filename=filename,yoff=1,/inch,ysize=9,xsize=7.5,xoff=0.5,/times
  return

endelse


end
