pro adot,doprint=doprint,plotter=plotter

n = 9000

a = (findgen(n)+1.)/n
z = 1./a - 1.
time = wmap3_lookback(z)
time = reverse(time)/1d9

hubble = wmap3_hubble(z)

adot = hubble*a

!p.multi=[0,1,3]
!p.charsize=2 
if keyword_set(plotter) then paperplot3,filename='adot.ps'
plot,time,a,xtit='time',ytit='a'

plot,time,adot,xtit='time',ytit='da/dz',yra=[0,1000];/ylog


deriva = deriv(time/38.759428,a)
  plot,time,deriva,xtit='time',ytit='da/dt'
if keyword_set(plotter) then paperplot3,/close


if keyword_set(doprint) then begin
    openw,lun,'adot.dat',/get_lun
    printf,lun,'  a        dadz       dadt '
    for i=0,n-1 do printf,lun,a[i],adot[i],deriva[i]
    close,lun
    free_lun,lun
endif

end
