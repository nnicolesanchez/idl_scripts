pro luminosity,plotter=plotter



readcol,'out.dm',w,bhoird,w,delta,w,time,w,dm,w,dE,format='a,l,a,f,a,f,a,f,a,f'

tunit=1.223128407d18/3.1556926d7/1d9 ; Gyr
lunit=50000.
munit=1.84793d16
Gcgs = 6.67d-8
kpccm = 3.0856776d21
msolg = 1.9891d33

ergunit = Gcgs*(msolg*munit)^2./(kpccm*lunit)
time=time*tunit
dE=dE*ergunit
delta = delta*tunit * 3.1556926d7* 1d9 ; seconds

!p.multi=[0,1,3]
!p.charsize=2

if keyword_set(plotter) then paperplot3,filename='luminosity.ps'

plot,time,dE,/ylog,xtit='Time(Gyr)',ytit='Energy (ergs)',title='h258 primary BH'

dedt = dE/delta

plot,time,dedt,xtit='Time(Gyr)',ytit='dE/dt (ergs/sec)',/ylog


; let's try to smooth and bin this plot.
nbins=300
dedt = congrid(dedt,nbins)
time=congrid(time,nbins)


plot,time,dedt,xtit='Time(Gyr)',ytit='dE/dt (ergs/sec)',/ylog,title='(binned)'


if keyword_set(plotter) then paperplot3,/close

;stop


end
