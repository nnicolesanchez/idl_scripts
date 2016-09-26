pro dttest

; examine dt problem/bug from June 2013 code

readcol,'out.mdot',w,w,w,time,w,w,mdot,w,w,mdotedd,w,a,/silent,format='a,a,a,f,a,a,f,a,a,f,a,f'


readcol,'out.dm',w,bhiord3,w,delta,w,time1,w,dm,w,dE,format='a,l,a,f,a,f,a,f,a,f',/silent
stop
!p.multi=[0,1,2]

dmdt = dm/delta

plot,time,mdot,/ylog,psym=3
  oplot,time1,dmdt,color=240,psym=3



; let's do cumulative

dt= time-shift(time,1)
n=n_elements(mdot)
dt=dt[1:n-1]
mdot=mdot[1:n-1]
dm=dm[1:n-1]
time=time[1:n-1]
time1=time1[1:n-1]
 
totalmdot= total(mdot*dt,/cum)
totaldm = total(dm,/cum)

plot,time,totalmdot,/ylog
  oplot,time1,totaldm,color=240


moo=mdot/dmdt[1:n-1]

plot,time,moo

stop
end
