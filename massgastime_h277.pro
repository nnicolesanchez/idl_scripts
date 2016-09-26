pro massgastime_h277,plotter=plotter

restore,'../bhorbits.sav'
tunit=1.223128407d18/3.1556926d7/1d9 ; Gyr
munit=1.84793d16
restore,'allgas.sav'
readcol,'out.distance',w,bhiords,w,eattime,w,gasiords,dsquare,w,w,smooth,w,eatenmass,format='a,l,a,f,a,l,f,a,a,f,a,f',/silent
;readcol,'out.mdot',w,bhiords,w,time,w,w,mdot,w,w,mdotedd,w,a, $
;  format='a,l,a,f,a,a,f,a,a,f,a,f',/silent
;readcol,'gasiords.17399719',gasid,meaten,format='l,f',/silent
readcol,'out.dm',w,bhiords,w,delta,w,dmtime,w,dm,format='a,l,a,f,a,f,a,f',/silent
loadct,39 ;color table
;the rest of the columns do not interest me.

; praise the gods, dmtime and eattime are equal arrays.
; this means that dm and gasiords correspond one to one
; YAY
; match dm with gasiords and I'll have clumpy and smooth.

; need to sync up times and redshifts
;match,time,eattime,t1,t2
;time=time[t1]
;a = a[t1]
;z = 1./a-1. ; redshift.
;eattime=eattime[t2]
;gasiords=gasiords[t2]
; this will miss a few but I will let that go for now.

; times and redshift, for the z axis on the plots
readcol,'../times.list',gyr,z,/silent
zplotnames=['6','4','2','1','0.5','0']
zplotnum=[6,4,2,1,0.5,0]
tickvalues=interpol(gyr,z,zplotnum) 
tickvalues[n_elements(tickvalues)-1]=gyr[n_elements(gyr)-1]

; redshift that gas was accreted by galaxy
zs = allgas.zaccr
ngas=n_elements(allgas.iords)

!p.multi=[0,1,1]
!p.charsize=1.5
bin=0.3


;  ### divide into clumpy and smooth for BH accretion ###
; clumpy
match_multi,allgas.clumpy,gasiords,c11
clumpybhtimes = eattime[c11]
clumpybhmass = eatenmass[c11]
; shocked
match_multi,allgas.shocked,gasiords,sh11
shockedbhtimes = eattime[sh11]
shockedbhmass = eatenmass[sh11]
if total(sh11) le -0.5 then numshocked = 0 else numshocked=n_elements(sh11)
; cold
match_multi,allgas.cold,gasiords,cold11
coldbhtimes = eattime[cold11]
coldbhmass = eatenmass[cold11]
; early
match_multi,allgas.early,gasiords,early1
earlybhtimes = eattime[early1]
earlybhmass = eatenmass[early1]
stop

; #### sort so they are in chronological order ###
; clumpy
sortclumpy=sort(clumpybhtimes)
clumpybhtimes=clumpybhtimes[sortclumpy]
clumpybhmass=clumpybhmass[sortclumpy]
; cold
sortcold=sort(coldbhtimes)
coldbhtimes=coldbhtimes[sortcold]
coldbhmass=coldbhmass[sortcold]
; shocked
sortshocked=sort(shockedbhtimes)
shockedbhtimes=shockedbhtimes[sortshocked]
shockedbhmass=shockedbhmass[sortshocked]
; early
sortearly=sort(earlybhtimes)
earlybhtimes=earlybhtimes[sortearly]
earlybhmass=earlybhmass[sortearly]

; ### then add for cumulative totals ###
coldtotal=total(coldbhmass,/cumulative)*munit
shockedtotal=total(shockedbhmass,/cumulative)*munit
clumpytotal=total(clumpybhmass,/cumulative)*munit
alltotal = total(eatenmass,/cumulative)*munit
earlytotal= total(earlybhmass,/cumulative)*munit

; ### for mass vs time plots,  need to be bin data. ###
nbins = 100
maxtime = max(eattime)
mintime = min(eattime)
bindivision = (maxtime-mintime)/nbins
tinds = indgen(nbins)*bindivision+mintime  ; indices for time array.
bindmtime = tinds
; now bin masses.  take the sum of the mass in each bin set.
bindm = fltarr(nbins)
bincoldbhmass =  fltarr(nbins)
binclumpybhmass =  fltarr(nbins)
binshockedbhmass =  fltarr(nbins)

for i=0,nbins-2 do begin
    inds = where(eattime ge (i*bindivision+mintime) AND eattime lt ((i+1)*bindivision)+mintime,ninds)
    if ninds ne 0 then bindm[i] = total(eatenmass[inds])
    coldinds = where(coldbhtimes ge (i*bindivision+mintime) AND coldbhtimes lt ((i+1)*bindivision)+mintime,ncoldinds)
    if ncoldinds ne 0 then bincoldbhmass[i] = total(coldbhmass[coldinds])
    clumpyinds = where(clumpybhtimes ge (i*bindivision+mintime) AND clumpybhtimes lt ((i+1)*bindivision)+mintime,nclumpyinds)
    if nclumpyinds ne 0 then binclumpybhmass[i] = total(clumpybhmass[clumpyinds])
    shockedinds = where(shockedbhtimes ge (i*bindivision+mintime) AND shockedbhtimes lt ((i+1)*bindivision)+mintime,nshockedinds)
    if nshockedinds ne 0 then binshockedbhmass[i] = total(shockedbhmass[shockedinds])
endfor

;### mass vs time plot ###
if keyword_set(plotter) then paperplot,filename='massvstime.ps'

plot,bindmtime*tunit,bindm*munit,xtit='time (Gyr)',ytit='mass accreted by central BH',/ylog,yra=[1,1e6]
   oplot,bindmtime*tunit,binclumpybhmass*munit,color=140;,psym=3
   oplot,bindmtime*tunit,bincoldbhmass*munit,color=90;,psym=3
   oplot,bindmtime*tunit,binshockedbhmass*munit,color=240;,psym=3

legend,['all','cold','shocked','clumpy'],lines=[0,0,0,0],color=[0,90,240,140],charsize=1,/right,/bottom
if keyword_set(plotter) then paperplot,/close

;stop


if keyword_set(plotter) then paperplot,filename='h277_allmassgas.ps'

;### cumulative mass plot ###

; use orbit file to get max BH mass
;bhorbitmass = fltarr(n_elements(bhorbits[1].x))
;for i=0,n_elements(bhorbits[35].x)-1 do bhorbitmass[i] =
;max(bhorbits.mass[i])*munit
;stop
bhorbitmass = bhorbits[223].mass*munit
bhorbittime = bhorbits[223].time*tunit
;bhorbittime=bhorbittime[0:36]
;bhorbitmass=bhorbitmass[0:36]
;print,bhorbittime
;print,bhorbitmass

bhorbittime=bhorbittime[where( bhorbittime ne 0)]
bhorbitmass=bhorbitmass[where( bhorbitmass ne 0)]

;print,bhorbittime

plot,dmtime*tunit,alltotal,xtit='Time (Gyr)',ytit='Cumulative BH Mass',/ylog,yra=[100,1e7],xstyle=8
   oplot,clumpybhtimes*tunit,clumpytotal,color=95;,psym=3
   oplot,bhorbittime,bhorbitmass,lines=1,color=20
   oplot,coldbhtimes*tunit,coldtotal,color=70;,psym=3
    if numshocked gt 0 then  oplot,shockedbhtimes*tunit,shockedtotal,color=210;,psym=3
;   oplot,earlybhtimes*tunit,earlytotal,color=50;,psym=3

help,bhorbittime
help,bhorbitmass 

legend,['total BH mass','total gas','clumpy','cold','shocked'],lines=[2,0,0,0,0],color=[20,0,95,70,210],charsize=1,/right,/bottom

axis,xaxis=1,xtit='Redshift',xtickv=tickvalues,xticks=n_elements(zplotnames)-1,xtickname=zplotnames

xyouts,0.1,4e6,' h277',charsize=1.5,/data

if keyword_set(plotter) then paperplot,/close

; save this info for future use in other plot!
bhaccretemass={dmtime:dmtime,alltotal:alltotal,clumpybhtimes:clumpybhtimes,coldbhtimes:coldbhtimes,shockedbhtimes:shockedbhtimes,clumpytotal:clumpytotal,coldtotal:coldtotal,shockedtotal:shockedtotal,earlybhtimes:earlybhtimes,earlytotal:earlytotal}
save,bhaccretemass,filename='bhaccretemass.sav'

stop

if keyword_set(plotter) then paperplot,filename='allmassgas_NSF.ps'

;### cumulative mass plot for NSF proposal ###

plot,dmtime*tunit,alltotal,xtit='Time (Gyr)',ytit='Cumulative BH Mass (M'+sunsymbol()+')',/ylog,yra=[1,1e7];,title='h277';,psym=3
   oplot,clumpybhtimes*tunit,clumpytotal,color=140,linestyle=1
   oplot,coldbhtimes*tunit,coldtotal,color=90,linestyle=2

legend,['Total Gas','Cold','Clumpy'],lines=[0,2,1],color=[0,90,140],charsize=1.3

if keyword_set(plotter) then paperplot,/close


;stop

; ### test that the parts add up to the whole. ###

nclumpy = n_elements(clumpytotal)
ncold = n_elements(coldtotal)
nshocked = n_elements(shockedtotal)
ndm = n_elements(alltotal)
nearly = n_elements(early1)

print,'does',alltotal[ndm-1],' equal',clumpytotal[nclumpy-1]+coldtotal[ncold-1]+shockedtotal[nshocked-1]+earlytotal[nearly-1]


print,'gas fractions '
print,'cold/total = ',total(coldbhmass)/total(dm)
print,'clumpy/total = ',total(clumpybhmass)/total(dm)
print,'shocked/total = ',total(shockedbhmass)/total(dm)
print,'early/total = ',total(earlybhmass)/total(dm)



;stop


end
