pro gastime_jillian,plotter=plotter

; study time gas is accreted onto halo, time BH eats it, and the
; difference.


tunit=1.223128407d18/3.1556926d7/1d9 ; Gyr
restore,'allgas.sav'
readcol,'out.distance',w,bhiords,w,eattime,w,gasiords,dsquare,w,w,smooth,w,eatenmass,format='a,l,a,f,a,l,f,a,a,f,a,f',/silent
readcol,'out.mdot',w,bhiords,w,time,w,w,mdot,w,w,mdotedd,w,a, $
  format='a,l,a,f,a,a,f,a,a,f,a,f',/silent
; need to sync up times and redshifts
; let's do an interpolation 
eata = interpol(a,time,eattime)
eatz = 1./eata-1. ; redshift.

;z = eatz
;time=eattime



; must find halo accretion history for ALL particles.

z = mrdfits('../grp1.accrz.fits',0)
tracestep = mrdfits('../tracedtostep.fits',0)
smooth=mrdfits('../smooth.accr.iord.fits',0)
all = mrdfits('../grp1.allgas.iord.fits',0)
early = mrdfits('../early.iord.fits',0)
cold = mrdfits('../unshock.iord.fits',0)
shocked = mrdfits('../shocked.iord.fits',0)
clumpy = mrdfits('../clumpy.accr.iord.fits',0)

; find early inds and remove them
match,all,early,e1,e2
notearly=all
remove,e1,notearly
; now notearly = all without early
;nnotearly = n_elements(notearly)

; redshift that gas was accreted by galaxy
zs = allgas.zaccr
notearlyallgas = where(zs lt 49.,nnotearly)
notearlyzs = zs[notearlyallgas]
ngas=n_elements(allgas.iords)
haloaccretedtime = fltarr(nnotearly)
firstaccretetime = fltarr(nnotearly)
lastaccretetime = fltarr(nnotearly)
deltat = fltarr(nnotearly)
haloaccretedtime = 13.72 - wmap3_lookback(notearlyzs)/1d9

readcol,'../times.list',redshift,time,/silent

; find the z values for smooth particles
match,notearly,smooth,sm1,sm2
smoothz = z[sm1]

match,smooth,cold,c1,c2
match,smooth,shocked,sh1,sh2

coldz = smoothz[c1]
shockedz = smoothz[sh1]

alltime = 13.72 - wmap3_lookback(z)/1d9 ; remember, z is all but early
coldhalotime = alltime[sm1[c1]]
shockedhalotime = alltime[sm1[sh1]]
match,clumpy,notearly,cl1,cl2
clumpyhalotime = alltime[cl2]

for i=0L,nnotearly-1 do begin
        ; find the redshift where particle enters halo
;	findz = findmin(zs[i],z) ; good match.
        ; convert the redshift to time
;	haloaccretedtime2[i] = eattime[findz[0]] 
	thisparticle = where(gasiords eq allgas.iords[notearlyallgas[i]],nthisparticle)
	if nthisparticle eq 0 then print,'particle info',allgas.iords[i],' lost?'
	if nthisparticle eq 0 then continue
        ;if nthisparticle gt 1 then thisparticle=thisparticle[0]
	accretetimes = eattime[thisparticle]
	firstaccretetime[i] = eattime[thisparticle[0]]
	lastaccretetime[i] = eattime[thisparticle[nthisparticle-1]]
        deltat[i]  = firstaccretetime[i]-(haloaccretedtime[i]/tunit)
        if deltat[i] lt -0.05 then stop
;        if deltat[i]*tunit lt 0 then begin
;            print,zs[i],z[findz[0]],' should be equal'
;            print,firstaccretetime[i],haloaccretedtime[i]
;            stop
;        endif
endfor

!p.multi=[0,1,2]
!p.charsize=2
bin=0.1

; divide into clumpy and smooth for galaxy entry
match,allgas.iords,allgas.clumpy,c1,c2
clumpyaccretetimes = haloaccretedtime[c1]
match,allgas.iords,allgas.shocked,sh1,sh2
shockedaccretetimes = haloaccretedtime[sh1]
match,allgas.iords,allgas.cold,cold1,cold2
coldaccretetimes = haloaccretedtime[cold1]

; divide into clumpy and smooth for BH accretion
; *** some of these times are > universe ** !!
match,gasiords,allgas.clumpy,c11,c22
clumpybhtimes = eattime[c11]
match,gasiords,allgas.shocked,sh11,sh22
shockedbhtimes = eattime[sh11]
match,gasiords,allgas.cold,cold11,cold22
coldbhtimes = eattime[cold11]

;if keyword_set(plotter) then psplot3,filename='allgastimes.ps'
;if keyword_set(plotter) then paperplot2,filename='allgastimes_halo.ps'

plothist,alltime,xall,yall,bin=bin,xtit='Time (Gyr)',title='Time all accreted gas enters halo',xra=[0.2,1.2],ytit='Number of particles',/noplot
   plothist,clumpyhalotime,xcl,ycl,/over,bin=bin,color=140,/noplot
   plothist,coldhalotime,xc,yc,/over,bin=bin,color=80,/noplot
   plothist,shockedhalotime,xsh,ysh,/over,bin=bin,color=240,/noplot

gasmass = 9.e5
plot,xall,yall*gasmass,xtit='Time (Gyr)',title='Time accreted gas enters halo',xra=[0.2,1.2],ytit='Gas Mass',/ylog,yra=[1e8,1e11]
  oplot,xcl,ycl*gasmass,color=140
  oplot,xc,yc*gasmass,color=80
  oplot,xsh,ysh*gasmass,color=240

;plothist,haloaccretedtime,bin=bin,xtit='time (Gyr)',title='Time all accreted gas enters 
;halo',xra=[0.2,1.2],ytit='Number $
;   plothist,clumpyaccretetimes,/over,bin=bin,color=140
;   plothist,coldaccretetimes,/over,bin=bin,color=80
;   plothist,shockedaccretetimes,/over,bin=bin,color=240


legend,['total','cold','shocked','clumpy'],lines=[0,0,0,0],color=[0,90,240,140],charsize=1

;if keyword_set(plotter) then paperplot,/close

;stop

bin=0.05

;if keyword_set(plotter) then paperplot,filename='allgastimes_BH.ps'

   plothist,firstaccretetime*tunit,bin=bin,xtit='time (Gyr)',title='Time gas accreted by BH',xra=[0.2,1.2],ytit='Number of accretion events'
   plothist,clumpybhtimes*tunit,/over,bin=bin,color=140
   plothist,coldbhtimes*tunit,/over,bin=bin,color=80
   plothist,shockedbhtimes*tunit,/over,bin=bin,color=240

;legend,['all gas eaten','cold ','shocked ','clumpy '],lines=[0,0,0,0],color=[0,80,240,140],charsize=1

;plothist,lastaccretetime*tunit,bin=bin,/over,color=80

;if keyword_set(plotter) then paperplot2,/close

;stop

!p.multi=[0,1,1]

;#### delta T plot ####

; what's the free fall time?
; tff = 0.25 * sqrt(3pi/2Grho)
; need to estimate uniform density rho.  at virial radius?
; of course this will change with time

; according to freefalltime.pro, the max tff = 0.68 Gyr (high z) and
; the min is 0.52 Gyr (z=0)

   tff = -20.;0.2/tunit; -20.;0.685/tunit
   ; to not cut out small/negative times, set tff = -20.

;if keyword_set(plotter) then paperplot,filename='allgastimes_deltat.ps'

bin=.03
plothist,deltat[where(deltat ge tff)]*tunit,ex,why,xtit='dt',ytit='N',bin=bin;xtit=textoidl( '\Delta')+'t',ytit='N',bin=bin
  plothist,deltat[c22[where(deltat[c22] ge tff)]]*tunit,exc,whyc,color=140,/over,bin=bin
  plothist,deltat[sh22[where(deltat[sh22] ge tff)]]*tunit,exsh,whysh,color=240,/over,bin=bin
  plothist,deltat[cold22[where(deltat[cold22] ge tff)]]*tunit,excold,whycold,color=80,/over,bin=bin

;legend,['total','cold','shocked','clumpy'],lines=[0,0,0,0],color=[0,90,240,140],charsize=1,/right

;if keyword_set(plotter) then paperplot,/close

;   oplot,replicate(0.685,2000),findgen(2000)
;stop

   peak  = max(why,w)
   clumpypeak = max(whyc,wc)
   shockedpeak = max(whysh,wsh)
   coldpeak = max(whycold,wcold)
   peaktime = ex[w]
   clumpypeaktime = exc[wc]
   shockedpeaktime = exsh[wsh]
   coldpeaktime = excold[wcold]

   xx = -.3

 ;  xyouts,xx,700,'peak = '+trim(peaktime)+' Gyr',charsize=1,/data
 ;     xyouts,xx,500,'peak = '+trim(clumpypeaktime)+' Gyr',charsize=1,/data,color=140
 ;  xyouts,xx,300,'peak = '+trim(shockedpeaktime)+' Gyr',charsize=1,/data,color=240
 ;  xyouts,xx,100,'peak = '+trim(coldpeaktime)+' Gyr',charsize=1,/data,color=80

;!p.multi=[0,1,2]
!p.multi=[0,1,1]
;if keyword_set(plotter) then paperplot2,filename='hz3.deltatnormalized.ps'
if keyword_set(plotter) then paperplot,filename='hz2_deltatnormalized.ps'

plothist,deltat[where(deltat ge tff)]*tunit,ex,why,xtit='dt',ytit='Normalized dN/dt',peak=1,bin=bin,title='hz2'
;xtit=textoidl( '\Delta')+'t',ytit='Normalized dN/dt',peak=1,bin=bin,title='hz2'
  plothist,deltat[c22[where(deltat[c22] ge tff)]]*tunit,exc,whyc,color=140,/over,bin=bin,peak=1
  plothist,deltat[sh22[where(deltat[sh22] ge tff)]]*tunit,exsh,whysh,color=240,/over,bin=bin,peak=1
  plothist,deltat[cold22[where(deltat[cold22] ge tff)]]*tunit,excold,whycold,color=80,/over,bin=bin,peak=1

legend,['total','cold','shocked','clumpy'],lines=[0,0,0,0],color=[0,90,240,140],charsize=1.5,/right

;if keyword_set(plotter) then paperplot,filename='h239.deltatovert.ps'

if keyword_set(plotter) then paperplot,/close
;stop

;if keyword_set(plotter) then paperplot,filename='hz3.deltatovert.ps'

plothist,(deltat[where(deltat ge tff)]*tunit)/haloaccretedtime[where(deltat ge tff)],ex,why,xtit='dt/t',ytit='N',peak=1,bin=bin,xra=[0,2] 
;xtit=textoidl( '\Delta')+'t/t',ytit='N',peak=1,bin=bin,xra=[0,2]
 plothist,(deltat[c22[where(deltat[c22] ge tff)]]*tunit)/haloaccretedtime[c22[where(deltat[c22] ge tff)]],$
        color=140,/over,bin=bin,peak=1
 plothist,(deltat[sh22[where(deltat[sh22] ge tff)]]*tunit)/haloaccretedtime[sh22[where(deltat[sh22] ge tff)]],$
        color=240,/over,bin=bin,peak=1
 plothist,(deltat[cold22[where(deltat[cold22] ge tff)]]*tunit)/haloaccretedtime[cold22[where(deltat[cold22] $
	ge tff)]],color=80,/over, bin=bin,peak=1

;legend,['total','cold','shocked','clumpy'],lines=[0,0,0,0],color=[0,90,240,140],charsize=1,/right


;if keyword_set(plotter) then paperplot2,/close

;stop
!p.multi=[0,1,1]
;if keyword_set(plotter) then paperplot,filename='hz3.deltatovert.ps'

plothist,(deltat[where(deltat ge tff)]*tunit)/haloaccretedtime[where(deltat ge tff)],ex,why,xtit='dt/t',ytit='N',peak=1,bin=bin,xra=[0,2]
 plothist,(deltat[c22[where(deltat[c22] ge tff)]]*tunit)/haloaccretedtime[c22[where(deltat[c22] ge tff)]],$
        color=140,/over,bin=bin,peak=1
 plothist,(deltat[sh22[where(deltat[sh22] ge tff)]]*tunit)/haloaccretedtime[sh22[where(deltat[sh22] ge tff)]],$
        color=240,/over,bin=bin,peak=1
 plothist,(deltat[cold22[where(deltat[cold22] ge tff)]]*tunit)/haloaccretedtime[cold22[where(deltat[cold22] $
        ge tff)]],color=80,/over, bin=bin,peak=1

legend,['total','cold','shocked','clumpy'],lines=[0,0,0,0],color=[0,90,240,140],charsize=1,/right

;if keyword_set(plotter) then paperplot,/close




;%%%%%  delta t vs dynamical time %%%%%%
readcol,'../dynamicaltime.info',rdyn,mdyn,tdyn,zdyn,format='f',/silent  ; tdyn is in YEARS
ndyn = n_elements(tdyn)
notearlyzs = zs[where(zs lt 49.,nnez)]
;stop
; match tdyn with each deltat PROBLEM  FOR z=50 entries!!  early!! bad!
match_multi,zdyn,notearlyzs,matchz,hist=hist,reverse_indices=ri  ; zs[matchz]
dynamicaltime = fltarr(nnez)
for i=0,nnez-1 do begin $
    mintime = findmin(notearlyzs[i],zdyn)
    dynamicaltime[i] = tdyn[mintime[0]]
endfor
;plot,dynamicaltime,deltat*tunit,psym=3,xtit='Dynamical Time (years)',/xlog,ytit=textoidl( '\Delta')+'t',xra=[6e8,2e9],xstyle=1,yra=[0,.5]

if keyword_set(plotter) then paperplot,filename='dynamicaltime.ps'

xbinsize=.4
ybinsize=.1
contour_plus,(dynamicaltime/1d7),deltat*tunit,xra=[1,9],xtit=' Dynamical Time (10^7 years)', xbin=xbinsize,ybin=ybinsize,levels=[1,3,4,6,10,20,50,100,200,500],yra=[-0.5,1.5],xstyle=1,ytit='dt'

meansh = fltarr(ndyn)
meanc = fltarr(ndyn)
meancl = fltarr(ndyn)
meanall = fltarr(ndyn)
; need to do some matching to make sure things are OK
notearlyiords = allgas.iords[notearlyallgas]
match,notearlyiords,clumpy,clone,cltwo
match,notearlyiords,shocked,sone,stwo
match,notearlyiords,cold,cone,ctwo
for j=0,ndyn-1 do begin 
    tsh = where(dynamicaltime[sone] eq tdyn[j],ntsh)
    if ntsh ne 0 then meansh[j] = mean(deltat[sone[tsh]])*tunit
    tc = where(dynamicaltime[cone] eq tdyn[j],ntc)
    if ntc ne 0 then meanc[j] = mean(deltat[cone[tc]])*tunit
    tcl = where(dynamicaltime[clone] eq tdyn[j],ntcl)
    if ntcl ne 0 then meancl[j] = mean(deltat[clone[tcl]])*tunit
    tall = where(dynamicaltime eq tdyn[j],ntall)
    if ntall ne 0 then meanall[j] = mean(deltat[tall])*tunit
     ; stop
endfor

shnz = where(meansh ne 0.0)
cnz = where(meanc ne 0.0)
clnz = where(meancl ne 0.0)
anz = where(meanall ne 0.0)

 oplot,(tdyn[shnz])/1d7,meansh[shnz],psym=2;,color=240
 oplot,(tdyn[cnz])/1d7,meanc[cnz],psym=1;,color=80
 oplot,(tdyn[clnz])/1d7,meancl[clnz],psym=4;,color=125
 oplot,(tdyn[anz])/1d7,meanall[anz],psym=5

legend,['cold','shocked','clumpy','total'],psym=[1,2,4,5],/right,charsize=1

if keyword_set(plotter) then paperplot,/close


stop


cumdeltat = total(deltat[where(deltat ge tff)]*tunit,/cum)/total(deltat[where(deltat ge tff)]*tunit)
cumclumpydeltat = total(deltat[c22[where(deltat[c22] ge tff)]]*tunit,/cum)/$
	total(deltat[c22[where(deltat[c22] ge tff)]]*tunit)
cumshockeddeltat = total(deltat[sh22[where(deltat[sh22] ge tff)]]*tunit,/cum)/$
	total(deltat[sh22[where(deltat[sh22] ge tff)]]*tunit)
cumcolddeltat = total(deltat[cold22[where(deltat[cold22] ge tff)]]*tunit,/cum)/$
	total(deltat[cold22[where(deltat[cold22] ge tff)]]*tunit)

plot,haloaccretedtime[sort(haloaccretedtime)],cumdeltat,xtit='halo accreted time (Gyr)',ytit='cumulative delta t'
  oplot,clumpyaccretetimes[sort(clumpyaccretetimes)],cumclumpydeltat,color=140
  oplot,shockedaccretetimes[sort(shockedaccretetimes)],cumshockeddeltat,color=240
  oplot,coldaccretetimes[sort(coldaccretetimes)],cumcolddeltat,color=80


;!p.multi=[0,1,2]
;plot,firstaccretetime*tunit,deltat*tunit,xtit='time first accreted',ytit='delta t',psym=3
;  oplot,firstaccretetime[cold1]*tunit,deltat[cold1]*tunit,psym=3,color=80
;  oplot,firstaccretetime[sh1]*tunit,deltat[sh1]*tunit,psym=3,color=250
;  oplot,firstaccretetime[c1]*tunit,deltat[c1]*tunit,psym=3,color=140
;  oplot,firstaccretetime*tunit,replicate(0.,3000),linestyle=2

;plot,haloaccretedtime*tunit,deltat*tunit,xtit='time accreted onto halo ',ytit='delta t',psym=3
;  oplot,haloaccretedtime[cold1]*tunit,deltat[cold1]*tunit,psym=3,color=80
;  oplot,haloaccretedtime[sh1]*tunit,deltat[sh1]*tunit,psym=3,color=250
;  oplot,haloaccretedtime[c1]*tunit,deltat[c1]*tunit,psym=3,color=140
;  oplot,haloaccretedtime*tunit,replicate(0.,3000),linestyle=2


;if keyword_set(plotter) then psplot3,/close

stop



end
