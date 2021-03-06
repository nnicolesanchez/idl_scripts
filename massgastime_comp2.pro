pro massgastime_comp2,plotter=plotter


; LOW RES h258 PART OF PLOT
restore,'../lowresh258/bhorbits.sav'
tunit = 1.223128407d18/3.1556926d7/1d9 ; Gyr
munit = 1.84793d16
restore,'../lowresh258/centralbh/allgas.sav'
readcol,'../lowresh258/centralbh/out.distance',w,bhiords,w,eattime,w,gasiords,dsquare,w,w,smooth,w,eatenmass,format='a,l,a,f,a,l,f,a,a,f,a,f',/silent
readcol,'../lowresh258/centralbh/out.dm',w,bhiords,w,delta,w,dmtime,w,dm,format='a,l,a,f,a,f,a,f',/silent
loadct,39  ;#Jillian's Colors
; the rest of the columns do not interest me.


; praise the gods, dmtime and eattime are equal arrays.
; this means that dm and gasiords correspond one to one
; YAY
; match dm with gasiords and I'll have clumpy and smooth.

; times and redshift, for the z axis on the plots
readcol,'../lowresh258/times.list',gyr,z,/silent
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
coldtotal=total(coldbhmass,/cumulative)
shockedtotal=total(shockedbhmass,/cumulative)
clumpytotal=total(clumpybhmass,/cumulative)
alltotal = total(eatenmass,/cumulative)
earlytotal= total(earlybhmass,/cumulative)

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

if keyword_set(plotter) then paperplot,filename='allmassgas.ps'

;### cumulative mass plot ###

; use orbit file to get max BH mass
;bhorbitmass = fltarr(n_elements(bhorbits[1].x))
;for i=0,n_elements(bhorbits[35].x)-1 do bhorbitmass[i] =
;max(bhorbits.mass[i])*munit
;#### for h258: bhorbits[39] &&& for h277: bhorbits[223] (represents
;index of bh with largest mass in final timestep)
;######### [0] for up to redshift = 2 BUT [223] for full range
bhorbitmass = bhorbits[0].mass*munit
bhorbittime = bhorbits[0].time*tunit

bhorbittime=bhorbittime[where( bhorbittime ne 0)]
bhorbitmass=bhorbitmass[where( bhorbitmass ne 0)]

plot,dmtime*tunit,alltotal,xtit='Time (Gyr)',ytit='Cumulative BH Mass (M'+sunsymbol()+')',/ylog,yra=[100,1e8],xra=[0,14],xstyle=8,linestyle=2 ;xstyle=1 forces range to my range
oplot,dmtime*tunit,alltotal*munit,color=240
legend,['Gasoline total BH mass','ChaNGa total BH mass'],lines=[0,3],color=[240,90],charsize=1,/right,/bottom
axis,xaxis=1,xtit='Redshift',xtickv=tickvalues,xticks=n_elements(zplotnames)-1,xtickname=zplotnames
stop



; HIGH RES h258 PART OF PLOT
;# Constants
tunit=1.223128407d18/3.1556926d7/1d9
;tunit = 1.223128407d18/3.1556926d7/1d9  ; Gyr 
munit = 1.84793d16  ; M_sun

;# Read above listed quantities
bhiords   = mrdfits('../highresh258/centralbh/centralbh_bhiords.fits',0)
gasiords  = mrdfits('../highresh258/centralbh/centralbh_gasiords.fits',0)  
bhmass    = mrdfits('../highresh258/centralbh/centralbh_bhmass.fits',0)    ; M_sun
eatenmass = mrdfits('../highresh258/centralbh/centralbh_eatenmass.fits',0) ; M_sun
time_z    = mrdfits('../highresh258/centralbh/centralbh_redshift.fits',0)  ; In chronological order
time_afac = (1/(time_z + 1))
time_gyr  = mrdfits('../highresh258/centralbh/centralbh_timeGyr.fits',0)   ; Gyr

;# Calculate time in Gyr
; times and redshift, for the z axis on the plots
;readcol,'../times.list',gyr,z,/silent
zplotnames=['6','4','2','1','0.5','0']
zplotnum=[6,4,2,1,0.5,0]
tickvalues=interpol(time_gyr,time_z,zplotnum) 
tickvalues[n_elements(tickvalues)-1]=time_gyr[n_elements(time_gyr)-1]
;print,time_gyr[0]
;stop

;# Read in times to create the x-axis for our plots later
;# Following is from Jillian's massgastime.pro
zplotnames=['6','4','2','1','0.5','0']
zplotnum=[6,4,2,1,0.5,0]
tickvalues=interpol(time_gyr,time_z,zplotnum)
tickvalues[n_elements(tickvalues)-1]=time_gyr[n_elements(gyr)-1]

;# Read in my .fits files (created by highres_gasiords.pro)
restore,'../highresh258/centralbh/allgas.sav' ; Gas is in units of M_sun
;help,allgas
zs = allgas.zaccr
ngas = n_elements(allgas.iords)

file       = mrdfits('../highresh258/BHAccLog.fits',1)
;help,file

uniq_ind = uniq(file.bhiord)
;print,'Number of unique BH indeces',n_elements(uniq_ind)
uniq_ids = file.bhiord[UNIQ(file.bhiord, SORT(file.bhiord))]
print,'Number of unique BH IDs',n_elements(uniq_ids)

;print,min(file.bhiord)
;print,max(file.bhiord)
gasuniq_ids = uniq(file.gasiord)
uniq_gas = file.gasiord[UNIQ(file.gasiord, SORT(file.gasiord))]
;print,'Number of unique gas iords,',n_elements(uniq_gas),'
;should not equal all of them',n_elements(file.gasiord),' since
;file.gasiord is across all time'

uniq_gasiords = gasiords[UNIQ(gasiords, SORT(gasiords))]
print,'Number of unique gas IDs in gasiords',uniq_gasiords[0:100]


match_multi,allgas.iords,gasiords,all1
print,'Number of gas IDS in gasiords',n_elements(gasiords)
; We want to use gasiords because it directly compares to eatenmass
; and redshift and we know the BH can accrete from a gas particle
; multiple times by taking "bites" out of it
print,'Number of unique gas IDs in allgas.sav',ngas
print,'Number of unique gas IDs in BHAccLog',n_elements(uniq_gas)
; This number SHOULD be larger than the uniq gas Ids value for both
; allgas.sav and gasiords since they both only consider halo 1 and
; BHAccLog includes ALL BHs
print,'Number of matched gas IDs between allgas.sav & BHAccLog.fits',n_elements(all1)
print,'Total BH accreted particles',n_elements(gasiords)
print,'BH accreted categorized  particles',n_elements(allgas.iords)

;stop

;# Divide Gas into Smooth and Clumpy
;# CLUMPY:
match_multi,allgas.clumpy,gasiords,c11
clumpybhtimes = time_z[c11]
clumpybhyears = time_Gyr[c11]
clumpybhmass  = eatenmass[c11]
;# SHOCKED
match_multi,allgas.shocked,gasiords,sh11
shockedbhtimes = time_z[sh11]
shockedbhyears = time_Gyr[sh11]
shockedbhmass  = eatenmass[sh11]
;# COLD
match_multi,allgas.cold,gasiords,cold11
coldbhtimes = time_z[cold11]
coldbhyears = time_Gyr[cold11]
coldbhmass  = eatenmass[cold11]
;# EARLY
match_multi,allgas.early,gasiords,early1
earlybhtimes = time_z[early1]
earlybhyears = time_Gyr[early1]
earlybhmass = eatenmass[early1]
print,'How many early gas in allgas?',n_elements(allgas.early)
print,'How many shocked gas in allgas?',n_elements(allgas.shocked)
print,'How many clumpy gas in allgas?',n_elements(allgas.clumpy)
print,'How many cold gas in allgas?',n_elements(allgas.cold)
;stop

;# Put back in chronological order
;# CLUMPY
sortclumpy    = sort(clumpybhtimes)
clumpybhtimes = clumpybhtimes[sortclumpy]
clumpybhyears = clumpybhyears[sortclumpy]
clumpybhmass  = clumpybhmass[sortclumpy]
clumpybhtimes = reverse(clumpybhtimes)
clumpybhyears = reverse(clumpybhyears)
clumpybhmass  = reverse(clumpybhmass)
;# SHOCKED
sortshocked    = sort(shockedbhtimes)
shockedbhtimes = shockedbhtimes[sortshocked]
shockedbhyears = shockedbhyears[sortshocked]
shockedbhmass  = shockedbhmass[sortshocked]
shockedbhtimes = reverse(shockedbhtimes)
shockedbhyears = reverse(shockedbhyears)
shockedbhmass  = reverse(shockedbhmass)
;# COLD
sortcold    = sort(coldbhtimes)
coldbhtimes = coldbhtimes[sortcold]
coldbhyears = coldbhyears[sortcold]
coldbhmass  = coldbhmass[sortcold]
coldbhtimes = reverse(coldbhtimes)
coldbhyears = reverse(coldbhyears)
coldbhmass  = reverse(coldbhmass)

;# EARLY
sortearly    = sort(earlybhtimes)
earlybhtimes = earlybhtimes[sortearly]
earlybhyears = earlybhyears[sortearly]
earlybhmass  = earlybhmass[sortearly]
earlybhtimes = reverse(earlybhtimes)
earlybhyears = reverse(earlybhyears)
earlybhmass  = reverse(earlybhmass)

;# TOTAL
eatenmass_rev = reverse(eatenmass)
time_z_rev    = reverse(time_z)
time_Gyr_rev    = reverse(time_Gyr)
;# BH TOTAL
;bhmass_rev    = reverse(bhmass)
bhmass_rev    = bhmass

;# Add and get cumulative totals
coldtotal    = total(coldbhmass,/cumulative)
shockedtotal = total(shockedbhmass,/cumulative)
clumpytotal  = total(clumpybhmass,/cumulative)
earlytotal   = total(earlybhmass,/cumulative)
alltotal     = total(eatenmass_rev,/cumulative)
bhtotal      = total(bhmass_rev,/cumulative)



;print,clumpybhtimes[0:50]
;print,shockedbhtimes[0:50]
;print,coldbhtimes[0:50]
;print,earlybhtimes
;print,time_z_rev[0:50]


eatenmasstotal = total(eatenmass,/cumulative)
print,'clumpy + cold + shocked + early ids',n_elements(clumpytotal),n_elements(coldtotal),n_elements(shockedtotal),n_elements(earlytotal),' equals ',n_elements(clumpytotal) + n_elements(coldtotal) + n_elements(shockedtotal) + n_elements(earlytotal),' should equal all gas ids',n_elements(gasiords)
print,'clumpy + cold + shocked + early times',n_elements(clumpybhtimes),n_elements(coldbhtimes),n_elements(shockedbhtimes),n_elements(earlybhtimes),' equals ',n_elements(clumpybhtimes) + n_elements(coldbhtimes) + n_elements(shockedbhtimes) +  n_elements(earlybhtimes),' should equal all gas ids',n_elements(time_z_rev)

;stop

;if keyword_set(plotter) then paperplot,filename='cumumassvstime.ps'
oplot,time_Gyr,reverse(alltotal),linestyle=3,color=90;,xtit='Time [Gyr]',ytit='Cumulative BH Mass (M'+sunsymbol()+')',yra=[1000,1e9],/ylog,xra=[0,time_gyr[0]],xstyle=9;,color=0
;oplot,time_Gyr,bhmass,linestyle=2
;oplot,coldbhyears,coldtotal,color=90,linestyle=3
;oplot,earlybhtimes,earlytotal,linestyle=2,color=255
;oplot,clumpybhyears,clumpytotal,color=140
;oplot,shockedbhyears,shockedtotal,color=240,linestyle=1
;legend,['Total BH Mass','Total
;Gas','Cold','Clumpy','Shocked'],lines=[2,0,0,3,0]
;,color=[0,0,80,140,240],charsize=1,/top,/left
;axis,xaxis=1,xtit='Redshift',xtickv=tickvalues,xticks=n_elements(zplotnames)-1,x
;tickname=zplotnames

if keyword_set(plotter) then paperplot,/close







end
