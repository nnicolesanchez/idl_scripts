;# This script is a redux of Jillian's massgastime.pro (find
;# in Vanderbilt VPAC40 /home/sanchenn/IDL/massgastime.pro
;# which works for low res, Gasoline simulation outputs)
;# This script works for high res, CHANGA files instead.

;# What we need:
;#      BHIORDS           (ids of bh)
;#      EATTIME/TIME_Z    (accretion time)
;#      GASIORDS          (ids of gas particles) 
;#      DM/DMQ/EATENMASS  (mass accreted from gas particles)

;# Our steps:
;# 1. Read in .fits for all above info
;#    NOTE: This is chronological order

;# Vanderbilt Univ.  -- VPAC40:
;#                      /home/sanchenn/IDL/highres_massgastime.pro
;# N. Nicole Sanchez -- Created Oct. 2015; Last Edit April 5, 2016
pro highres_massgastime,plotter=plotter
;loadct,4
loadct,39  ;#Jillian's Colors

;# Constants
tunit=1.223128407d18/3.1556926d7/1d9
;tunit = 1.223128407d18/3.1556926d7/1d9  ; Gyr 
munit = 1.84793d16  ; M_sun

;# Read above listed quantities
bhiords   = mrdfits('centralbh_bhiords.fits',0)
gasiords  = mrdfits('centralbh_gasiords.fits',0)  
bhmass    = mrdfits('centralbh_bhmass.fits',0)    ; M_sun
eatenmass = mrdfits('centralbh_eatenmass.fits',0) ; M_sun
time_z    = mrdfits('centralbh_redshift.fits',0)  ; In chronological order
time_afac = (1/(time_z + 1))
time_gyr  = mrdfits('centralbh_timeGyr.fits',0)   ; Gyr

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
restore,'allgas.sav' ; Gas is in units of M_sun
;help,allgas
zs = allgas.zaccr
ngas = n_elements(allgas.iords)

file       = mrdfits('../BHAccLog.fits',1)
;help,file

uniq_ind = uniq(file.bhiord)
;print,'Number of unique BH indeces',n_elements(uniq_ind)
uniq_ids = file.bhiord[UNIQ(file.bhiord, SORT(file.bhiord))]
print,'Number of unique BH IDs',n_elements(uniq_ids)

;print,min(file.bhiord)
;print,max(file.bhiord)
gasuniq_ids = uniq(file.gasiord)
uniq_gas = file.gasiord[UNIQ(file.gasiord, SORT(file.gasiord))]
;print,'Number of unique gas iords,',n_elements(uniq_gas),' should not equal all of them',n_elements(file.gasiord),' since file.gasiord is across all time'

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
stop

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
print,'clumpy + cold + shocked + early times',n_elements(clumpybhtimes),n_elements(coldbhtimes),n_elements(shockedbhtimes),n_elements(earlybhtimes),' equals ',n_elements(clumpybhtimes) + n_elements(coldbhtimes) + n_elements(shockedbhtimes) + n_elements(earlybhtimes),' should equal all gas ids',n_elements(time_z_rev)

;stop



if keyword_set(plotter) then paperplot,filename='cumumassvstime.ps'
plot,time_Gyr,reverse(alltotal),xtit='Time [Gyr]',ytit='Cumulative BH Mass (M'+sunsymbol()+')',yra=[1000,1e9],/ylog,xra=[0,time_gyr[0]],xstyle=9;,color=0
oplot,time_Gyr,bhmass,linestyle=2
oplot,coldbhyears,coldtotal,color=90,linestyle=3
;oplot,earlybhtimes,earlytotal,linestyle=2,color=255
oplot,clumpybhyears,clumpytotal,color=140
oplot,shockedbhyears,shockedtotal,color=240,linestyle=1
legend,['Total BH Mass','Total Gas','Cold','Clumpy','Shocked'],lines=[2,0,0,3,0],color=[0,0,80,140,240],charsize=1,/top,/left
axis,xaxis=1,xtit='Redshift',xtickv=tickvalues,xticks=n_elements(zplotnames)-1,xtickname=zplotnames
if keyword_set(plotter) then paperplot,/close
;stop

;REDSHIFT vs mass
if keyword_set(plotter) then paperplot,filename='cumumassvsz.ps'
plot,time_z,reverse(alltotal),xtit='Redshift',ytit='Cumulative BH Mass (M'+sunsymbol()+')',yra=[1000,1e8],/ylog,xra=[14,0];,color=0
oplot,time_z,bhmass,linestyle=2
oplot,coldbhtimes,coldtotal,color=80,linestyle=3
oplot,clumpybhtimes,clumpytotal,color=140
oplot,shockedbhtimes,shockedtotal,color=240,linestyle=1
legend,['Total BH Mass','Total Gas','Cold','Clumpy','Shocked'],lines=[2,0,0,3,0],color=[0,0,80,140,240],charsize=1,/top,/left
if keyword_set(plotter) then paperplot,/close

axis,xaxis=1,xtit='Redshift',xtickv=tickvalues,xticks=n_elements(zplotnames)-1,xtickname=z
;stop


;# Create mass vs time plots: need to bin data
;# BIN TIME
nbins = 100
maxtime = max(time_z)
mintime = min(time_z)
bindivision = (maxtime-mintime)/nbins
tinds = indgen(nbins)*bindivision+mintime  ; indices for time array.
bindmtime = tinds
;# BIN MASS, use sum of mass in each bin set
bindm = fltarr(nbins)
bincoldbhmass =  fltarr(nbins)
binclumpybhmass =  fltarr(nbins)
binshockedbhmass =  fltarr(nbins)
for i=0,nbins-2 do begin
    inds = where(time_z ge (i*bindivision+mintime) AND time_z lt ((i+1)*bindivision)+mintime,ninds)
    if ninds ne 0 then bindm[i] = total(eatenmass[inds])
    coldinds = where(coldbhtimes ge (i*bindivision+mintime) AND coldbhtimes lt ((i+1)*bindivision)+mintime,ncoldinds)
    if ncoldinds ne 0 then bincoldbhmass[i] = total(coldbhmass[coldinds])
    clumpyinds = where(clumpybhtimes ge (i*bindivision+mintime) AND clumpybhtimes lt ((i+1)*bindivision)+mintime,nclumpyinds)
    if nclumpyinds ne 0 then binclumpybhmass[i] = total(clumpybhmass[clumpyinds])
    shockedinds = where(shockedbhtimes ge (i*bindivision+mintime) AND shockedbhtimes lt ((i+1)*bindivision)+mintime,nshockedinds)
    if nshockedinds ne 0 then binshockedbhmass[i] = total(shockedbhmass[shockedinds])
endfor
stop
;# Mass vs Redshift
if keyword_set(plotter) then paperplot,filename='massvstime.ps'
;print,bindmtime
plot,bindmtime,bindm,xtit='Redshift',ytit='Mass Accreted by Central BH',/ylog,yra=[1,1e6],xra=[20,0]
oplot,bindmtime,binclumpybhmass,color=140;,psym=3
oplot,bindmtime,bincoldbhmass,color=90   ;,psym=3
oplot,bindmtime,binshockedbhmass,color=240 ;,psym=3
legend,['all','cold','shocked','clumpy'],lines=[0,0,0,0],color=[0,90,240,140],charsize=1,/right,/bottom
if keyword_set(plotter) then paperplot,/close
stop

;# Cumulative Mass Plot
; THIS PLOT NEEDS FIXING
;# Need to convert all times to a factor to multiply by tunit
;if keyword_set(plotter) then paperplot,filename='allmassgas.ps'
;plot,time_z,reverse(alltotal),xtit='Redshift',ytit='Cumulative BH Mass (M'+sunsymbol()+')',/ylog,yra=[100,1e8],xra=[20,0];,xstyle=8  ;puts on upper x-axis
;oplot,time_z,clumpytotal,color=125
;oplot,time_z,bhmass,lines=2,color=0
;oplot,time_z,coldtotal,lines=0,color=80
;oplot,time_z,shockedtotal,lines=0,color=240
;oplot,time_z,earlytotal,color=5

;legend,['total BH mass','total gas','cold','shocked','clumpy'],lines=[2,0,0,0,0],color=[0,0,80,240,125],charsize=1,/top,/left
;# Fix following line for upper x-axis
;axis,xaxis=1,xtit='Redshift',xtickv=tickvalues,xticks=n_elements(zplotnames)-1,xtickname=zplotnames
;if keyword_set(plotter) then paperplot,/close
;stop

; ### test that the parts add up to the whole. ###

nclumpy = n_elements(clumpytotal)
ncold = n_elements(coldtotal)
nshocked = n_elements(shockedtotal)
ndm = n_elements(alltotal)
nearly = n_elements(earlytotal)

print,'Does the final mass of alltotal',alltotal[ndm-1],' equal the final masses of all the others?',clumpytotal[nclumpy-1]+coldtotal[ncold-1]+shockedtotal[nshocked-1]+earlytotal[nearly-1]

print,'gas mass fractions '
print,'cold/total = ',coldtotal[ncold-1]/alltotal[ndm-1]
print,'clumpy/total = ',clumpytotal[nclumpy-1]/alltotal[ndm-1]
print,'shocked/total = ',shockedtotal[nshocked-1]/alltotal[ndm-1]
print,'early/total = ',earlytotal[nearly-1]/alltotal[ndm-1]

; These bottom three things need fixing because all arrays have been made
; cumulative arrays and therefore total()s can NOT be taken of them
; and be representative.

;print,'gas mass fractions (weird?)'
;print,'cold/total = ',total(coldtotal)/total(alltotal)
;print,'clumpy/total = ',total(clumpytotal)/total(alltotal)
;print,'shocked/total = ',total(shockedtotal)/total(alltotal)
;print,'early/total = ',total(earlytotal)/total(alltotal)

;print,'does the total mass of alltotal',alltotal[ndm-1],' equal the added sums of all the others?','cold',total(coldtotal),' clumpy',total(clumpytotal),' shocked',total(shockedtotal),' early',total(earlytotal)
;print,'total =',coldtotal[ncold-1]+clumpytotal[nclumpy-1]+shockedtotal[nshocked-1]+earlytotal[nearly-1]

;print,'gas number fractions '
;print,'cold/total = ',n_elements(coldbhmass)
;/n_elements(alltotal)
;print,'clumpy/total = ',n_elements(clumpybhmass)
;/n_elements(alltotal)
;print,'shocked/total = ',n_elements(shockedbhmass)
;/n_elements(alltotal)
;print,'early/total = ',n_elements(earlybhmass)/n_elements(alltotal)

end
