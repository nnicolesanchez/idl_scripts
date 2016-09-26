pro galnumfraction_changa,plotter=plotter,plottitle=plottitle,dolegend=dolegend

; determine the mass fraction of gas in halo 1 vs time that is clumpy,
; cold, shocked.

; plottitle = if I want a title on the top of my plot.  "h258" or
; whatever.
; dolegend = do I want a legend?

; units
lunit = 50000. ; kpc
munit = 1.84793d16 ; msun
loadct,39
;loadct,4    ; #For Nicole's colors
grey = cgcolor('Dark Gray')
if not keyword_set(plottitle) then plottitle=''

; files
haloidfile = file_search('*haloid.dat')
readcol,haloidfile[0],file,haloid,format='a,i',/silent
file='../highresh258'+reverse(file)
nfiles=n_elements(file)
print,nfiles
haloid=reverse(haloid)          ; these are backwards.
clumpy = mrdfits('clumpy.accr.iord.fits',0)
cold = mrdfits('unshock.iord.fits',0)
shocked = mrdfits('shocked.iord.fits',0)
;accrz = mrdfits('early.iord.fits',0)
accrz = mrdfits('grp1.accrz.fits',0)

redshifts = accrz[uniq(accrz,sort(accrz))]
;print,n_elements(redshifts)
 
; only analyze to z=0 and where tracing begins.
redshifts = redshifts[where (redshifts ge 0. AND redshifts le 4.5)]
redshifts = [4.5,reverse(redshifts)]
grp1 = mrdfits('grp1.allgas.iord.fits',0)
early = mrdfits('early.iord.fits',0)
; sanity check.
print,n_elements(redshifts),' must equal',n_elements(file)
;stop
; has this already been run?  if so, read in the table.  if not, make
; a table.

;JAN.22.2015: currently off to test what is wrong with h277
gasfractiontable = file_search('gasfraction.dat')
if gasfractiontable[0] eq '' then begin

;stop
; eliminate early to get accrz for gas I care about
   notearly = grp1
   match,early,notearly,e1,e2
   remove,e2,notearly

   coldmass = fltarr(nfiles)
   shockedmass = fltarr(nfiles)
   clumpymass = fltarr(nfiles)
   earlymass = fltarr(nfiles)
   coldfrac = fltarr(nfiles)
   shockedfrac = fltarr(nfiles)
   clumpyfrac = fltarr(nfiles)
   earlyfrac = fltarr(nfiles)
   coldmassfrac = fltarr(nfiles)
   shockedmassfrac = fltarr(nfiles)
   clumpymassfrac = fltarr(nfiles)
   earlymassfrac = fltarr(nfiles)


   for i=1,nfiles-1 do begin
      rtipsy,file[i],h,g,d,s
      d=''
      s=''
      iord = read_lon_array(file[i]+'.iord')
      iord= iord[0:h.ngas-1]
      ingal = where(accrz ge redshifts[i],ningal)
      print,ningal,' in gal'

      match,clumpy,notearly[ingal],cl1,cl2
 ;        print,n_elements(cl2)/ningal
 ;        print,cl2
      match,shocked,notearly[ingal],s1,s2
      match,cold,notearly[ingal],c1,c2
    ; ID with iords and get masses
      match,clumpy[cl1],iord,icl1,icl2
      match,shocked[s1],iord,is1,is2
      match,cold[c1],iord,ic1,ic2
      match,early,iord,e1,e2
      coldmass[i] = total(g[ic2].mass)*munit
      shockedmass[i] = total(g[is2].mass)*munit
      clumpymass[i] = total(g[icl2].mass)*munit
      earlymass[i] = total(g[e1].mass)*munit
      totalmass=coldmass[i]+shockedmass[i]+clumpymass[i]+earlymass[i]
      
      coldmassfrac[i]=coldmass[i]/totalmass
      	coldmassfrac[0]=1.0000
      clumpymassfrac[i]=clumpymass[i]/totalmass
      shockedmassfrac[i]=shockedmass[i]/totalmass
      earlymassfrac[i]=earlymass[i]/totalmass

      totaln=n_elements(c2)+n_elements(s2)+n_elements(cl2)+n_elements(e1)
      print,'Total particles in this timestep including early =',totaln
      coldfrac[i]=float(n_elements(c2))/totaln
      	coldfrac[0]=1.0000
      shockedfrac[i]=float(n_elements(s2))/totaln
      clumpyfrac[i]=float(n_elements(cl2))/totaln
      earlyfrac[i]=float(n_elements(e1))/totaln
		

print,'Does cold gas fraction',float(n_elements(c2))/totaln,' match gas fraction',coldmass[i]/totalmass
;print,shockedmass[i]/totalmass
print,'Does clumpy gas fraction',float(n_elements(cl2))/totaln,' match gas fraction', clumpymass[i]/totalmass
;print,earlymass[i]/totalmass

;    stat = mrdfits(file[i]+'.amigastat.fits',1)
;stop
      endfor
stop

openw,lun,'gasfraction.dat',/get_lun,width=500
printf,lun,'cold mass  ','clumpy mass  ','shocked mass   ','early mass   '
for i=0,nfiles-1 do printf,lun,coldmass[i],clumpymass[i],shockedmass[i],earlymass[i]
close,lun
free_lun,lun

openw,lun,'gasfraction.dat',/get_lun,width=500
printf,lun,'cold mass frac ','clumpy mass frac ','shocked mass frac  ','early mass frac  ','cold # frac   ','clumpy # frac   ','shocked # frac   ','early # frac   '
for i=0,nfiles-1 do printf,lun,coldmassfrac[i],clumpymassfrac[i],shockedmassfrac[i],earlymassfrac[i],coldfrac[i],clumpyfrac[i],shockedfrac[i],earlyfrac[i]
close,lun
free_lun,lun

;stop

endif else readcol,'gasfraction.dat',coldmass,clumpymass,shockedmass,earlymass,coldfrac,clumpyfrac,shockedfrac,earlyfrac,format='d',/silent

;totaln = n_elements(clumpyfrac)+n_elements(coldfrac)+n_elements(earlyfrac)+n_elements(shockedfrac)
totalmass=clumpymass+coldmass+shockedmass+earlymass
; normalize to total gas mass in halo 1?
; good sanity check to see if it adds up 

if keyword_set(plotter) then paperplot,filename=plottitle+'numfraction.ps'

plot,redshifts,clumpyfrac,xtit='Redshift',ytit='Gas Fraction',xra=[4.5,0],xstyle=1,yra=[0,1],linestyle=0,thick=4,title=plottitle
  oplot,redshifts,clumpyfrac,color=159,linestyle=0,thick=3;,psym=-1
  oplot,redshifts,coldfrac,color=60,linestyle=0,thick=3
;  oplot,redshifts,shockedfrac,color=170,linestyle=0,thick=3
;  oplot,redshifts,earlyfrac,color=5,linestyle=0,thick=3

;print,clumpyfrac[1]
;print,coldfrac[1]
print,clumpyfrac[43]
print,coldfrac[43]
;stop

if keyword_set(dolegend) then legend,['cold','clumpy'],colors=[60,159],lines=[0,0],thick=4,/right,charsize=1,position=[0.1,0.98]


if keyword_set(dolegend) then legend,['BH','galaxy'],colors=[0,0],lines=[0,0],thick=[8,3],charsize=1,position=[2.8,0.98]


; now.  It's likely I want to put the fraction of gas accreted by the
; BH here as well.  But only if the necessary files exist.



bhfile = file_search('centralbh/bhaccretemass.sav')
if n_elements(bhfile) ne 0 then begin
    restore,bhfile[0]
    ; these times need to be coverted to redshift, ugh.
    readcol,'centralbh/out.mdot',w,w,w,mdottime,w,w,w,w,w,w,w,a,format='a,a,a,f,a,a,a,a,a,a,a,f',/silent
    redshift = 1./a-1.
    dmtimez = interpol(redshift,mdottime,bhaccretemass.dmtime)
    ; OK now match times to each gas phase
    ; interpolate all times.
    clumpytotal = interpol(bhaccretemass.clumpytotal,bhaccretemass.clumpybhtimes,bhaccretemass.dmtime)
    coldtotal = interpol(bhaccretemass.coldtotal,bhaccretemass.coldbhtimes,bhaccretemass.dmtime)
    shockedtotal = interpol(bhaccretemass.shockedtotal,bhaccretemass.shockedbhtimes,bhaccretemass.dmtime)
    earlytotal = interpol(bhaccretemass.earlytotal,bhaccretemass.earlybhtimes,bhaccretemass.dmtime)
    ; interpolation is sucky for early at the end.  fix it.
    earlytotal[where(earlytotal ge max(bhaccretemass.earlytotal))] =  max(bhaccretemass.earlytotal)
    clumpyfraction = clumpytotal/bhaccretemass.alltotal
    coldfraction = coldtotal/bhaccretemass.alltotal
    shockedfraction = shockedtotal/bhaccretemass.alltotal
    earlyfraction = earlytotal/bhaccretemass.alltotal

;print,clumpyfraction
;print,coldfraction[260396]
;print,shockedfraction[260396]
;print,earlyfraction[260396]
;colors for Jillian: clumpy,125,cold,80,shocked,240
    oplot,dmtimez,clumpyfraction,color=159,linestyle=0,thick=8
    oplot,dmtimez,coldfraction,color=60,linestyle=0,thick=8
;    oplot,dmtimez,shockedfraction,color=170,linestyle=2
;    oplot,dmtimez,earlyfraction,color=5,linestyle=2

; testing
;plot,bhaccretemass.dmtime,bhaccretemass.alltotal,/ylog,xra=[.01,.015]
;oplot,bhaccretemass.shockedbhtimes,bhaccretemass.shockedtotal,color=240
;oplot,bhaccretemass.coldbhtimes,bhaccretemass.coldtotal,color=80
;oplot,bhaccretemass.clumpybhtimes,bhaccretemass.clumpytotal,color=80


endif

if keyword_set(plotter) then paperplot,/close

;stop




end
