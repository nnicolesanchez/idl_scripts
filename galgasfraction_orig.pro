pro galgasfraction,plotter=plotter,plottitle=plottitle,dolegend=dolegend

; determine the mass fraction of gas in halo 1 vs time that is clumpy,
; cold, shocked.

; plottitle = if I want a title on the top of my plot.  "h258" or
; whatever.
; dolegend = do I want a legend?

; units
lunit = 50000. ; kpc
munit = 1.84793d16 ; msun
grey = cgcolor('Dark Gray')

if not keyword_set(plottitle) then plottitle=''

; files
haloidfile = file_search('*haloid.dat')
readcol,haloidfile[0],file,haloid,format='a,i',/silent
file=reverse(file)
nfiles=n_elements(file)
haloid=reverse(haloid)          ; these are backwards.
clumpy = mrdfits('clumpy.accr.iord.fits',0)
cold = mrdfits('unshock.iord.fits',0)
shocked = mrdfits('shocked.iord.fits',0)
accrz = mrdfits('grp1.accrz.fits',0)
redshifts = accrz[uniq(accrz,sort(accrz))]
; only analyze to z=3.89 and where tracing begins.
redshifts = redshifts[where(redshifts ge 3.85 AND redshifts le 16.)]
redshifts = [50,reverse(redshifts)]
grp1 = mrdfits('grp1.allgas.iord.fits',0)
early = mrdfits('early.iord.fits',0)
; sanity check.
print,n_elements(redshifts),' must equal',n_elements(file)

; has this already been run?  if so, read in the table.  if not, make
; a table.
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

   for i=0,nfiles-1 do begin
      rtipsy,file[i],h,g,d,s
      d=''
      s=''
      iord = read_lon_array(file[i]+'.iord')
      iord= iord[0:h.ngas-1]
      ingal = where(accrz ge redshifts[i],ningal)
      print,ningal,' in gal'

      match,clumpy,notearly[ingal],cl1,cl2
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
;    stat = mrdfits(file[i]+'.amigastat.fits',1)
;stop
endfor

openw,lun,'gasfraction.dat',/get_lun,width=500
printf,lun,'cold mass  ','clumpy mass  ','shocked mass   ','early mass   '
for i=0,nfiles-1 do printf,lun,coldmass[i],clumpymass[i],shockedmass[i],earlymass[i]
close,lun
free_lun,lun

endif else readcol,'gasfraction.dat',coldmass,clumpymass,shockedmass,earlymass,format='d',/silent


totalmass = coldmass+shockedmass+clumpymass+earlymass

; normalize to total gas mass in halo 1?
; good sanity check to see if it adds up 

if keyword_set(plotter) then paperplot,filename='gasfraction.ps'

plot,redshifts,clumpymass/totalmass,xtit='Redshift',ytit='Gas Mass Fraction',xra=[11,3],xstyle=1,yra=[0,1],linestyle=0,title=plottitle
  oplot,redshifts,clumpymass/totalmass,color=125,linestyle=0;,psym=-1
  oplot,redshifts,coldmass/totalmass,color=80,linestyle=0
  oplot,redshifts,shockedmass/totalmass,color=240,linestyle=0
  oplot,redshifts,earlymass/totalmass,color=grey,linestyle=0

if keyword_set(dolegend) then legend,['clumpy','unshocked','shocked','early'],colors=[125,80,240,grey],lines=[0,0,0,0],/right,charsize=1


if keyword_set(dolegend) then legend,['BH','galaxy'],colors=[0,0],lines=[2,0],thick=[0,4],charsize=1,position=[8.3,0.95]


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


    oplot,dmtimez,clumpyfraction,color=125,linestyle=2
    oplot,dmtimez,coldfraction,color=80,linestyle=2
    oplot,dmtimez,shockedfraction,color=240,linestyle=2
    oplot,dmtimez,earlyfraction,color=grey,linestyle=2
; testing
;plot,bhaccretemass.dmtime,bhaccretemass.alltotal,/ylog,xra=[.01,.015]
;oplot,bhaccretemass.shockedbhtimes,bhaccretemass.shockedtotal,color=240
;oplot,bhaccretemass.coldbhtimes,bhaccretemass.coldtotal,color=80
;oplot,bhaccretemass.clumpybhtimes,bhaccretemass.clumpytotal,color=80


endif

if keyword_set(plotter) then paperplot,/close

;stop




end
