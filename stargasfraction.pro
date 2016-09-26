pro stargasfraction,plotter=plotter,plottitle=plottitle,dolegend=dolegend

; determine the mass fraction of stars in halo 1 vs time that came from clumpy,
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
haloid=reverse(haloid) ; these are backwards.
readcol,'times.list',time,redshifts,format='f',/silent
; sanity check.
   print,n_elements(redshifts),' must equal',n_elements(file)
   if n_elements(redshifts) ne n_elements(file) then begin
      difference = n_elements(redshifts) - n_elements(file)
      redshifts = redshifts[difference:n_elements(redshifts)-1]
      print,n_elements(redshifts),' must equal',n_elements(file)
   endif

istherefile = file_search('starfraction.dat')
if istherefile eq '' then begin

   clumpy = mrdfits('clumpy.accr.iord.fits',0)
   cold = mrdfits('unshock.iord.fits',0)
   shocked = mrdfits('shocked.iord.fits',0)
;   accrz = mrdfits('grp1.accrz.fits',0)
;   redshifts = accrz[uniq(accrz,sort(accrz))]
; only analyze to z=3.89 and where tracing begins.
;   redshifts = redshifts[where(redshifts ge 3.85 AND redshifts le 16.)]
;   redshifts = [50,reverse(redshifts)]
   early = mrdfits('early.iord.fits',0)


   coldstarmass = fltarr(nfiles)
   shockedstarmass = fltarr(nfiles)
   clumpystarmass = fltarr(nfiles)
   earlystarmass = fltarr(nfiles)
   coldBHmass = fltarr(nfiles)
   shockedBHmass = fltarr(nfiles)
   clumpyBHmass = fltarr(nfiles)
   earlyBHmass = fltarr(nfiles)

   for i=0,nfiles-1 do begin
      rtipsy,file[i],h,g,d,s
      d=''
      g=''
      gasiord= read_lon_array(file[i]+'.igasorder')
      gasiord = gasiord[h.ngas+h.ndark:h.n-1]
      grp =  read_lon_array(file[i]+'.amiga.grp')
      grp = grp[h.ngas+h.ndark:h.n-1]

      ingal = where(grp eq haloid[i],ningal)
      print,ningal,' in gal'

      match,clumpy,gasiord,cl1,cl2
      match,shocked,gasiord,s1,s2
      match,cold,gasiord,c1,c2
      match,early,gasiord,e1,e2

      clumpystarmass[i] = total(s[cl2].mass)*munit
      coldstarmass[i] = total(s[c2].mass)*munit
      shockedstarmass[i] = total(s[s2].mass)*munit
      earlystarmass[i] = total(s[e2].mass)*munit
      
   endfor

   openw,lun,'starfraction.dat',/get_lun,width=500
   printf,lun,'total STELLAR MASSES'
   printf,lun,'cold mass  ','clumpy mass  ','shocked mass   ','early mass   '
   for i=0,nfiles-1 do printf,lun,coldstarmass[i],clumpystarmass[i],shockedstarmass[i],earlystarmass[i]
   close,lun
   free_lun,lun

endif else  readcol,istherefile[0],coldstarmass,clumpystarmass,shockedstarmass,earlystarmass,format='d',/silent

totalstarmass = coldstarmass+shockedstarmass+clumpystarmass+earlystarmass
nn=n_elements(coldstarmass)

; print the fractions
print,'cold fraction= ',coldstarmass[nn-1]/totalstarmass[nn-1]
print,'clumpy fraction= ',clumpystarmass[nn-1]/totalstarmass[nn-1]
print,'hot fraction= ',shockedstarmass[nn-1]/totalstarmass[nn-1]
print,'early fraction= ',earlystarmass[nn-1]/totalstarmass[nn-1]



; normalize to total gas mass in halo 1?
; good sanity check to see if it adds up 

if keyword_set(plotter) then paperplot,filename='stargasfraction.ps'

plot,redshifts,clumpystarmass/totalstarmass,xtit='Redshift',ytit='Cumulative Stellar Mass Fraction',xra=[10,3],xstyle=1,yra=[0,1],title=plottitle
  oplot,redshifts,clumpystarmass/totalstarmass,color=125
  oplot,redshifts,coldstarmass/totalstarmass,color=80
  oplot,redshifts,shockedstarmass/totalstarmass,color=240
  oplot,redshifts,earlystarmass/totalstarmass,color=grey

if keyword_set(dolegend) then legend,['clumpy','unshocked','shocked','early'],colors=[125,80,240,grey],lines=[0,0,0,0],/right,charsize=1

; now.  I want to put the fraction of gas accreted by the
; BH here as well. 


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

endif

if keyword_set(plotter) then paperplot,/close

;stop




end
