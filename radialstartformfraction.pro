pro radialstartformfraction,plotter=plotter,plottitle=plottitle,dolegend=dolegend

; determine the mass fraction of stars in halo 1 vs time that came from clumpy,
; cold, shocked.
; do it by radius!

; divide it up by recently accreted BH gas and recently formed stars.
; use time bins that equal the output intervals at first.

; units
lunit = 50000. ; kpc
munit = 1.84793d16 ; msun
tunit=1.223128407d18/3.1556926d7/1d9 ; Gyr
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
      time = time[difference:n_elements(redshifts)-1]
      redshifts = redshifts[difference:n_elements(redshifts)-1]
      print,n_elements(redshifts),' must equal',n_elements(file)
   endif

istherefile = file_search('radialstartformfraction.dat')
if istherefile eq '' then begin


   clumpy = mrdfits('clumpy.accr.iord.fits',0)
   cold = mrdfits('unshock.iord.fits',0)
   shocked = mrdfits('shocked.iord.fits',0)
   accrz = mrdfits('grp1.accrz.fits',0)
   early = mrdfits('early.iord.fits',0)
   bhfile = file_search('centralbh/bhaccretemass.sav')
   restore,bhfile[0]
;   readcol,'centralbh/out.mdot',w,w,w,mdottime,w,w,w,w,w,w,w,a,format='a,a,a,f,a,a,a,a,a,a,a,f',/silent



   coldstarmass = fltarr(nfiles)
   shockedstarmass = fltarr(nfiles)
   clumpystarmass = fltarr(nfiles)
   earlystarmass = fltarr(nfiles)
   coldBHmass = fltarr(nfiles)
   shockedBHmass = fltarr(nfiles)
   clumpyBHmass = fltarr(nfiles)
   earlyBHmass = fltarr(nfiles)
; initialize time
   t0 = 0.0
   radiusofinterest = 0.5 ; physical kpc
   for i=0,nfiles-1 do begin
      rtipsy,file[i],h,g,d,s
      d=''
      g=''
      gasiord= read_lon_array(file[i]+'.igasorder')
      gasiord = gasiord[h.ngas+h.ndark:h.n-1]
      grp =  read_lon_array(file[i]+'.amiga.grp')
      grp = grp[h.ngas+h.ndark:h.n-1]
      tform = read_ascii_array(file[i]+'.timeform')
      tform = tform[h.ngas+h.ndark:h.n-1]
      readcol,file[i]+'.shrinkcenters',haloidc,xx,yy,zz,format='i,d,d,d',/silent
      location = where(haloidc eq haloid[i])
      xc = xx[location[0]]
      yc = yy[location[0]]
      zc = zz[location[0]]
      ; recenter stars
      sx = s.x-xc
      sy = s.y-yc
      sz = s.z-zc
      radius = sqrt(sx*sx+sy*sy+sz*sz)*h.time*lunit ; physical kpc

    ; find stars in tform range
      ingal = where(grp eq haloid[i] AND tform gt t0 AND tform le time[i]/tunit AND radius le radiusofinterest,ningal)
      print,ningal,' in gal from t = ',t0*tunit,time[i]
  ;    stop
      if ningal le 10 then begin
         clumpystarmass[i] = 0.
         coldstarmass[i] = 0.
         shockedstarmass[i] = 0.
         earlystarmass[i] = 0.
      endif else begin

         match,clumpy,gasiord[ingal],cl1,cl2
         match,shocked,gasiord[ingal],s1,s2
         match,cold,gasiord[ingal],c1,c2
         match,early,gasiord[ingal],e1,e2

         clumpystarmass[i] = total(s[ingal[cl2]].mass)*munit
         coldstarmass[i] = total(s[ingal[c2]].mass)*munit
         shockedstarmass[i] = total(s[ingal[s2]].mass)*munit
         earlystarmass[i] = total(s[ingal[e2]].mass)*munit
      endelse
    ; now BH sums
      clumpybhtime = where(bhaccretemass.clumpybhtimes gt t0 AND bhaccretemass.clumpybhtimes le time[i]/tunit,nbhclumpy)
      if nbhclumpy ne 0 then clumpyBHmass[i] = total(bhaccretemass.clumpytotal[clumpybhtime])*munit else clumpyBHmass[i] = 0.

      coldbhtime = where(bhaccretemass.coldbhtimes gt t0 AND bhaccretemass.coldbhtimes le time[i]/tunit,nbhcold)
      if nbhcold ne 0 then coldBHmass[i] = total(bhaccretemass.coldtotal[coldbhtime])*munit else coldBHmass[i] = 0.

      shockedbhtime = where(bhaccretemass.shockedbhtimes gt t0 AND bhaccretemass.shockedbhtimes le time[i]/tunit,nbhshock)
      if nbhshock ne 0 then shockedBHmass[i] = total(bhaccretemass.shockedtotal[shockedbhtime])*munit else shockedBHmass[i] = 0.

      earlybhtime = where(bhaccretemass.earlybhtimes gt t0 AND bhaccretemass.earlybhtimes le time[i]/tunit,nbhearly)
      if nbhearly ne 0 then earlyBHmass[i] = total(bhaccretemass.earlytotal[earlybhtime])*munit else earlyBHmass[i] = 0.

      t0 = time[i]/tunit
   endfor

   openw,lun,'radialstartformfraction.dat',/get_lun,width=500
   printf,lun,'DIFFERENTIAL STELLAR MASSES '
   printf,lun,'time  ','cold mass  ','clumpy mass  ','shocked mass   ','early mass'
   for i=0,nfiles-1 do printf,lun,time[i],coldstarmass[i],clumpystarmass[i],shockedstarmass[i],earlystarmass[i]
   close,lun
   free_lun,lun

   openw,lun,'radialBHtformfraction.dat',/get_lun,width=500
   printf,lun,'DIFFERENTIAL BH accretion MASSES '
   printf,lun,'time  ','cold mass  ','clumpy mass  ','shocked mass   ','early mass'
   for i=0,nfiles-1 do printf,lun,time[i],coldBHmass[i],clumpyBHmass[i],shockedBHmass[i],earlyBHmass[i]
   close,lun
   free_lun,lun

endif else begin
   readcol,istherefile[0],time,coldstarmass,clumpystarmass,shockedstarmass,earlystarmass,format='d',/silent
   readcol,'radialBHtformfraction.dat',time,coldbhmass,clumpybhmass,shockedbhmass,earlybhmass,format='d',/silent
endelse

totalstarmass = coldstarmass+shockedstarmass+clumpystarmass+earlystarmass

; normalize to total gas mass in halo 1?
; good sanity check to see if it adds up 

if keyword_set(plotter) then paperplot,filename='radialstartformgasfraction.ps'

plot,redshifts,clumpystarmass/totalstarmass,xtit='Redshift',ytit='Differential Stellar Mass Fraction',xra=[8,3.8],xstyle=1,yra=[0,1],title=plottitle
  oplot,redshifts,clumpystarmass/totalstarmass,color=125
  oplot,redshifts,coldstarmass/totalstarmass,color=80
  oplot,redshifts,shockedstarmass/totalstarmass,color=240
  oplot,redshifts,earlystarmass/totalstarmass,color=grey

if keyword_set(dolegend) then begin
   legend,['clumpy','unshocked','shocked','early'],colors=[125,80,240,grey],lines=[0,0,0,0],/right,charsize=1,position=[4.5,0.55]
   legend,['gas forming stars','gas feeding MBH'],colors=0,lines=[0,2],charsize=1;,position=[7.7,0.55]
endif

; now. I want to put the fraction of gas accreted by the
; BH here as well. 

totalBHfraction = clumpyBHmass+coldBHmass+earlyBHmass+shockedBHmass

    clumpyfraction = clumpyBhmass/totalBHfraction
    coldfraction = coldBHmass/totalBHfraction
    shockedfraction = shockedBHmass/totalBHfraction
    earlyfraction = earlyBHmass/totalBHfraction


    oplot,redshifts,clumpyfraction,color=125,linestyle=2
    oplot,redshifts,coldfraction,color=80,linestyle=2
    oplot,redshifts,shockedfraction,color=240,linestyle=2
    oplot,redshifts,earlyfraction,color=grey,linestyle=2

if keyword_set(plotter) then paperplot,/close

;stop




end
