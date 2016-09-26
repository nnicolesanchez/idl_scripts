pro radialmetalfraction,plotter=plotter,plottitle=plottitle,dolegend=dolegend

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

istherefile = file_search('radialmetalfraction.dat')
if istherefile eq '' then begin


   clumpy = mrdfits('clumpy.accr.iord.fits',0)
   cold = mrdfits('unshock.iord.fits',0)
   shocked = mrdfits('shocked.iord.fits',0)
   accrz = mrdfits('grp1.accrz.fits',0)
   early = mrdfits('early.iord.fits',0)
   bhfile = file_search('centralbh/allgas.sav')
   restore,bhfile[0]
   readcol,'centralbh/out.distance',w,w,w,accretetime,w,gasiord,format='a,a,a,f,a,d',/silent

   coldmetals = fltarr(nfiles)
   shockedmetals = fltarr(nfiles)
   clumpymetals = fltarr(nfiles)
   earlymetals = fltarr(nfiles)
   totalmetals = fltarr(nfiles)
   coldBHmetals = fltarr(nfiles)
   shockedBHmetals = fltarr(nfiles)
   clumpyBHmetals = fltarr(nfiles)
   earlyBHmetals = fltarr(nfiles)
   totalBHmetals = fltarr(nfiles)
   temperature = fltarr(nfiles)
   BHtemperature = fltarr(nfiles)
   totalmetalerr = fltarr(nfiles)
   totalBHmetalerr = fltarr(nfiles)

; initialize time
   t0 = 0.0
   radiusofinterest = 1.5 ; physical kpc
   for i=0,nfiles-1 do begin
      rtipsy,file[i],h,g,d,s
      d=''
      s=''
      grp =  read_lon_array(file[i]+'.amiga.grp')
      grp = grp[0:h.ngas-1]
      iord =  read_lon_array(file[i]+'.iord')
      iord = iord[0:h.ngas-1]
      readcol,file[i]+'.shrinkcenters',haloidc,xx,yy,zz,format='i,d,d,d',/silent
      location = where(haloidc eq haloid[i])
      xc = xx[location[0]]
      yc = yy[location[0]]
      zc = zz[location[0]]
      ; recenter stars
      gx = g.x-xc
      gy = g.y-yc
      gz = g.z-zc
      radius = sqrt(gx*gx+gy*gy+gz*gz)*h.time*lunit ; physical kpc

    ; find stars in tform range
      ingal = where(grp eq haloid[i] AND radius le radiusofinterest,ningal)
      print,ningal,' in gal from t = ',t0*tunit,time[i]
  ;    stop
      if ningal le 10 then begin
         clumpymetals[i] = 0.
         coldmetals[i] = 0.
         shockedmetals[i] = 0.
         earlymetals[i] = 0.
      endif else begin

         match,clumpy,iord[ingal],cl1,cl2
         match,shocked,iord[ingal],s1,s2
         match,cold,iord[ingal],c1,c2
         match,early,iord[ingal],e1,e2

         if n_elements(cl1) ge 1 then clumpymetals[i] = mean(g[ingal[cl2]].zmetal)
         if n_elements(s1) ge 1 then coldmetals[i] = mean(g[ingal[c2]].zmetal)
         if n_elements(c1) ge 1 then shockedmetals[i] = mean(g[ingal[s2]].zmetal)
         if n_elements(e1) ge 1 then earlymetals[i] = mean(g[ingal[e2]].zmetal)
      endelse


      totalmetals[i] = mean(g[ingal].zmetal)
      temperature[i] = mean(g[ingal].tempg)

      ; what's the districution of metals (and temp)?
      startp = [0,.1,1.]
      plothist,g[ingal].zmetal,exx,whyx,bin=.005,/noplot
      errors = replicate(1.0,n_elements(exx))
      gaussx = mpfitfun("MYGAUSS",exx,whyx,errors,startp,/quiet)
      ; gaussx[1] is sigma
      totalmetalerr[i] = gaussx[1]

    ; now BH sums
      timerange= where(accretetime gt t0 AND accretetime le time[i]/tunit,ntimerange)
      accretediords = gasiord[timerange]
      ; match gas accreted during this step with gas IDs
      match,accretediords,iord,i1,i2
      iordbh = iord[i2]
      print,n_elements(accretediords)-n_elements(i1),' gas particles lost'
      
      match,allgas.clumpy,iordbh,bhcl1,bhcl2
      match,allgas.cold,iordbh,bhc1,bhc2
      match,allgas.shocked,iordbh,bhs1,bhs2
      match,allgas.early,iordbh,bhe1,bhe2

      if n_elements(bhcl1) ne 0 then clumpyBHmetals[i] = mean(g[i2[bhcl2]].zmetal)
      if n_elements(bhc1) ne 0 then coldBHmetals[i] = mean(g[i2[bhc2]].zmetal)
      if n_elements(bhs1) ne 0 then shockedBHmetals[i] = mean(g[i2[bhs2]].zmetal)
      if n_elements(bhe1) ne 0 then earlyBHmetals[i] = mean(g[i2[bhe2]].zmetal)

      totalbhmetals[i] = mean(g[i2].zmetal)
      Bhtemperature[i] = mean(g[i2].tempg)

      plothist,g[i2].zmetal,exx,whyx,bin=.005,/noplot
      errors = replicate(1.0,n_elements(exx))
      gaussx = mpfitfun("MYGAUSS",exx,whyx,errors,startp,/quiet)
                                ; gaussx[1] is sigma                                                                                                          
      totalBHmetalerr[i] = gaussx[1]

      ; end of BH metal calculation
      t0 = time[i]/tunit
   endfor

   ; print files
   openw,lun,'radialmetalfraction.dat',/get_lun,width=500
   printf,lun,'DIFFERENTIAL STELLAR METALSES '
   printf,lun,'time  ','cold metals  ','clumpy metals  ','shocked metals   ','early metals', ' temperature', '  metals error'
   for i=0,nfiles-1 do printf,lun,time[i],coldmetals[i],clumpymetals[i],shockedmetals[i],earlymetals[i],totalmetals[i],temperature[i],totalmetalerr[i]
   close,lun
   free_lun,lun

   openw,lun,'radialBHmetalfraction.dat',/get_lun,width=500
   printf,lun,'DIFFERENTIAL BH accretion METALSES '
   printf,lun,'time  ','cold metals  ','clumpy metals  ','shocked metals   ','early metals', '  temperature','   metals error'
   for i=0,nfiles-1 do printf,lun,time[i],coldBHmetals[i],clumpyBHmetals[i],shockedBHmetals[i],earlyBHmetals[i],totalbhmetals[i],BHtemperature[i],totalBHmetalerr[i]
   close,lun
   free_lun,lun

endif else begin
   readcol,istherefile[0],time,coldmetals,clumpymetals,shockedmetals,earlymetals,totalmetals,temperature,totalmetalerr,format='d',/silent
   readcol,'radialBHmetalfraction.dat',time,coldbhmetals,clumpybhmetals,shockedbhmetals,earlybhmetals,totalbhmetals,BHtemperature, totalbhmetalerr,format='d',/silent
endelse


; normalize to total gas metals in halo 1?
; good sanity check to see if it adds up 

if keyword_set(plotter) then paperplot,filename='radialmetalfraction.ps'

plot,redshifts,totalmetals,xtit='Redshift',ytit='Metallicity (Z)',xra=[8,3.8],xstyle=1,title=plottitle,yra=[0,.05]
  oplot,redshifts,totalmetals,thick=8
  oplot,redshifts,totalmetalerr/2.+totalmetals,linestyle=1
  oplot,redshifts,totalmetals-totalmetalerr/2.,linestyle=1

;  oplot,redshifts,clumpymetals,color=125
;  oplot,redshifts,coldmetals,color=80
;  oplot,redshifts,shockedmetals,color=240
;  oplot,redshifts,earlymetals,color=grey
 
if keyword_set(dolegend) then begin
;   legend,['clumpy','unshocked','shocked','early','total'],colors=[125,80,240,grey,0],lines=[0,0,0,0,0],/right,charsize=1,position=[4,0.01]
   legend,['gas near BH','gas feeding MBH'],colors=0,lines=[0,2],charsize=1;,position=[7.7,0.55]
endif

; now. I want to put the fraction of gas accreted by the
; BH here as well. 


;    oplot,redshifts,clumpybhmetals,color=125,linestyle=2
;    oplot,redshifts,coldbhmetals,color=80,linestyle=2
;    oplot,redshifts,shockedbhmetals,color=240,linestyle=2
;    oplot,redshifts,earlybhmetals,color=grey,linestyle=2
    oplot,redshifts,totalbhmetals,linestyle=2,thick=8
    oplot,redshifts,totalBHmetalerr/2.+totalBHmetals,linestyle=3
    oplot,redshifts,totalBHmetals-totalBHmetalerr/2.,linestyle=3



if keyword_set(plotter) then paperplot,/close

;stop


if keyword_set(plotter) then paperplot,filename='radialtempfraction.ps'

plot,redshifts,temperature,xtit='Redshift',ytit='Temp (K)',xra=[8,3.8],xstyle=1,title=plottitle,/ylog,yra=[1e5,1e7]

if keyword_set(dolegend) then   legend,['gas near BH','gas feeding MBH'],colors=0,lines=[0,2],charsize=1;,position=[7.7,0.55]                          

    oplot,redshifts,bhtemperature,linestyle=2
if keyword_set(plotter) then paperplot,/close

end



FUNCTION MYGAUSS, X, P
  RETURN,   GAUSS1(X, P[0:2])
END
