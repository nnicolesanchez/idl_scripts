; #### Created: April 2015 #### Last edit: April 14, 2015 ####
; This program creates a histogram of which halos the gas particles
; accreted by the black hole immediately originated from prior to the 
; merger.

pro halohistogram_clumpy,plotter=plotter,plottitle=plottitle
loadct,39
tunit=1.223128407d18/3.1556926d7/1d9 ; Gyr

;#### We want to read in iords of all particles 
;#### Get into each timestep and find iord list
haloidfile = file_search('*haloid.dat')
readcol,haloidfile[0],file,haloid,format='a,i',/silent
file	   = reverse(file)
nfiles     = n_elements(file)
print,'This should be your number of timesteps: ',nfiles

;#### We need to read in the iords for all the gas that has been
;#### accreted by the central black hole
readcol,'centralbh/out.distance',w,bhiords,w,eattime,w,bhgasiords,dsquare,w,w,smooth,w,eatenmass,format='a,l,a,f,a,l,f,a,a,f,a,f',/silent

;#### So we also only want to look at the gas being accreted after
;#### z=0.8 so that we are looking only at gas from the merger
;#### Need eat times in reasonable units
eattime_Gyr       = eattime * tunit
;print,eattime_Gyr[1000:1500]
index_gtz0_8 = where(eattime_Gyr gt 6.4)
print,n_elements(bhgasiords_gtz0_8)

test = eattime_Gyr[index_gtz0_8]
print,test[0:100]
eattime_gtz0_8    = eattime[index_gtz0_8]
print,eattime_gtz0_8[0:100]
bhgasiords_gtz0_8 = bhgasiords[index_gtz0_8]
print,bhgasiords_gtz0_8[0:100]

;#### Using i = 12 (00192) for testing
for i = 12,nfiles-1 do begin
	rtipsy,file[i],h,g,d,s
	
;#### We only want gas particle iords and haloids
;#### 'iords' : particles IDs for ALL particles: gas,dark matter,stars (in that order)
;#### 'amiga' : particle halo IDs for ALL particles
	iords 	   = read_lon_array(file[i]+'.iord')
	amiga	   = read_lon_array(file[i]+'.amiga.grp')

;#### So to cut off the excess star and dark matter IDs
	gasiords   = iords[0:h.ngas-1]
	gashalos   = amiga[0:h.ngas-1]

;#### But we only want to look at the gas being accreted through mergers
;#### (clumpy gas) Need indices of gasiords that match clumpy
	match,bhgasiords_gtz0_8,gasiords,cl1,cl2
        ;print,n_elements(cl2)
        ghist      = histogram(gashalos[cl2])
        ;print,ghist[0:15]
        print,i
        stop
	max_ghist  = max(ghist)
		annoying   = i + 1
	timestep   = strtrim(annoying,1)	

;#### For plotting
	if keyword_set(plotter) then paperplot,filename='clumpyhistogram_gtz0_8_'+timestep+'.ps'

;#### DO I WANT TO REMOVE ALL 0's and 1's?????
;#### Since gasiords index matches gashalos
	plothist,gashalos[cl2],/ylog,yra=[1,10000],xra=[1,50],$
	title='Histogram of BH Accretion in Halos at Timestep (gt z=0.8): '+timestep,  $
	xtitle='Halo IDs',thick=1,color=90;,/fill,fcolor=90;,fcolor=170,/fill;,/noplot,/fill

        ;x_tot = n_elements(xhist)-1
        ;y_tot = n_elements(yhist)-1
        ;plot,xhist[2:x_tot],yhist[2:y_tot],yra=[1,155000],xra=[1,100],/ylog,histogram=1
	if keyword_set(plotter) then paperplot,/close

;trace halo id of gas particles accreted by BH back in time not just clumpy
endfor
end


; remove everything before merger event at z=0
; define when merger event happens
; time when halo id's switch is about where merger happens
; define that tiem (maybe a snapshot before that)
; that point onward halo 2 is fueling
; determine which histogram is the important one
