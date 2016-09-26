;##### Stellar Fractions vs. Time #####
;This program plots the stellar fraction of a galaxy across time and redshift.

pro galstarfraction,plotter=plotter,galaxyname=galaxyname

;plotter    = plot saved to .ps file rather than to screen
;galaxyname = name of galaxy apears on plot title
;ASK HOW TO DO BELOW
;plottitle  = 'Star Fraction of '+galaxyname
if not keyword_set(galaxyname) then galaxyname='Star Fraction of Galaxy'

;##### Units #####
lunit = 50000.     ; [kpc]
munit = 1.84793d16 ; M_sun
grey  = cgcolor('Dark Gray')
loadct,4

;###### Read In Files #####
haloidfile = file_search('*haloid.dat')
readcol,haloidfile[0],file,haloid,format='a,i',/silent
file   = reverse(file)
nfiles = n_elements(file) 

clumpy = mrdfits('clumpy.accr.iord.fits',0)
cold = mrdfits('unshock.iord.fits',0)
shocked = mrdfits('shocked.iord.fits',0)
grp1 = mrdfits('grp1.allgas.iord.fits',0)
early = mrdfits('early.iord.fits',0)

;##### Determine Redshifts for X-Axis #####
haloid    = reverse(haloid) ;backwards list of haloids
accrz     = mrdfits('grp1.accrz.fits',0)
redshifts = accrz[uniq(accrz,sort(accrz))]
;We only want to analyze between z=0 and where tracing begins
;!!!!!!! Do I want to hardcode in 4.5???
redshifts = redshifts[where(redshifts ge 0. AND redshifts le 4.5)] 
redshifts = [4.5,reverse(redshifts)]

print,'Do I have the same number of timesteps',nfiles,' as I do red shifts',n_elements(redshifts),'?'
stop


;##### Star Fraction Table #####
;Checks to see if table has already been created.
;starfractiontable = file_search('starfraction.dat')
;if starfractiontable[0] eq '' then begin

;Determine which particles are early to find the gas in each timestep
notearly = grp1
match,early,notearly,e1,e2
remove,e2,notearly

starfractiontable = file_search('starfraction.dat')
if starfractiontable[0] eq '' then begin

;For each timestep, we must read in the gas iord for the parent of the star particles. 
;Then match them to the type of gas. 
;Count the number of stars of that type and calculate the fraction compared to the total gas in that timestep.

coldfrac    = fltarr(nfiles)
clumpyfrac  = fltarr(nfiles)
shockedfrac = fltarr(nfiles)
earlyfrac   = fltarr(nfiles)

for i=0,nfiles-1 do begin
	starparentiord = read_lon_array(file[i]+'.igasorder')
	print,'Number of stars in galaxy at timestep',i,' is',n_elements(starparentiord)
	
	match,cold,starparentiord,c1,c2
	match,clumpy,starparentiord,cl1,cl2
	match,shocked,starparentiord,s1,s2
	match,early,starparentiord,e1,e2
	coldstars=n_elements(c1)
	clumpystars=n_elements(cl1)
	shockedstars=n_elements(s1)
	earlystars=n_elements(e1)
	totalstars=coldstars+clumpystars+shockedstars+earlystars

	coldfrac[i]    = float(coldstars)/totalstars
	clumpyfrac[i]  = float(clumpystars)/totalstars
	shockedfrac[i] = float(shockedstars)/totalstars
	earlyfrac[i]   = float(earlystars)/totalstars

endfor

	openw,lun,'starfraction.dat',/get_lun,width=500
	printf,lun,'cold # frac   ','clumpy # frac   ','shocked # frac   ','early * frac   '
	for i=0,nfiles-1 do printf,lun,coldfrac[i],clumpyfrac[i],shockedfrac[i],earlyfrac[i]
	close,lun
	free_lun,lun

	stop

;##### Read in Fraction Table if Already Created #####
endif else readcol,'starfraction.dat',coldfrac,clumpyfrac,shockedfrac,earlyfrac,format='d',/silent


;##### Plot Star Fraction vs. Redshift
if keyword_set(plotter) then paperplot,filename='starfraction.ps'

plot,redshifts,coldfrac,xtit='Redshift',ytit='Star Fraction',xra=[4.5,0],xstyle=1,yra=[0,1],linestyle=0,thick=3,title=galaxyname
	oplot,redshifts,clumpyfrac,color=90,linestyle=0,thick=3
	oplot,redshifts,coldfrac,color=60,linestyle=0,thick=3
	oplot,redshifts,shockedfrac,color=170,linestyle=0,thick=3
	oplot,redshifts,earlyfrac,color=0,linestyle=0,thick=3

legend,['cold','clumpy','shocked','early'],colors=[60,90,170,0],lines=[0,0,0,0],/right,charsize=1,position=[0.1,0.98]

if keyword_set(plotter) then paperplot,/close

end
