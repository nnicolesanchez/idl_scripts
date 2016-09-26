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

;##### Files to be Read In #####
haloidfile = file_search('*haloid.dat')
readcol,haloidfile[0],file,haloid,format='a,i',/silent
file   = reverse(file)
nfiles = n_elements(file) 

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
starfractiontable = file_search('starfraction.dat')
if starfractiontable[0] eq '' then begin

;Creates table which will hold the total number of stars in each timestep.
stars = fltarr(nfiles)

for i=0,nfiles-1 do begin
	rtipsy,file[i],h,g,d,s
	stars[i] = n_elements(s)
	print,stars[i]
	;stop
endfor

openw,lun,'starfraction.dat',/get_lun,width=500
printf,lun,'Timestep #    ','Number of Stars'
for i=0,nfiles-1 do printf,lun,i+1,stars[i]
close,lun
free_lun,lun

endif else readcol,'starfraction.dat',timestep,stars,format='d',/silent

;##### Total Number of Stars in Galaxy #####
for i=0,nfiles-1 do begin
	start      = 0
	starstotal = start + stars[i]
endfor
print,'The total number of stars in this galaxy is:',starstotal
;stop

;##### Star Fraction ##### 
starfraction = stars/starstotal

;##### Plot Star Fraction vs. Redshift
if keyword_set(plotter) then paperplot,filename='starfraction.ps'

plot,redshifts,starfraction,xtit='Redshift',ytit='Star Fraction',xra=[4.5,0],xstyle=1,yra=[0,1],linestyle=0,thick=3,title=galaxyname

if keyword_set(plotter) then paperplot,/close

end
