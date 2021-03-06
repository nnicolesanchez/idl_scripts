pro histcombine

; read in fits files for central BH
; combine them so I have the histories of all the gas particles in a more simple format.

; gas iords I care about - central BH has eaten these.
; this file is the combined gasiords file for all central BH IDs
readcol,'gasiords.good.all',gasiord,meaten,format='l,f',/silent
ngas = n_elements(gasiord)
 
; for particles in halo 1

h1clumpy = mrdfits('../clumpy.accr.iord.fits',0)
h1shocked = mrdfits('../shocked.iord.fits',0)
h1cold = mrdfits('../unshock.iord.fits',0)
h1early = mrdfits('../early.iord.fits',0)
h1smooth = mrdfits('../smooth.accr.iord.fits',0)
h1tracestep = mrdfits('../tracedtostep.fits',0)
h1accrz = mrdfits('../grp1.accrz.fits',0)
h1allgas = mrdfits('../grp1.allgas.iord.fits',0)


; pare down these arrays so only the iords from the BH accretion are involved.
match,h1early,h1allgas,e1,e2
notearly=h1allgas
lremove,e2,notearly   ; all gas that isn't early.
; notearly should match the redshift array.
match,notearly,gasiord,ne1,ne2  ; WEIRD PROBLEM STEP
notearly = notearly[ne1]
goodredshifts = h1accrz[ne1]
; now pare down the smooth and tracestep arrays
match,h1smooth,gasiord,sm1,sm2
smooth=h1smooth[sm1]
tracestep = h1tracestep[sm1]
; pare down cold and shocked and clumpy
match,h1shocked,gasiord,sh1,sh2
match,h1cold,gasiord,cold1,cold2
match,h1clumpy,gasiord,cl1,cl2
shocked=h1shocked[sh1]
cold=h1cold[cold1]
clumpy=h1clumpy[cl1]
; and also early and allgas
match,h1early,gasiord,e11,e22
early = h1early[e11]
;;; ALL IS WELL THUS FAR  ;;;;
match,h1allgas,gasiord,all1,all2
yayallgas = h1allgas[all1]
nyay = n_elements(all1)
; set accrz array to have high z values for early particles.
; first get early indices right
match,yayallgas,early,e111,e222
goodinds = lindgen(nyay)
remove,e111,goodinds  ; inds without early particles
zaccr = fltarr(n_elements(yayallgas))
zaccr[e111] = 50. ; for early particles
; now for the rest of the indices
poo = yayallgas
remove,e111,poo
match,yayallgas,poo,y1,y2

zaccr[goodinds] = goodredshifts; h1accrz does not match h1allgas

accrz = zaccr

allidentified = [shocked,cold,clumpy,early]

if n_elements(allidentified) ne n_elements(gasiords) then  print,n_elements(gasiord)-n_elements(allidentified)  ,' lost particles'


stop

iords = gasiord

; sweet. save into one nice file.
allgas = {clumpy:clumpy,zaccr:accrz,early:early,cold:cold,shocked:shocked,tracestep:tracestep,smooth:smooth,iords:yayallgas}

save,allgas,filename='allgas.sav'

;stop
end
