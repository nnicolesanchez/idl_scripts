pro coldgasfraction_quies

; what fraction of this galaxy is 
 ; - cold gas
 ; - clumpy gas
 ; - shocked gas

smooth = mrdfits('../smooth.accr.iord.fits',0)
shocked = mrdfits('../shocked.iord.fits',0)
clumpy = mrdfits('../clumpy.accr.iord.fits',0)
allgas = mrdfits('../grp1.allgas.iord.fits',0)
cold = mrdfits('../unshock.iord.fits',0)
early = mrdfits('../early.iord.fits',0)

ngas = n_elements(allgas)
ncold = n_elements(cold)
nearly = n_elements(early)
nshocked = n_elements(shocked)
nclumpy = n_elements(clumpy)

print,'for galaxy h277'

print,'cold/total fraction = ',float(ncold)/float(ngas)

print,'clumpy/total fraction = ',float(nclumpy)/float(ngas)

print,'shocked/total fraction = ',float(nshocked)/float(ngas)

print,'early/total fraction = ',float(nearly)/float(ngas)




end
