;******************************************************************************
;Plot the rotation curve of an exponential disk.
;******************************************************************************

pro wmap3_vcdisk

Rd = 4.0
disk_mass = 1e10

;Define Newton's constant in units of kpc * km^2 / (Msun * s^2)
G = 4.467e-6

r = findgen(1000)

y = r / (2.0 * Rd)

central_sb = disk_mass / (2.0 * !PI * Rd^2.0)

vc = sqrt((4.0 * !PI * G * central_sb * Rd * y^2.0) * (beseli(y, 0) * beselk(y,0) - beseli(y, 1) * beselk(y,1)))

norm = sqrt(G * disk_mass / Rd)

myps, filename = 'vcdisk.ps'

plot, r / Rd, vc / norm, xtit = textoidl('r / R_d'), ytit = textoidl('V_c / (GM/R_d)^{1/2}'), xrange = [0, 10], yrange = [0, 0.8], /xstyle, /ystyle

psclose

end
