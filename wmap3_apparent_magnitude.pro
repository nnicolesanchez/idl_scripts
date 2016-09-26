function wmap3_apparent_magnitude, M, z

;COMPUTE THE APPARENT MAGNITUDE, GIVEN AN ABSOLUTE MAGNITUDE AND A
;REDSHIFT

return, M + 5.0 * alog10(wmap3_lum_dist(z)) - 5.0

end
