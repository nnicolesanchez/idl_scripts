function wmap3_dist_modulus, z

;COMPUTE DISTANCE MODULUS AT A GIVEN REDSHIFT

return, 5.0 * alog10(wmap3_lum_dist(z)) - 5.0

end
