function wmap3_absolute_magnitude, m, z, K

;COMPUTES THE ABSOLUTE MAGNITUDE, GIVEN AN APPARENT MAGNITUDE AND A
;REDSHIFT AND (OPTIONAL) K CORRECTION


if (n_params() eq 2) then K=0

return, m - 5.0 * alog10(wmap3_lum_dist(z)) + 5.0 - K

end
