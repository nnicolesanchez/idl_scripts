function wmap3_proper_size, arcsecs, z

;COMPUTES THE PROPER SIZE OF AN OBJECT IN PARSECS, GIVEN ITS ANGULAR EXTENT
;IN ARCSECONDS, AND ITS REDSHIFT

return, wmap3_ang_diam_dist(z) * arcsecs * 3.14159628 / (60.0 * 60.0 * 180.0)

end
