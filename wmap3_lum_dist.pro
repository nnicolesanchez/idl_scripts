function wmap3_lum_dist, z

;COMPUTES THE LUMINOSITY DISTANCE IN PARSECS TO AN OBJECT AT REDSHIFT Z

return, wmap3_prop_motion_dist(z) * (1.0 + z)

end
