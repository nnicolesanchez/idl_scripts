function wmap3_ang_diam_dist, z

;COMPUTES THE ANGULAR DIAMETER DISTANCE IN PARSECS TO AN OBJECT AT REDSHIFT Z

return, wmap3_prop_motion_dist(z) / (1.0 + z)

end
