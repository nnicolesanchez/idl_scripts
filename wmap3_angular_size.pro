function wmap3_angular_size, proper_size, z

;COMPUTES THE ANGULAR SIZE OF AN OBJECT IN ARSECONDS, GIVEN ITS PROPER
;SIZE IN PARSECS, AND ITS REDSHIFT

return, proper_size * 60.0 * 60.0 * 180.0 / (wmap3_ang_diam_dist(z) * !PI)

end
