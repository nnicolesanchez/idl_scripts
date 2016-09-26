function wmap3_angular_size2, proper_size, ang_diam_distance

;COMPUTES THE ANGULAR SIZE OF AN OBJECT IN ARSECONDS, GIVEN ITS PROPER
;SIZE IN PARSECS, AND ITS ANGULAR DIAMETER DISTANCE IN PARSECS.

return, proper_size * 60.0 * 60.0 * 180.0 / (ang_diam_distance * !PI)

end
