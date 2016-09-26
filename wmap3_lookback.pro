function wmap3_lookback, z

;COMPUTES THE LOOKBACK TIME IN YEARS TO AN OBJECT AT REDSHIFT Z

return, qsimp('wmap3_lbfunc', 1.0/(1.0 + z), 1.0)

end
