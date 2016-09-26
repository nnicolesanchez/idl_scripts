;******************************************************************************
;Calculate rho_crit as a function of z, in units of solar masses
;versus cubic kiloparsecs.
;******************************************************************************

function wmap3_rho_crit, z

return, 136.05 * (wmap3_Hubble(z) / 73.0)^2.0

end
