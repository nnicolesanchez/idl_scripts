;******************************************************************************
;What is the Hubble parameter as a function of z?
;******************************************************************************

function wmap3_hubble, z

wmap3_set_cosmology,H_0,omega_r,omega_m,omega_v,omega_k

a = 1.0 / (1.0 + z)

return, H_0 * (omega_r / a^4 + omega_m / a^3 + omega_v + omega_k / a^2)^0.5

end
