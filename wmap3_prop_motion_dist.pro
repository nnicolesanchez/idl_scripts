function wmap3_prop_motion_dist, z

;COMPUTES THE PROPER MOTION DISTANCE IN PARSECS TO AN OBJECT AT REDSHIFT Z
;SEE 'THE COSMOLOGICAL CONSTANT' BY SEAN CARROLL, EQUATION 42

wmap3_set_cosmology,H_0,omega_r,omega_m,omega_v,omega_k

H = H_0 * 1000.0 * 100.0 / ((10.0^6.0) * (3.09 * 10.0^18.0)) ;in sec^-1

integ = qsimp('wmap3_prop_motion_dist_func', 1.0/(1.0+z), 1.0)

if (omega_k eq 0) then x=integ else begin
    R_0 = 1.0 / (H * sqrt(abs(omega_k)))
    if (omega_k gt 0) then x=R_0*sinh(integ/R_0) else x=R_0*sin(integ/R_0)
endelse

return, (3.0 * 10.0^10.0 * x) / (3.09 * 10.0^18.0)

end
