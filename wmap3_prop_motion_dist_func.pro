function wmap3_prop_motion_dist_func, a

;THIS FUNCTION DEFINES THE INTEGRAND IN THE INTEGRAL THAT HAS TO BE COMPUTED
;IN ORDER TO COMPUTE MOST INTERESTING COSMOLOGICAL PARAMETERS.  SEE
;'THE COSMOLOGICAL CONSTANT' BY SEAN CARROLL, EQUATION 42 FOR THE PROPER
;MOTION DISTANCE.

wmap3_set_cosmology,H_0,omega_r,omega_m,omega_v,omega_k

H = H_0 * 1000.0 * 100.0 / ((10.0^6.0) * (3.09 * 10.0^18.0)) ;in sec^-1

H_a = H * (omega_r / a^4 + omega_m / a^3 + omega_v + omega_k / a^2)^0.5

return, 1.0 / (a^2.0 * H_a)

end
