function wmap3_lbfunc, a

;THIS FUNCTION DEFINES THE INTEGRAND IN THE INTEGRAL THAT HAS TO BE COMPUTED
;IN ORDER TO COMPUTE THE LOOKBACK TIME.  SEE 'THE COSMOLOGICAL CONSTANT' BY
;SEAN CARROLL, EQUATION 36.

wmap3_set_cosmology,H_0,omega_r,omega_m,omega_v,omega_k

H = H_0 * 1.02528 * 10.0^(-12.0) ;in years^-1

H_a = H * (omega_r / a^4 + omega_m / a^3 + omega_v + omega_k / a^2)^0.5

return, 1.0 / (a * H_a)

end
