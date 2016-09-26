pro wmap3_set_cosmology,H_0,omega_r,omega_m,omega_v,omega_k

;THIS PROCEDURE SETS THE CONSTANTS FOR A GIVEN COSMOLOGY
;ALL OTHER FUNCTIONS IN THIS DIRECTORY CALL THIS PROCEDURE
;AS NEEDED, SO CONSTANTS NEED ONLY BE CHANGED ONCE

H_0 = 73.0 ;km/s
omega_r = 0.0
omega_m = 0.24
omega_v = 0.76
;omega_m = 0.2
;omega_v = 0.0
omega_k = 1-omega_r-omega_m-omega_v

end
