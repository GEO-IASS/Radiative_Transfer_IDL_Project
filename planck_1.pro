function planck_1, temp, wav    ;planck function in terms of temp and wavelength
h= 6.62607D-27 ;planck constant in erg*s
c=2.99792D10  ;speed of light in cm*s^-1
k= 1.38065D-16 ; Boltzmann constant in ergK^-1
const= 2*h*c^2
econst= (h*c)/k 
evar= temp*wav
pfunc= (const/(wav^5))*1/(exp(econst/(evar))-1)
return, pfunc

end