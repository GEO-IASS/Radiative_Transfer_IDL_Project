set_plot, 'ps'
device, filename='3-3schusterprofile.ps'
Ts=5700                 ; solar surface temperature
Tl=4200                 ; solar T-min temperature = ‘reversing layer’
a=0.1                   ; damping parameter
wav=5000.D-8            ; wavelength in cm
tau0=.95                  ; reversing layer thickness at line center
u=indgen(201)/10.-10.   ; u = -10 to 10 in 0.1 steps
int=u                   ; declare array
for i=0,200 do begin
  tau=tau0 * voigt(a,abs(u[i]))
  int[i]=planck_1(Ts,wav) * exp(-tau) + planck_1(Tl,wav)*(1.-exp(-tau))
endfor
cgplot,u,int, aspect=.82, xtitle='u-dimensionless wavelength', ytitle='emergent rad. Intensity with Voigt function', color = 'green'

tau0=[0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100]
for itau=0,8 do begin
  for i=0,200 do begin
    tau=tau0[itau] * voigt(a,abs(u[i]))
    int[i]=planck_1(Ts,wav) * exp(-tau) + planck_1(Tl,wav)*(1.-exp(-tau))
  endfor
  cgplot,u,int, color='blue', /overplot
endfor
device,/close
filenmame='3-3schusterprofile.ps'
end
