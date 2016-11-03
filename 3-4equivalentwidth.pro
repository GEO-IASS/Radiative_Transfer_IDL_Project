function profile,a,tau0,u
  ; return a Schuster-Schwarzschild profile
  ; input: a = damping parameter
  ; tau0 = SS layer thickness at line center
  ; u = wavelength array in Doppler units
  ; output: int = intensity array
  Ts=5700
  Tl=4200
  wav=5000.E-8
  int=u
  usize=SIZE(u)   ; IDL SIZE returns array type and dimensions
  for i=0,usize[1]-1 do begin
    tau=tau0 * voigt(a,abs(u[i]))
    int[i]=planck_1(Ts,wav)*exp(-tau) + planck_1(Tl,wav)*(1.-exp(-tau))
  endfor
  return,int 
end
