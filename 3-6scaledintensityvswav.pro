set_plot, 'ps'
device, filename='3-3scaledintensityvswav.ps'
u=indgen(201)/10.-10.   ; u = -10 to 10 in 0.1 steps
int=u                   ; declare array
for iwav=1,3 do begin
  wav=(iwav^2+1)*1.D-5      ; wav = 2000, 5000, 10000 Angstrom
  for itau=0,8 do begin
    for i=0,200 do begin
      tau=tau0[itau] * voigt(a,abs(u[i]))
      int[i]=planck_1(Ts,wav) * exp(-tau) + planck_1(Tl,wav)*(1.-exp(-tau))
    endfor
    int=int/int(0)          ; convert into relative intensity
    if (iwav eq 1 and itau eq 0) then cgplot,u,int, xtitle='wavelength', ytitle='Ilambda/Icontinuum-relative I', color='pink'
    if (iwav eq 1 and itau gt 0) then cgplot,u,int, /overplot, color='blue'
    if (iwav eq 2) then cgplot,u,int,/overplot,linestyle=1, color='orange' ; dotted
    if (iwav eq 3) then cgplot,u,int,/overplot,linestyle=4, color='red' ; dash dot dot dot
  endfor
endfor
device,/close
filenmame='3-3scaledintensityvswav.ps'
end
