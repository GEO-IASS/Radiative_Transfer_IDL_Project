set_plot, 'ps'
device, filename='3-2emergentintensitylogxy.ps'
B=2.
tau=indgen(101)/10.+0.01
; set array tau = 0.01-10 in steps 0.01
int=tau
; declare float array of the same size
for I0=4,0,-1 do begin
  ; step down from I0=4 to I0=0
  for i=0,100 do int[i]=I0 * exp(-tau[i]) + B*(1-exp(-tau[i]))
  if (i0 eq 4) then plot,tau,int,/xlog, /ylog,$
    xtitle='tau',ytitle='Intensity',charsize=1.3
  if (i0 ne 4) then oplot,tau,int
endfor
device,/close
filenmame='3-2emergentintensitylogxy.ps'
end
