set_plot, 'ps'
device, filename='3-5.ps'
tau0=10^(indgen(61)/10.-2.) ; 10^-2 to 10^4, 0.1 steps in the log
eqw=tau0
; same size array
for i=0,60 do begin
  int=profile(a,tau0[i],u)
  reldepth=(int[0]-int)/int[0]
  eqw[i]=total(reldepth)*0.4
endfor
plot,tau0,eqw,xtitle='tau0',ytitle='equivalent width',/xlog,/ylog

device,/close
filenmame='3-5.ps'
end
