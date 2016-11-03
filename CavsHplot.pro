set_plot, 'ps'
device, filename='CavsHplot.ps'
temp=indgen(191)*100.+1000  ; T = 1000-20000 in delta T = 100
CaH = temp    ; declare ratio array
Caabund=2.E-6 ; A_Ca = N_Ca / N_H
for i=0,190 do begin
  NCa = sahabolt_Ca(temp[i],1e2,2,1)
  NH = sahabolt_H(temp[i],1e2,2)
  CaH[i]=NCa*Caabund/NH
endfor
plot,temp,CaH,/ylog,$
  xtitle=’temperature’,ytitle=’Ca II K / H alpha’
  device,/close
  filenmame='CavsHplot.ps'
  end