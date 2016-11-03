; temperature sensitivity CaIIK and Halpha
set_plot, 'ps'
device, filename='2-9CavsHtemp.ps'
temp=indgen(101)*100.+2000.   ; T = 2000-12000, delta T = 100
dNCadT = temp                 ; declare array
dNHdT = temp                   ; declare array
dT=1.
for i=0,100 do begin
  NCa = sahabolt_ca(temp[i],1e2,2,1)      ; Ca ion ground state
  NCa2 = sahabolt_ca(temp[i]-dT,1e2,2,1)  ; idem dT cooler
  dNCadT[i]= (NCa - NCa2)/dT/NCa          ; fractional diff quotient
  NH = sahabolt_H(temp[i],1e2,2)          ; H atom 2nd level
  NH2 = sahabolt_H(temp[i]-dT,1e2,2)      ; idem dT cooler
  dNHdT[i] = (NH-NH2)/dT/NH               ; fractional diff quotient
endfor
plot,temp,abs(dNHdT),/ylog,yrange=[1E-5,1],$
  xtitle=’temperature’,ytitle=’abs d n(r,s) / n(r,s)’
oplot,temp,abs(dNCadT),linestyle=2        ; Ca curve dashed
device,/close
filenmame='2-9CavsHtemp.ps'
end