; recompute as arrays and overplot relative
set_plot, 'ps'
device, filename='2-9dipplot.ps'
NCa=temp      ; declare array
NH=temp       ; declare array
for i=0,100 do begin
  NCa[i] = sahabolt_ca(temp[i],1e2,2,1) ; Ca ion ground state
  NH[i] = sahabolt_H(temp[i],1e2,2)     ; H atom 2nd level
endfor
oplot,temp,NH/max(NH)
oplot,temp,NCa/max(NCa),linestyle=2     ; Ca curve again dashed
device,/close
filenmame='2-9dipplot.ps'
end