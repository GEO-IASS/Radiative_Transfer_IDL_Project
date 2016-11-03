set_plot, 'ps'
device, filename='2-10plot.ps'
temp=indgen(191)*100.+1000. ; array 1000 - 20 000 in steps 1000
nH=temp                     ; declare same size array
for i=0,190 do nH[i]=sahabolt_H(temp[i],1e2,1)
plot,temp,nH,$
  xtitle=’temperature’,ytitle=’neutral hydrogen fraction’
  device,/close
  filenmame='2-10plot.ps'
end