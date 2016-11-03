set_plot, 'ps'
device, filename='resultsplot.ps'
temp=1000*indgen(31) ; make array 0,...,30000 in steps of 1000 K
print,temp ; check
pop=fltarr(5,31) ; declare float array for n(r,T)

for T=1,30 do $ ; $ continues statement to next line
for r=1,4 do pop[r,T]=sahabolt_E(temp[T],131.,r,1)

cgPlot, temp, pop[1,*], /ylog,yrange=[1E-3,1.1], xrange=[0,32000], /window,$
BACKGROUND=white, TITLE='R = 1, S = 0,2,4 Populations n_r,1/N',$
COLOR='Grn3', XTITLE='temp',YTITLE='population', LEGENDS='r1=grn3'
cgPlot, temp, pop[2,*], /WINDOW, /OVERPLOT,$
COLOR='Red3', LEGENDS= 'r2=red3'
cgPlot, temp, pop[3,*], /WINDOW, /OVERPLOT,$
COLOR='Blu3', LEGENDS= 'r3=blu3'
cgPlot, temp, pop[4,*], /WINDOW, /OVERPLOT,$
COLOR='ORG3', LEGENDS= 'r4=org3' 

for T=1,30 do $; repeat for s=2 (excitation energy = 1 eV)
  for r=1,4 do pop[r,T]=sahabolt_E(temp[T],131.,r,2)
  cgPlot, temp, pop[1,*], /window,/OVERPLOT, $
    COLOR='GrN5', XTITLE='temp',YTITLE='population', LEGENDS='r1=grn5''
  cgPlot, temp, pop[2,*], /WINDOW, /OVERPLOT,$
    COLOR='Red5', LEGENDS= 'r2=red5'
  cgPlot, temp, pop[3,*], /WINDOW, /OVERPLOT,$
    COLOR='Blu5', LEGENDS= 'r3=blu5'
  cgPlot, temp, pop[4,*], /WINDOW, /OVERPLOT,$
    COLOR='org5', LEGENDS= 'r4=org5'

for T=1,30 do $
  ; repeat for s=4 (excitation energy = 3 eV)
  for r=1,4 do pop[r,T]=sahabolt_E(temp[T],131.,r,4)
  cgPlot, temp, pop[1,*], /window,/OVERPLOT, $
    COLOR='GrN7', XTITLE='temp',YTITLE='population', LEGENDS='r1=grn7'
  cgPlot, temp, pop[2,*], /WINDOW, /OVERPLOT,$
    COLOR='Red7', LEGENDS= 'r2=red7'
  cgPlot, temp, pop[3,*], /WINDOW, /OVERPLOT,$
    COLOR='Blu7', LEGENDS= 'r3=blu7'
  cgPlot, temp, pop[4,*], /WINDOW, /OVERPLOT,$
    COLOR='org7', LEGENDS= 'r4=org7'
device,/close
filenmame='resultsplot.ps'
end 