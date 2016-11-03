set_plot, 'ps'
device, filename='plancktempwavLOGY.ps'
wav=indgen(100)*200.+1000.  ; produces wav[0,...99] = 1000 - 20800
print,wav                   ; check that
b=wav                       ; declare float array of the same size
c=wav
for i=0,99 do b[i]=planck_1(8000.,wav[i]*1.E-8)
cgplot,wav,b,/WINDOW, ASPECT =0.8,/ylog,xtitle='wavelength (Angstrom)',ytitle='Planck function', $
  xmargin = [45,5],$
  charsize=1.1              ; bigger characters
for T=8000.,5000.,-200. do begin ; step from 8000 K down to 5000 K
  for i=0,99. do b[i]=planck_1(T,wav[i]*1.E-8)
 cgplot,wav,b,/window, /OVERPLOT              ; overplots extra curves in existing graph
endfor                      ; begin...end sequences can’t go on command line
device,/close
filenmame='plancktempwavLOGY.ps'
end