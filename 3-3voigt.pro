set_plot, 'ps'
device, filename='3-3voigtalog.ps'
u=indgen(201)/10.-10.                     ;u = -10 to 10 in 0.1 steps
vau=u                                     ;declare same-size array
a=5                                     ;damping parameter
for i=0,200 do vau[i]=voigt(a,abs(u[i]))  ;taking abs corrects IDL errors
plot,u,vau,/ylog                  ;yrange fixed to compare plots
print, vau
device,/close
filenmame='3-3voigtalog.ps'
end
