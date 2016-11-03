set_plot, 'ps'
device, filename='3-4plot1.ps'
u=indgen(1001)/2.5-200.
a=0.1
tau0=1e2
int=profile(a,tau0,u)
cgplot,u,int

reldepth=(int[0]-int)/int[0]  ;line depth in relative units
cgplot,u,reldepth
eqw=total(reldepth)*0.4       ;integral = TOTAL times interval
print,eqw

device,/close
filenmame='3-4plot1.ps'
end
