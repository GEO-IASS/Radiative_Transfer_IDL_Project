function partfunc_E, temp
; partition functions Schadee element
; input: temp (K)
; output: fltarr(4) = partition functions U1,..,U4
u=fltarr(4)
chiion=[7,16,31,51]
k=8.61734D-5
for r = 0,3 do $
for s = 0, chiion[r]-1 do u[r]=u[r] + exp(-s/(k*temp))
return,u
end

function boltz_E,temp,r,s
  ; compute Boltzmann population for level r,s of Schadee element E
  ; input: temp (temperature, K)
  ;r (ionization stage nr, 1 - 4 where 1 = neutral E)
  ;s (level nr, starting at s=1)
  ; output: relative level population n_(r,s)/N_r
  u=partfunc_E(temp)
  keV=8.61734D-5
  ; Boltzmann constant in ev/deg
  relnrs = 1./u[r-1]*exp(-(s-1)/(keV*temp))
  return, relnrs
end

function saha_E,temp,elpress,ionstage
  ; compute Saha population fraction N_r/N for Schadee element E
  ; input: temperature, electron pressure, ion stage
  ; physics constants
  kerg=1.380658D-16     ; Boltzmann constant (erg K; double precision)
  kev=8.61734D-5        ; Boltzmann constant (eV/deg)
  h=6.62607D-27         ; Planck constant (erg s)
  elmass=9.109390D-28 ; electron mass (g)
  ; kT and electron density
  kevT=kev*temp
  kergT=kerg*temp
  eldens=elpress/kergT
  chiion=[7,16,31,51]     ; ionization energies for element E
  u=partfunc_E(temp)      ; get partition functions U[0]...u[3]
  u=[u,2]                 ; add estimated fifth value to get N_4 too
  sahaconst=(2*!pi*elmass*kergT/(h*h))^1.5 * 2./eldens
  nstage=dblarr(5)        ; double-precision float array
  nstage[0]=1.            ; relative fractions only (no abundance)
  for r=0,3 do $
    nstage[r+1] = nstage[r]*sahaconst*u[r+1]/u[r]*exp(-chiion(r)/kevT)
  ntotal=total(nstage)    ; sum all stages = element density
  nstagerel=nstage/ntotal ; fractions of element density
  return,nstagerel[ionstage-1] ; ion stages start at 1, IDL at 0
end

function sahabolt_E,temp,elpress,ion,level
  ; compute Saha-Boltzmann populaton n_(r,s)/N for level r,s of E
  ; input: temperature, electron pressure, ionization stage, level nr
  return, saha_E(temp,elpress,ion) * boltz_E(temp,ion,level)
end

function sahabolt_H,temp,elpress,level
  ; compute Saha-Boltzmann population n_(1,s)/N_H for hydrogen level
  ; input: temperature, electron pressure, level number
  ; physics constants
  kerg=1.380658D-16   ; Boltzmann constant (erg K; double precision)
  kev=8.61734D-5      ; Boltzmann constant (eV/deg)
  h=6.62607D-27       ; Planck constant (erg s)
  elmass=9.109390D-28 ; electron mass (g)
  ; kT and electron density
  kevT=kev*temp
  kergT=kerg*temp
  eldens=elpress/kergT        ; energy levels and weights for hydrogen
  nrlevels=100                ; reasonable partition function cut-off value
  g=intarr(2,nrlevels)        ; declaration weights (too many for proton)
  chiexc=fltarr(2,nrlevels)   ; declaration excitation energies (idem)
  for s=0,nrlevels-1 do begin ; enclose multiple lines with begin...end
    g[0,s]=2*(s+1)^2                  ; statistical weights
    chiexc[0,s]=13.598*(1-1./(s+1)^2) ; excitation energies
  endfor                      ; begin...end cannot go on command line!
  g[1,0]=1                    ; statistical weight free proton
  chiexc[1,0]=0.              ; excitation energy proton ground state
  ; partition functions
  u=fltarr(2)
  u[0]=0
  for s=0,nrlevels-1 do u[0]=u[0]+ g[0,s]*exp(-chiexc[0,s]/kevT)
  u[1]=g[1,0]
  ; Saha
  sahaconst=(2*!pi*elmass*kergT/(h*h))^1.5 * 2./eldens
  nstage=dblarr(2)        ; double-precision float array
  nstage[0]=1.            ; relative fractions only
  nstage[1] = nstage[0] * sahaconst * u[1]/u[0] * exp(-13.598/kevT)
  ntotal=total(nstage)    ; sum both stages = total hydrogen density
  ; Boltzmann
  nlevel = nstage[0]*g[0,level-1]/u[0]*exp(-chiexc[0,level-1]/kevT)
  nlevelrel=nlevel/ntotal ; fraction of total hydrogen density
  ;stop                   ; in for parameter inspection
  return,nlevelrel
end

function partfunc_Ca, temp
  ; partition functions Ca+
  ; input: temp (K)
  ; output: fltarr(4) = partition functions U1,..,U4
  u=fltarr(4)
  chiion=[6.113,11.871,50.91,67.15]
  k=8.61734D-5
  for r = 0,3 do $
    for s = 0, chiion[r]-1 do u[r]=u[r] + exp(-s/(k*temp))
  return,u
end

function boltz_Ca,temp,r,s
  ; compute Boltzmann population for level r,s of Schadee element E
  ; input: temp (temperature, K)
  ;r (ionization stage nr, 1 - 4 where 1 = neutral E)
  ;s (level nr, starting at s=1)
  ; output: relative level population n_(r,s)/N_r
  u=partfunc_E(temp)
  keV=8.61734D-5
  ; Boltzmann constant in ev/deg
  relnrs = 1./u[r-1]*exp(-(s-1)/(keV*temp))
  return, relnrs
end

function saha_Ca,temp,elpress,ionstage
  ; compute Saha population fraction N_r/N for Schadee element E
  ; input: temperature, electron pressure, ion stage
  ; physics constants
  kerg=1.380658D-16     ; Boltzmann constant (erg K; double precision)
  kev=8.61734D-5        ; Boltzmann constant (eV/deg)
  h=6.62607D-27         ; Planck constant (erg s)
  elmass=9.109390D-28 ; electron mass (g)
  ; kT and electron density
  kevT=kev*temp
  kergT=kerg*temp
  eldens=elpress/kergT
  chiion=[6.113,11.871,50.91,67.15]     ; ionization energies for element E
  u=partfunc_Ca(temp)      ; get partition functions U[0]...u[3]
  u=[u,2]                 ; add estimated fifth value to get N_4 too
  sahaconst=(2*!pi*elmass*kergT/(h*h))^1.5 * 2./eldens
  nstage=dblarr(5)        ; double-precision float array
  nstage[0]=1.            ; relative fractions only (no abundance)
  for r=0,3 do $
    nstage[r+1] = nstage[r]*sahaconst*u[r+1]/u[r]*exp(-chiion(r)/kevT)
  ntotal=total(nstage)    ; sum all stages = element density
  nstagerel=nstage/ntotal ; fractions of element density
  return,nstagerel[ionstage-1] ; ion stages start at 1, IDL at 0
end

function sahabolt_Ca,temp,elpress,ion,level
  ; compute Saha-Boltzmann populaton n_(r,s)/N for level r,s of E
  ; input: temperature, electron pressure, ionization stage, level nr
  return, saha_Ca(temp,elpress,ion) * boltz_Ca(temp,ion,level)
end

