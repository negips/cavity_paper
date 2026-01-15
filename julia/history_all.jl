#!/bin/julia

println("Nonlinear Ω for open cavity.")

using CSV
using DataFrames
using PyPlot
using FFTW
using Statistics
using Peaks
using Printf
using JLD2

# close("all")

lafs   = 16
ifsave = true

fols   = Vector{String}(undef,11)
pref   = "../../NekTools/cavity/nonlinear/"


fols[1]   = "re4140/"
fols[2]   = "re4200/"
fols[3]   = "re4300/"
fols[4]   = "re4400/"
fols[5]   = "re4500/"
fols[6]   = "re4550/"
fols[7]   = "re4600/"
fols[8]   = "re4700/"
fols[9]   = "re4800/"
fols[10]  = "re4900/"
fols[11]  = "re5000/"

ncases = length(fols)

Reall  = Vector{Float64}(undef,ncases)
Reall[1]   = 4140.0
Reall[2]   = 4200.0
Reall[3]   = 4300.0
Reall[4]   = 4400.0
Reall[5]   = 4500.0
Reall[6]   = 4550.0
Reall[7]   = 4600.0
Reall[8]   = 4700.0
Reall[9]   = 4800.0
Reall[10]  = 4900.0
Reall[11]  = 5000.0


Ω     = zeros(Float64,ncases)

for i in 1:ncases
  global Ω

  fol   = fols[i]
  
  local fname = "cavity.his"
  
  file  = pref*fol*fname
  
  f     = CSV.File(file,delim=' ',header=false,skipto=4,stripwhitespace=true,ignorerepeated=true)
  df    = DataFrame(f)
  
  l     = length(df[:,1])
  
  ind1  = 1:2:l
  ind2  = 2:2:l
  
  # Point 1
  time1 = df[ind1,1]
  vx1   = df[ind1,2]
  vy1   = df[ind1,3]
  
  # Point 2
  time2 = df[ind2,1]
  vx2   = df[ind2,2]
  vy2   = df[ind2,3]
  
  # h1  = figure(num=1)
  # ax1 = h1.subplots()
  # ax1.plot(time1,vx1)
  
  L     = length(vx1)
  i1    = L - 2^13 + 1
  i2    = L
  S     = vx1[i1:i2]
  S0    = mean(S)
  S     = S .- S0
  T     = time1[i1:i2]
  DT    = T[2] - T[1]
  fs    = 1.0/DT
  
  # h2  = figure(num=2)
  # ax2 = h2.subplots()
  # ax2.plot(T,S)
  
  peak_ind    = argmaxima(S)
  peak_vals   = S[peak_ind]
  peak_times  = T[peak_ind]
  delta_times = diff(peak_times)
  ω           = 2.0*π./delta_times
  npeaks      = length(ω)
  
  mpeaks      = 9
  if npeaks > mpeaks
   ωend        = ω[npeaks-mpeaks+1:npeaks]
  else
   ωend        = ω
  end
  Ω[i]     = mean(ωend)

  @printf("%s:, %.3e:, Ω: %.4f \n", fol, Reall[i], Ω[i])
end  

h3 = figure()
ax3 = h3.subplots()
ax3.plot(Reall,Ω,linestyle="none",marker="o",label="Full-System")
ax3.set_ylabel(L"Ω",fontsize=lafs)
ax3.set_xlabel(L"Re",fontsize=lafs)
legend(loc="center left")

if ifsave
   svname = "cavity_nonlinear_omega.jld2"
   save(svname,"Reall",Reall,"Ω",Ω);
   println(svname*" saved.")
end  





