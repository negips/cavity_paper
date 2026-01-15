#println("Time-Stepping for the reduced Center-Manifold Approximation")

include("ReadAmplitudeFile.jl")
include("$(JULIAHOME)/center_manifold/Module_CenterManifold/CenterManifold.jl")

using PyPlot
using Peaks
using Statistics
using Printf
using Random

include("$JULIACOMMON/SetZero.jl")

lafs   = 18
lgfs   = 16

Rec    = 4.131331106555E+03
#Rec    = 1.09503844263237E+02
#Rec    = 46.30
Reci   = 1.0/Rec
# ν      = 0.02
# Rei    = Reci*(1.0 - ν)
# Re     = 1.0/Rei

Re     = 4500.0
Rei    = 1.0/Re
ν      = 1.0 - Rei/Reci 

pdir   = "plots/"
cm     = get_cmap("tab10")

if Re<4400.0
  LineSt = "-"
else
  LineSt = "--"
end  

ifsave = true

# close("all")
#FOL    = "cyl_n9_x40_O6_p1/"
#FOL    = "cavity_custom4/"
#FOL    = "cyl_3modes_O6/"
#FOL    = "cyl_4modes/"
#FOL    = "cyl_medium_3modes_O8/"
#FOL    = "cyl_subd_5modes/"

#FOL    = "cavity_5modes/"
#FOL    = "cavity_5modes_O4_PERT/"
FOL    = "cavity_normal_form_pert/"

nmodes = 6

fnames = Vector{String}(undef,nmodes)

for i in 1:nmodes
  fnames[i] = @sprintf("%s%s%.2d",FOL,"AmplitudeTerms.",i-1)
end  

Coeffs = Vector{Matrix{ComplexF64}}(undef,nmodes)
for i in 1:nmodes
  AMat      = ReadAmplitudeFile(fnames[i])
  #Coeffs[i] = AMat[1:35,1:3]
  #r,c = size(AMat)
  #for j in 1:c
  #  for k in 11:r
  #    AMat[k,j] = 0.0
  #  end
  #end
  Coeffs[i] = copy(AMat)
end  

#Coeffs[2][2,1] = Coeffs[2][2,1] - 9.1854401740e-02
#Coeffs[3][3,1] = Coeffs[3][3,1] - 9.1854401740e-02


dt     = 0.001
nsteps = 400000


rng    = Xoshiro(124);

amp    = 1.00e-3
v0     = amp*rand(rng,ComplexF64,nmodes)
v0[1]  = Complex(ν)
v0[2]  = Complex(-9.1854401740e-02)

#v0[2]  = Complex(0.05)
#v0[3]  = Complex(0.0001)
#v0[4]  = v0[3]'

# Do Complex conjugations of modes
v0[4] = v0[3]'
if (nmodes == 6)
  #v0[4] = 0.0
  v0[6] = v0[5]'
end  


hist   = zeros(Float64,nsteps,nmodes)
histt  = zeros(Float64,nsteps)

time   = 0.0
v      = copy(v0)

for i in 1:nsteps
  global time, v

  EvolveSystemN!(Coeffs,v,dt,1)

  v[4]      = v[3]'
  v[6]      = v[5]'

  hist[i,:] = real.(v)
  time      = time + dt
  histt[i]  = time
end

MPeaks       = 20
ind2         = 5
peak_ind2    = argmaxima(hist[:,ind2])
peak_amp2    = hist[peak_ind2,ind2]
peak_times2  = histt[peak_ind2]
delta_times2 = diff(peak_times2)
ω2           = 2.0*π./delta_times2
npeaks       = length(ω2)
if npeaks > MPeaks
 ωend        = ω2[npeaks-MPeaks:npeaks]
 A2end       = peak_amp2[npeaks-MPeaks:npeaks]
else
 ωend        = ω2
 A2end       = peak_amp2
end
Ω2        = mean(ωend)
Ω2var     = std(ωend)
A2        = mean(A2end)
@printf("ν: %g ; Re = %g ; Ω: %.3e; |A|: %.3e \n", ν, Re, Ω2, A2)


ind4         = 3
peak_ind4    = argmaxima(hist[:,ind4])
peak_amp4    = hist[peak_ind4,ind4]
peak_times4  = histt[peak_ind4]
delta_times4 = diff(peak_times4)
ω4           = 2.0*π./delta_times4
npeaks       = length(ω4)
if npeaks > MPeaks
 ωend        = ω4[npeaks-MPeaks:npeaks]
 A4end       = peak_amp4[npeaks-MPeaks:npeaks]
else
 ωend        = ω4
 A4end       = peak_amp4
end
Ω4        = mean(ωend)
Ω4var     = std(ωend)
A4        = mean(A4end)
@printf("ν: %g ; Re = %g ; Ω: %.3e ; |A|: %.3e \n", ν, Re, Ω4, A4)

# 
# h1 = figure(num=1)
# plot(histt,histr,linestyle="-",label=L"Center-Manifold: \mathcal{O}(3)")

# h1.savefig("center_manifold_omega.eps")

h2  = figure(num=2,figsize=[8,6],layout="constrained")
h2.clf()
ax2 = h2.subplots()
ax2.plot(histt,hist[:,ind2],linestyle=LineSt,label=L"z_{1}",color=cm(0))
ax2.plot(histt,hist[:,ind4],linestyle=LineSt,label=L"z_{3}",color=cm(1))
ax2.set_ylabel(L"\mathfrak{R}(z)",fontsize=lafs)
ax2.set_xlabel(L"t",fontsize=lafs)
ax2.legend(fontsize=lgfs,loc="upper center",ncols=2)

if ifsave
  tlast = histt[end]
  tmax  = tlast +  0.0
  tmin  = tlast - 10.0

  ax2.set_xlim([tmin; tmax])

  ymax  = maximum([A2;A4])*1.3
  ymin  = -ymax
  ax2.set_ylim([ymin; ymax])

  svname = @sprintf("response_Re%.0f.eps",Re)
  h2.savefig(pdir*svname)
end


#h3  = figure(num=3,figsize=[8,6],layout="constrained")
#ax3 = h3.subplots()

lab0 = @sprintf(" - Re=%.0f",Re)
lab1 = "z_{1}"*lab0
lab3 = "z_{3}"*lab0

ax3.plot(peak_times2,peak_amp2,linestyle=LineSt,label=L"%$(lab1)",color=cm(0))
ax3.plot(peak_times4,peak_amp4,linestyle=LineSt,label=L"%$(lab3)",color=cm(1))
ax3.set_ylabel(L"\mathfrak{R}(z)^{pk}",fontsize=lafs)
ax3.set_xlabel(L"t",fontsize=lafs)
ax3.legend(fontsize=lgfs,loc="upper left",ncols=2)

if ifsave
  tlast = histt[end]
  tmax  = tlast +  0.0
  tmin  = 0.0
 
  ax3.set_ylim([-1.0e-3; 1.2e-2])
  ax3.set_xlim([tmin; tmax])
 
  svname = @sprintf("peak_response.eps")
  h3.savefig(pdir*svname)
end





