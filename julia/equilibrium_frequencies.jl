#println("Time-Stepping for the reduced Center-Manifold Approximation")

#include("read_amplitudes.jl")

include("amplitudes.jl")

#---------------------------------------------------------------------- 

using LinearAlgebra
using FileIO

freqsave    = true
plotedge    = true 

# Mode numbering starts from 0.
# A1  - Modes 5/6
# A2  - Modes 3/4

#σ      = -9.1854401740e-02
σ           = -0.092
cmap        = get_cmap("tab10")

VBase       = zeros(ComplexF64,nmodes)
VBase[1]    = ϵ3
VBase[2]    = σ0

ind1        = 5
ind1c       = Conj[ind]

ind2        = 3
ind2c       = Conj[ind2]
#---------------------------------------------------------------------- 


close("all")

cmap = get_cmap("tab10")

h10,ax10 = subplots(1,1,figsize=[7,5],gridspec_kw=Dict("width_ratios"=>[1]),layout="constrained",clear=true)
ax10.set_ylabel(L"Ω",fontsize=lafs)
ax10.set_xlabel(L"Re",fontsize=lafs)

h11,ax11 = subplots(1,1,figsize=[7,5],gridspec_kw=Dict("width_ratios"=>[1]),layout="constrained",clear=true)
ax11.set_ylabel(L"|z|",fontsize=lafs)
ax11.set_xlabel(L"Re",fontsize=lafs)

nsteps1   = 1000  
nsteps2   = 1000
nsteps3   = 1000
nsteps    = nsteps1+nsteps2+nsteps3
# λ1,2 - LCO
Rei_end   = 1.0/5000.0
ϵ_end     = 1.0 - Rei_end/Reci
ϵs1       = 1.0e-6
ϵs2       = ϵ3
ϵs3       = ϵ4
ϵs4       = ϵ_end
IND1      = 1:nsteps1
IND2      = (nsteps1+1):(nsteps1+nsteps2)
IND3      = (nsteps1+nsteps2+1):nsteps

#ϵrange    = LinRange(ϵ_start,ϵ_end,nsteps)
ϵrange1   = LinRange(ϵs1,ϵs2,nsteps1)
ϵrange2   = LinRange(ϵs2,ϵs3,nsteps2)
ϵrange3   = LinRange(ϵs3,ϵs4,nsteps3)
ϵrange    = [ϵrange1; ϵrange2; ϵrange3]

ReRange1  = Rec./(1.0 .- ϵrange)
Ω12       = zeros(Float64,nsteps)
A12       = zeros(Float64,nsteps)

ind1      = 5
ind1c     = Conj[ind1]
for i in 1:nsteps
  VB      = zeros(ComplexF64,nmodes)
  VB[1]   = ϵrange[i]
  VB[2]   = σ0
 
  A1eq,Ω1eq  = GetEquilibriumAmpOmega(G,VB,ind1,ind1c)
  Ω12[i] = Ω1eq
  A12[i] = A1eq
end
ax10.plot(ReRange1,Ω12,linestyle="-",linewidth=2,color=cmap(0),label=L"λ_{1,2}-LCO")

ax11.plot(ReRange1[IND1],A12[IND1],linestyle="-",linewidth=2,color=cmap(0),label=L"z_{1}")
ax11.plot(ReRange1[IND2],A12[IND2],linestyle="-",linewidth=2,color=cmap(0),label=L"z_{1}")
ax11.plot(ReRange1[IND3],A12[IND3],linestyle=":",linewidth=2,color=cmap(0),label=L"z_{1}")


# λ3,4 LCO
#-------------------------------------------------- 
ϵ_start   = ϵ2
Rei_end   = 1.0/5000.0
ϵ_end     = 1.0 - Rei_end/Reci
ϵs1       = ϵ2
ϵs2       = ϵ3
ϵs3       = ϵ4
ϵs4       = ϵ_end

ϵrange1   = LinRange(ϵs1,ϵs2,nsteps1)
ϵrange2   = LinRange(ϵs2,ϵs3,nsteps2)
ϵrange3   = LinRange(ϵs3,ϵs4,nsteps3)

ϵrange    = [ϵrange1; ϵrange2; ϵrange3]

ReRange2  = Rec./(1.0 .- ϵrange)
Ω34       = zeros(Float64,nsteps)
A34       = zeros(Float64,nsteps)

ind2        = 3
ind2c       = Conj[ind2]
for i in 1:nsteps
  VB       = zeros(ComplexF64,nmodes)
  VB[1]    = ϵrange[i]
  VB[2]    = σ0
 
  A2eq,Ω2eq = GetEquilibriumAmpOmega(G,VB,ind2,ind2c)
  Ω34[i]    = Ω2eq
  A34[i]    = A2eq
end
ax10.plot(ReRange2,Ω34,linestyle="-",linewidth=2,color=cmap(1),label=L"λ_{3,4}-LCO")

ax11.plot(ReRange2[IND1],A34[IND1],linestyle=":" ,linewidth=2,color=cmap(1),label=L"z_{3}")
ax11.plot(ReRange2[IND2],A34[IND2],linestyle="-" ,linewidth=2,color=cmap(1),label=L"z_{3}")
ax11.plot(ReRange2[IND3],A34[IND3],linestyle="-" ,linewidth=2,color=cmap(1),label=L"z_{3}")



# Edge-State
#---------------------------------------------------------------------- 
ϵ_start   = ϵ3
ϵ_end     = ϵ4

ϵrange    = LinRange(ϵ_start,ϵ_end,nsteps)
ReRange3  = Rec./(1.0 .- ϵrange)
Ω12_edge  = zeros(Float64,nsteps)
Ω34_edge  = zeros(Float64,nsteps)
A12_edge  = zeros(Float64,nsteps)
A34_edge  = zeros(Float64,nsteps)

for i in 1:nsteps
  VB       = zeros(ComplexF64,nmodes)
  VB[1]    = ϵrange[i]
  VB[2]    = σ0
 
  edge     = GetEdgeState(amplitude_z3,amplitude_z5,ϵrange[i],σ0)
  if isempty(edge)
    println("i: $i, ϵ:$(ϵrange[i])")
  else
    A34_edge[i] = edge[1]
    A12_edge[i] = edge[2]
  end  
  VB        = zeros(ComplexF64,nmodes)
  VB[1]     = ϵrange[i]
  VB[2]     = σ0
  VB[ind1]  = A12_edge[i]
  VB[ind1c] = A12_edge[i]

  VB[ind2]  = A34_edge[i]
  VB[ind2c] = A34_edge[i]
  #A1eq,Ω1eq = GetEquilibriumAmpOmega(G,VB,ind1,ind1c)
  #Ω12_edge[i] = Ω1eq
  λe1       = GetEffectiveEigenvalue(G,VB,ind1,ind1c)
  Ω12_edge[i] = imag(λe1)

  #A2eq,Ω2eq = GetEquilibriumAmpOmega(G,VB,ind2,ind2c)
  #Ω34_edge[i] = Ω2eq
  λe2       = GetEffectiveEigenvalue(G,VB,ind2,ind2c)
  Ω34_edge[i] = imag(λe2)

end

if plotedge
  ax10.plot(ReRange3,Ω12_edge,color=cmap(0),linestyle="--")
  ax10.plot(ReRange3,Ω34_edge,color=cmap(1),linestyle="--")
end  


ax10.plot([Re2; Re2],[6.0; 13.0],linestyle="-.",color="gray",linewidth=1)
ax10.plot([Re3; Re3],[6.0; 13.0],linestyle="-.",color="gray",linewidth=1)
ax10.plot([Re4; Re4],[6.0; 13.0],linestyle="-.",color="gray",linewidth=1)
ax10.set_ylim([7.0; 11.0])
ax10.set_xlim([Rec; 5000.0])

ax11.plot(ReRange3,A12_edge,linestyle="--",linewidth=2,color=cmap(0),label=L"z_{1}")
ax11.plot(ReRange3,A34_edge,linestyle="--",linewidth=2,color=cmap(1),label=L"z_{3}")
ax11.plot([Re2; Re2],[0.0; 0.03],linestyle="-.",color="gray",linewidth=1)
ax11.plot([Re3; Re3],[0.0; 0.03],linestyle="-.",color="gray",linewidth=1)
ax11.plot([Re4; Re4],[0.0; 0.03],linestyle="-.",color="gray",linewidth=1)
ax11.set_ylim([0.0; 0.02])
ax11.set_xlim([Rec; 5000.0])

D = load("cavity_nonlinear_omega.jld2")
Ω_nonlinear = get(D,"Ω",[])
ReAll       = get(D,"Reall",[])
ax10.plot(ReAll,Ω_nonlinear,color="black",linestyle="none",marker="o",markersize=6)

if (freqsave)
  svname = "equilibrium_omega.eps"
  h10.savefig(pdir*svname)

  svname = "bifurcation_diagram.eps"
  h11.savefig(pdir*svname)
end  

@printf "Done."
















