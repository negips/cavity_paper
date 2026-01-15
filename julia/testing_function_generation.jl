function GetAmplitudeSq(ind::Int,var::String,G::Vector{Matrix{T}},Conj::Vector{Int};tol=1.0e-10) where {T<:ComplexF64}

  nv   = length(Conj)
  OrdM = length(G)

  H    = Vector{Matrix{T}}(undef,OrdM+1)
  H1   = zeros(T,nv,nv)
  H[1] = H1
  # Check for purely real ones.


  # Assuming we do have conjugates
  i1 = ind
  i2 = Conj[i1]
  for ord in 1:OrdM
    
    ord1 = ord+1
    nt1  = CenterManifold.NInteractionTerms(ord1,nv)
    Hi   = zeros(T,nv,nt1)
    h1   = zeros(T,nt1)
    h2   = zeros(T,nt1)

    nt   = CenterManifold.NInteractionTerms(ord,nv)


    
    g1   = G[ord][i1,:]
    ind1 = zeros(Int,ord+1)
    ind2 = zeros(Int,ord+1)
   
    for i in 1:nt
      ind0        = CenterManifold.GetInteractionTerms(i,ord,nv)

      copyto!(ind1,1,ind0,1,ord)
      ind1[ord1]  = i2-1
      n1          = CenterManifold.GetIndexNumber(ind1,nv)
      val1        = G[ord][i1,i]
      Hi[i1,n1]   = val1

      copyto!(ind2,1,ind0,1,ord)
      ind2[ord1]  = i1-1
      n2          = CenterManifold.GetIndexNumber(ind2,nv)
      val2        = G[ord][i2,i]
      Hi[i2,n2]   = val2
    end
    h1            = Hi[i1,:] .+ Hi[i2,:]

    Hi[i1,:]      = h1
    Hi[i2,:]      = h1

    H[ord1] = copy(Hi)
  end    


  return H
end
#----------------------------------------------------------------------

function GetAmpNullClines(H,par1,par2,ind1,ind1c,ind2,ind2c)

  MOrd = length(H)
  NV   = size(H[1],2)

  # |z3|^2
  #-------------------------------------------------- 
  # 1st order
  # |z3|^2
  i1_1  = [ind1; ind1c] .- 1
  n1_1  = CenterManifold.GetIndexNumber(i1_1,NV)
  v1_1  = real.(H[2][ind1,n1_1])
  
  
  # 2nd order
  # ϵ*|z3|^2
  i2_1  = [par1; ind1; ind1c] .- 1
  n2_1  = CenterManifold.GetIndexNumber(i2_1,NV)
  v2_1  = real.(H34[3][ind1,n2_1])
  
  # σ'*|z3|^2
  i2_2  = [par2; ind1; ind1c] .- 1
  n2_2  = CenterManifold.GetIndexNumber(i2_2,NV)
  v2_2  = real.(H34[3][ind1,n2_2])
  
  # 3rd order
  # ϵ^2*|z3|^2
  i3_1  = [par1; par1; ind1; ind1c] .- 1
  n3_1  = CenterManifold.GetIndexNumber(i3_1,NV)
  v3_1  = real.(H34[4][ind1,n3_1])
  
  # ϵ*σ'*|z3|^2
  i3_2  = [par1; par2; ind1; ind1c] .- 1
  n3_2  = CenterManifold.GetIndexNumber(i3_2,NV)
  v3_2  = real.(H34[4][ind1,n3_2])
  
  # σ'^2*|z3|^2
  i3_3  = [par2; par2; ind1; ind1c] .- 1
  n3_3  = CenterManifold.GetIndexNumber(i3_3,NV)
  v3_3  = real.(H34[4][ind1,n3_3])
  
  # LHS
  #------------------------------  
  # |z3|^2*|z3|^2
  i3_4  = [ind1; ind1c; ind1; ind1c] .- 1
  n3_4  = CenterManifold.GetIndexNumber(i3_4,NV)
  v3_4  = real.(H34[4][ind1,n3_4])
  L_11  = -v3_4         # L3_34
  
  # |z5|^2*|z3|^2
  i3_5  = [ind2; ind2c; ind1; ind1c] .- 1
  n3_5  = CenterManifold.GetIndexNumber(i3_5,NV)
  v3_5  = real.(H34[4][ind1,n3_5])
  L_12  = -v3_5         # L3_56
  
  amplitude_z1(z1,z2,ϵ,σ) = v1_1 + v2_1*ϵ + v2_2*σ + v3_1*ϵ^2 + v3_2*ϵ*σ + v3_3*σ^2 - (z1^2)*L_11 - (z2^2)*L_12
  # nullc_z1(ϵ,σ) = v1_1 + v2_1*ϵ + v2_2*σ + v3_1*ϵ^2 + v3_2*ϵ*σ + v3_3*σ^2
  # nullc_z1e(ϵ)  = nullc_z3(ϵ,σ0)

  return L_11,L_12,amplitude_z1 
end
#---------------------------------------------------------------------- 

using Roots

include("read_amplitudes.jl")

cmap  = get_cmap("tab10")

#-------------------------------------------------- 
Conj = [1; 2; 4; 3; 6; 5]     # Conjugate Indices
SInd = [5; 6; 3; 4; 1; 2]     # Swapped Indices

#σ0   = -9.1854401740e-02
σ0   = -9.2e-02

ind  = 5
ind2 = Conj[ind]
#CenterManifold.DisplayTerms(ind,"z",G)
CenterManifold.DisplayTermsSwapped(ind,"z",G,SInd)

#CenterManifold.DisplayTerms(ind2,"z",G)

ind  = 3
ind2 = Conj[ind]

H34 = GetAmplitudeSq(ind,"z",G,Conj)
CenterManifold.DisplayTerms(ind,"z",H34)
#CenterManifold.DisplayTermsSwapped(ind,"z",H34,SInd)

ind  = 5
H56 = GetAmplitudeSq(ind,"z",G,Conj)
#CenterManifold.DisplayTerms(ind,"z",H56)
#CenterManifold.DisplayTermsSwapped(ind,"z",H56,SInd)

#---------------------------------------------------------------------- 

# Positions of values
ind1  = 3
ind2  = Conj[ind1]
ind3  = 5
ind4  = Conj[ind3]
HH    = H34

NV    = 6         # No of variables

# # |z3|^2
# #-------------------------------------------------- 
# # 1st order
# # |z3|^2
# i1_1  = [ind1; ind2] .- 1
# n1_1  = CenterManifold.GetIndexNumber(i1_1,NV)
# v1_1  = real.(H34[2][ind1,n1_1])
# 
# 
# # 2nd order
# # ϵ*|z3|^2
# i2_1  = [1; ind1; ind2] .- 1
# n2_1  = CenterManifold.GetIndexNumber(i2_1,NV)
# v2_1  = real.(H34[3][ind1,n2_1])
# 
# # σ'*|z3|^2
# i2_2  = [2; ind1; ind2] .- 1
# n2_2  = CenterManifold.GetIndexNumber(i2_2,NV)
# v2_2  = real.(H34[3][ind1,n2_2])
# 
# # 3rd order
# # ϵ^2*|z3|^2
# i3_1  = [1; 1; ind1; ind2] .- 1
# n3_1  = CenterManifold.GetIndexNumber(i3_1,NV)
# v3_1  = real.(H34[4][ind1,n3_1])
# 
# # ϵ*σ'*|z3|^2
# i3_2  = [1; 2; ind1; ind2] .- 1
# n3_2  = CenterManifold.GetIndexNumber(i3_2,NV)
# v3_2  = real.(H34[4][ind1,n3_2])
# 
# # σ'^2*|z3|^2
# i3_3  = [2; 2; ind1; ind2] .- 1
# n3_3  = CenterManifold.GetIndexNumber(i3_3,NV)
# v3_3  = real.(H34[4][ind1,n3_3])
# 
# # LHS
# #------------------------------  
# # |z3|^2*|z3|^2
# i3_4  = [ind1; ind2; ind1; ind2] .- 1
# n3_4  = CenterManifold.GetIndexNumber(i3_4,NV)
# v3_4  = real.(H34[4][ind1,n3_4])
# L3_34 = -v3_4
# 
# # |z5|^2*|z3|^2
# i3_5  = [ind3; ind4; ind1; ind2] .- 1
# n3_5  = CenterManifold.GetIndexNumber(i3_5,NV)
# v3_5  = real.(H34[4][ind1,n3_5])
# L3_56 = -v3_5

L3_34,L3_56,amplitude_z3 = GetAmpNullClines(H34,1,2,ind1,ind2,ind3,ind4)
#amplitude_z3(z3,z5,ϵ,σ) = v1_1 + v2_1*ϵ + v2_2*σ + v3_1*ϵ^2 + v3_2*ϵ*σ + v3_3*σ^2 - (z3^2)*L3_34 - (z5^2)*L3_56
nullc_z3(ϵ,σ)  = amplitude_z3(0.0,0.0,ϵ,σ)
nullc_z3e(ϵ)   = amplitude_z3(0.0,0.0,ϵ,σ0)

#nullc_z3(ϵ,σ) = v1_1 + v2_1*ϵ + v2_2*σ + v3_1*ϵ^2 + v3_2*ϵ*σ + v3_3*σ^2
#nullc_z3e(ϵ)  = nullc_z3(ϵ,σ0)

# # |z5|^2
# #-------------------------------------------------- 
# # 1st order
# # |z5|^2
# i1_1  = [ind3; ind4] .- 1
# n1_1  = CenterManifold.GetIndexNumber(i1_1,NV)
# w1_1  = real.(H56[2][ind3,n1_1])
# 
# # 2nd order
# # ϵ*|z5|^2
# i2_1  = [1; ind3; ind4] .- 1
# n2_1  = CenterManifold.GetIndexNumber(i2_1,NV)
# w2_1  = real.(H56[3][ind3,n2_1])
# 
# # σ'*|z5|^2
# i2_2  = [2; ind3; ind4] .- 1
# n2_2  = CenterManifold.GetIndexNumber(i2_2,NV)
# w2_2  = real.(H56[3][ind3,n2_2])
# 
# # 3rd order
# # ϵ^2*|z5|^2
# i3_1  = [1; 1; ind3; ind4] .- 1
# n3_1  = CenterManifold.GetIndexNumber(i3_1,NV)
# w3_1  = real.(H56[4][ind3,n3_1])
# 
# # ϵ*σ'*|z5|^2
# i3_2  = [1; 2; ind3; ind4] .- 1
# n3_2  = CenterManifold.GetIndexNumber(i3_2,NV)
# w3_2  = real.(H56[4][ind3,n3_2])
# 
# # σ'^2*|z5|^2
# i3_3  = [2; 2; ind3; ind4] .- 1
# n3_3  = CenterManifold.GetIndexNumber(i3_3,NV)
# w3_3  = real.(H56[4][ind3,n3_3])
# 
# # LHS
# #------------------------------  
# # |z5|^2*|z5|^2
# i3_4  = [ind3; ind4; ind3; ind4] .- 1
# n3_4  = CenterManifold.GetIndexNumber(i3_4,NV)
# w3_4  = real.(H56[4][ind3,n3_4])
# L5_56 = -w3_4
# 
# # |z3|^2*|z5|^2
# i3_5  = [ind1; ind2; ind3; ind4] .- 1
# n3_5  = CenterManifold.GetIndexNumber(i3_5,NV)
# w3_5  = real.(H56[4][ind3,n3_5])
# L5_34 = -w3_5

# amplitude_z5(z3,z5,ϵ,σ) = w1_1 + w2_1*ϵ + w2_2*σ + w3_1*ϵ^2 + w3_2*ϵ*σ + w3_3*σ^2 - (z3^2)*L5_34 - (z5^2)*L5_56
# nullc_z5(ϵ,σ) = w1_1 + w2_1*ϵ + w2_2*σ + w3_1*ϵ^2 + w3_2*ϵ*σ + w3_3*σ^2
# nullc_z5e(ϵ)  = nullc_z5(ϵ,σ0)

L5_56,L5_34,amplitude_z5 = GetAmpNullClines(H56,1,2,ind3,ind4,ind1,ind2)
nullc_z5(ϵ,σ)  = amplitude_z5(0.0,0.0,ϵ,σ)
nullc_z5e(ϵ)   = amplitude_z5(0.0,0.0,ϵ,σ0)

# Plotting
#----------------------------------------------------------------------

close("all")

nullclinesave = true
pdir          = "plots/"

h2,(ax2,ax3,ax4) = subplots(1,3,sharey=true,figsize=[16,5],gridspec_kw=Dict("width_ratios"=>[1, 1, 1]),layout="constrained",clear=true)

xlims = [-0.0005; 0.020]
ylims = [-0.0005; 0.020]
xtks  = Vector(0:5)*0.004
ytks  = Vector(0:5)*0.004

ax2.set_xlim(xlims)
ax2.set_ylim(ylims)
ax3.set_xlim(xlims)
ax3.set_ylim(ylims)
ax4.set_xlim(xlims)
ax4.set_ylim(ylims)

ax2.set_xticks(xtks)
ax2.set_yticks(ytks)
ax3.set_xticks(xtks)
ax3.set_yticks(ytks)
ax4.set_xticks(xtks)
ax4.set_yticks(ytks)

ax2.set_ylabel(L"|z_{3}|",fontsize=lafs)


nstep       = 10000
z3_zeros    = zeros(Float64,nstep)
z5_zeros    = zeros(Float64,nstep)
z3_line     = LinRange(0.0,0.025,nstep)
z5_line     = LinRange(0.0,0.025,nstep)

#h2    = figure(num=2,figsize=[7.5,6.5])
#h2.clf()
#ax2   = h2.subplots()
# ax2.legend(fontsize=lgfs)
# ax2.set_xlabel(L"|z_{1}|",fontsize=lafs)
# ax2.set_ylabel(L"|z_{3}|",fontsize=lafs)
# ax2.set_xlim(xlims)
# ax2.set_ylim(ylims)

include("second_bifurcation.jl")
include("second_third_bifurcation.jl")
include("third_bifurcation.jl")

h3,(ax5,ax6,ax7) = subplots(1,3,sharey=true,figsize=[16,5],gridspec_kw=Dict("width_ratios"=>[1, 1, 1]),layout="constrained",clear=true)

ax5.set_ylabel(L"|z_{3}|",fontsize=lafs)

ax5.set_xlim(xlims)
ax5.set_ylim(ylims)
ax6.set_xlim(xlims)
ax6.set_ylim(ylims)
ax7.set_xlim(xlims)
ax7.set_ylim(ylims)

ax5.set_xticks(xtks)
ax5.set_yticks(ytks)
ax6.set_xticks(xtks)
ax6.set_yticks(ytks)
ax7.set_xticks(xtks)
ax7.set_yticks(ytks)

#include("third_fourth_bifurcation.jl")
include("re4600_bifurcation.jl")
include("fourth_bifurcation.jl")
include("fourth_fifth_bifurcation.jl")
#include("fifth_bifurcation.jl")

println("Done.")




