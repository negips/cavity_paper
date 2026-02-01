using Roots

#---------------------------------------------------------------------- 
function GetAmplitudeSq(ind::Int,G::Vector{Matrix{T}},Conj::Vector{Int};ifreal=false) where {T<:ComplexF64}

  nv   = length(Conj)
  OrdM = length(G)

  if (ifreal)
    H    = Vector{Matrix{Float64}}(undef,OrdM+1)
    H1   = zeros(Float64,1,nv)
    H[1] = H1
  else
    H    = Vector{Matrix{T}}(undef,OrdM+1)
    H1   = zeros(T,1,nv)
    H[1] = H1
  end  

  # Assuming we do have conjugates
  i1 = ind
  i2 = Conj[i1]
  for ord in 1:OrdM
    
    ord1 = ord+1
    nt1  = CenterManifold.NInteractionTerms(ord1,nv)
    Hi   = zeros(T,1,nt1)
    # Hc   = zeros(T,1,nt1)

    h1   = zeros(T,nt1)
    h2   = zeros(T,nt1)

    nt   = CenterManifold.NInteractionTerms(ord,nv)
    
    g1   = G[ord][i1,:]
    g2   = G[ord][i2,:]

    ind1 = zeros(Int,ord+1)
    ind2 = zeros(Int,ord+1)
   
    for i in 1:nt
      ind0        = CenterManifold.GetInteractionTerms(i,ord,nv)

      copyto!(ind1,1,ind0,1,ord)
      ind1[ord1]  = i2-1
      n1          = CenterManifold.GetIndexNumber(ind1,nv)
      val1        = g1[i]
      h1[n1]      = val1

      copyto!(ind2,1,ind0,1,ord)
      ind2[ord1]  = i1-1
      n2          = CenterManifold.GetIndexNumber(ind2,nv)
      val2        = g2[i]
      h2[n2]      = val2
    end
    h1            = h1 .+ h2

    if ifreal
      Hi[1,:]       = copy(real.(h1))
    else  
      Hi[1,:]       = copy(h1)
    end  
    # Hc[1,:]       = h1

    H[ord1]       = copy(Hi)
  end    


  return H
end
#----------------------------------------------------------------------
function GetAmplitudeAbs(ind::Int,G::Vector{Matrix{T}},Conj::Vector{Int}) where {T<:ComplexF64}

  # Remove one |z|^2 from the amplitude equations

  nv   = length(Conj)
  OrdM = length(G)

  # Coefficients for |z|^2
  Hsq  = GetAmplitudeSq(ind,G,Conj)

  i1   = ind
  i2   = Conj[i1]
  ind0 = [i1; i2] .- 1
  n0   = CenterManifold.GetIndexNumber(ind0,nv)
  h0   = Hsq[2][1,n0]         # Only Term at zeroth order.

  return h0,Hsq
end
# #----------------------------------------------------------------------
# function GetReducedAmpDT(x,ind,indc,ϵ,σ,Hi::Vector{Matrix{T}},nv::Int) where {T}
# 
#     z       = zeros(Float64,nv)
#     z[1]    = ϵ
#     z[2]    = σ
#     z[ind]  = x
#     if (indc>0)
#       z[indc] = x
#     end  
#     dx      = CenterManifold.GetAsymptoticField(z,Hi)[1]
# 
#     return dx
# end  
# #---------------------------------------------------------------------- 
# function GetEquilibriumAmp(Hi,ϵ,σ,ind,indc,nv)
# 
#   f(x)      = GetReducedAmpDT(x,ind,indc,ϵ,σ,Hi,nv)
# 
#   xi        = 1.0e-6
#   xe        = 0.1
#   x_zeros   = find_zeros(f,xi,xe)
# 
#   return x_zeros
# end
# #---------------------------------------------------------------------- 
function GetReducedAmpDT(x,ind,indc,x0,Hi::Vector{Matrix{T}},nv::Int) where {T}

    z       = copy(x0)
    z[ind]  = x
    if (indc>0)
      z[indc] = x
    end  
    dx      = CenterManifold.GetAsymptoticField(z,Hi)[1]

    return dx
end  
#---------------------------------------------------------------------- 
function GetEquilibriumAmp(Hi,ind,indc,z0,nv)

  f(x)      = GetReducedAmpDT(x,ind,indc,z0,Hi,nv)

  xi        = 1.0e-6
  xe        = 0.1
  x_zeros   = find_zeros(f,xi, xe)

  return x_zeros
end
#---------------------------------------------------------------------- 

include("read_amplitudes.jl")

#-------------------------------------------------- 
Conj = [1; 2; 4; 3; 6; 5]     # Conjugate Indices
SInd = [5; 6; 3; 4; 1; 2]     # Swapped Indices

#σ0   = -9.1854401740e-02
σ0   = -9.2e-02

ind  = 3
ind2 = Conj[ind]
#CenterManifold.DisplayTerms(ind,"z",G)
#CenterManifold.DisplayTermsSwapped(ind,"z",G,SInd)

#CenterManifold.DisplayTerms(ind2,"z",G)

ind  = 3
ind2 = Conj[ind]

H34 = GetAmplitudeSq(ind,G,Conj,ifreal=true)
CenterManifold.DisplayTerms(1,"z",H34)
#CenterManifold.DisplayTermsSwapped(1,"z",H34,SInd)

H34_0,H34_1 = GetAmplitudeAbs(ind,G,Conj)

ind  = 5
H56 = GetAmplitudeSq(ind,G,Conj,ifreal=true)
#CenterManifold.DisplayTerms(1,"z",H56)
#CenterManifold.DisplayTermsSwapped(1,"z",H56,SInd)

H56_0,H56_1 = GetAmplitudeAbs(ind,G,Conj)

#---------------------------------------------------------------------- 

# Positions of values
ind1  = 3
ind2  = Conj[ind1]
ind3  = 5
ind4  = Conj[ind3]

NV    = 6         # No of variables
NV0   = 2         # No. of parameter variables

ϵind   = 1
σind   = 2
NV0Ind = [ϵind; σind]

#ϵ       = 0.01
# r       = 4349.6
# ϵ       = GetEpsilon(r,Rec)
# σ       = -9.2e-02
# z0      = zeros(Float64,NV)
# z0[1]   = ϵ
# z0[2]   = σ
# z_zeros = GetEquilibriumAmp(H56,ind3,ind4,z0,NV)
#---------------------------------------------------------------------- 

# Plot Amplitudes/Reynolds

npts    = 1000
ϵrange  = LinRange(1.0e-6,0.20,npts)
Rerange = zeros(Float64,length(ϵrange))
EqAmp56 = zeros(Float64,length(ϵrange))
EqAmp34 = zeros(Float64,length(ϵrange))

for i in 1:length(ϵrange)

  local σ       = -9.2e-02
  local z0      = zeros(Float64,NV)
  local z0[1]   = ϵrange[i]
  local z0[2]   = σ
  Rerange[i]    = GetReynolds(ϵrange[i],Rec)
  zzeros        = GetEquilibriumAmp(H56,ind3,ind4,z0,NV)
  if !isempty(zzeros)
    EqAmp56[i]  = zzeros[1]
  else
    EqAmp56[i]  = -1.0e-6
  end

end  


#close("all")

h2,(ax2) = subplots(1,1,sharey=true,figsize=[6,5],gridspec_kw=Dict("width_ratios"=>[1]),layout="constrained",clear=true)

xlims = [-0.0005; 0.020]
ylims = [-0.0005; 0.020]
# xtks  = Vector(0:5)*0.004
# ytks  = Vector(0:5)*0.004
#ax2.set_xlim(xlims)
#ax2.set_ylim(ylims)

ax2.plot(Rerange,EqAmp56,linewidth=2)


for i in 1:length(ϵrange)

  local σ       = -9.2e-02
  local z0      = zeros(Float64,NV)
  local z0[1]   = ϵrange[i]
  local z0[2]   = σ
  # local z0[ind3]= EqAmp56[i]
  # local z0[ind4]= EqAmp56[i]
 
  Rerange[i]    = GetReynolds(ϵrange[i],Rec)
  zzeros        = GetEquilibriumAmp(H34,ind1,ind2,z0,NV)
  if !isempty(zzeros)
    EqAmp34[i]  = zzeros[1]
  else
    EqAmp34[i]  = -1.0e-6
  end
end  

ax2.plot(Rerange,EqAmp34,linewidth=2)

for i in 1:length(ϵrange)

  local σ       = -9.2e-02
  local z0      = zeros(Float64,NV)
  local z0[1]   = ϵrange[i]
  local z0[2]   = σ
  local z0[ind3]= EqAmp56[i]
  local z0[ind4]= EqAmp56[i]
 
  Rerange[i]    = GetReynolds(ϵrange[i],Rec)
  zzeros        = GetEquilibriumAmp(H34,ind1,ind2,z0,NV)
  if !isempty(zzeros)
    EqAmp34[i]  = zzeros[1]
  else
    EqAmp34[i]  = -1.0e-6
  end
end  

ax2.plot(Rerange,EqAmp34,linewidth=2)

#---------------------------------------------------------------------- 
println("NullCline function Generated.")







