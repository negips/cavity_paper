include("ReadAmplitudeFile.jl")
include("$(JULIAHOME)/center_manifold/Module_CenterManifold/CenterManifold.jl")

using PyPlot
using Peaks
using Statistics
using Printf
using Random

include("$JULIACOMMON/SetZero.jl")

#---------------------------------------------------------------------- 
function GetReynolds(ϵ::Float64,Rec::Float64)

  Reci = 1.0/Rec
  Rey  = Rec/(1.0 - ϵ)

  return Rey
end  
#---------------------------------------------------------------------- 
function GetEpsilon(Rey::Float64,Rec::Float64)

  ϵ    = 1.0 - Rec/Rey

  return ϵ
end  
#---------------------------------------------------------------------- 

lafs   = 16
lgfs   = 12
lw     = 3        # Linewidth
mew    = 3        # markeredgewidth
mks    = 10       # markersize

Rec    = 4.131331106555E+03
Reci   = 1.0/Rec

FOL    = "cavity_normal_form_pert/"
#FOL    = "cavity_normal_form_pert_O3/"
nmodes = 6
fnames = Vector{String}(undef,nmodes)
for i in 1:nmodes
  fnames[i] = @sprintf("%s%s%.2d",FOL,"AmplitudeTerms.",i-1)
end  

Coeffs = Vector{Matrix{ComplexF64}}(undef,nmodes)
for i in 1:nmodes
  AMat      = ReadAmplitudeFile(fnames[i])
  Coeffs[i] = copy(AMat)
end  

Tmax,OrdM = size(Coeffs[1])
Nt1       = CenterManifold.NInteractionTerms(1,nmodes)
Nt2       = CenterManifold.NInteractionTerms(2,nmodes)
Nt3       = CenterManifold.NInteractionTerms(3,nmodes)

G1        = zeros(ComplexF64,nmodes,Nt1)
G2        = zeros(ComplexF64,nmodes,Nt2)
G3        = zeros(ComplexF64,nmodes,Nt3)
G         = Vector{Matrix{ComplexF64}}(undef,OrdM)

for ord in 1:OrdM
  nt = CenterManifold.NInteractionTerms(ord,nmodes)
  GO = zeros(ComplexF64,nmodes,nt)
  for i in 1:nmodes
    GO[i,:] = copy(Coeffs[i][1:nt,ord])
  end
  G[ord] = copy(GO)
end

println("AmplitudeTerms.* Reading Done.")

