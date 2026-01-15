#println("Time-Stepping for the reduced Center-Manifold Approximation")

include("read_amplitudes.jl")

#---------------------------------------------------------------------- 
function LinearizedOP(v::Vector{T},v0::Vector{T},G::Vector{Matrix{T}}) where {T<:Number}

  OrdM = length(G)
  nv   = size(G[1],2)
  dv   = zeros(T,nv)

  # Linear Part
  ord  = 1
  n    = CenterManifold.NInteractionTerms(ord,nv)
  dv   = G[ord]*v

  # Second Order
  if (OrdM>=2)
    ord = 2
    n    = CenterManifold.NInteractionTerms(ord,nv)
    I2   = CenterManifold.AllInteractionIndices(ord,nv) .+ 1
    for i in 1:nv
      for j in 1:n
        j1    = I2[j,1]
        j2    = I2[j,2]
        dv[i] = dv[i] + G[ord][i,j]*(v0[j1]*v[j2] + v[j1]*v0[j2])
      end
    end
  end  
  
  # Third Order
  if (OrdM>=3)
    ord = 3
    n   = CenterManifold.NInteractionTerms(ord,nv)
    I3  = CenterManifold.AllInteractionIndices(ord,nv) .+ 1
    for i in 1:nv
      for j in 1:n
        j1    = I3[j,1]
        j2    = I3[j,2]
        j3    = I3[j,3]
        dv[i] = dv[i] + G[ord][i,j]*(v0[j1]*v0[j2]*v[j3] + v0[j1]*v[j2]*v0[j3] + v[j1]*v0[j2]*v0[j3])
      end
    end
  end  

  return dv
end
#---------------------------------------------------------------------- 
# 4th Order Runge-Kutta Steps for time dependent Operator
function FloquetRK4!(v::AbstractVector{T},G::Vector{Matrix{S}},v0::Vector{T},Ω::Vector{ComplexF64},time::Float64,dt::Float64) where {T,S<:Number}

  # localprec = eltype(v[1])
  two = T(2)
  six = T(6)

  ph1   = exp.(Ω*time)
  Vt1   = v0.*ph1
  k1    = LinearizedOP(v,Vt1,G)
  v1    = v .+ dt/two*k1

  ph2   = exp.(Ω*(time+dt/2))
  Vt2   = v0.*ph2
  k2    = LinearizedOP(v1,Vt2,G)
  v2    = v .+ dt/two*k2

  ph3   = exp.(Ω*(time+dt/2))
  Vt3   = v0.*ph3
  k3    = LinearizedOP(v2,Vt3,G)
  v3    = v .+ dt*k3

  ph4   = exp.(Ω*(time+dt))
  Vt4   = v0.*ph4
  k4    = LinearizedOP(v3,Vt4,G)
 
  v    .= v .+ dt/six*(k1 .+ two*k2 .+ two*k3 .+ k4)

  return nothing    
end  
#----------------------------------------------------------------------

using LinearAlgebra
include("$(JULIAHOME)/center_manifold/Module_StepperArnoldi/StepperArnoldi.jl")

include("FloquetStepper.jl")

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
A1eq,Ω1eq   = GetEquilibriumAmpOmega(G,VBase,ind1,ind1c)
@printf("Eq. |z1|: %.5f, Eq. Ω1: %.5f\n",A1eq,Ω1eq)

ind2        = 3
ind2c       = Conj[ind2]
A2eq,Ω2eq   = GetEquilibriumAmpOmega(G,VBase,ind2,ind2c)
@printf("Eq. |z2|: %.5f, Eq. Ω2: %.5f\n",A2eq,Ω2eq)

VB          = copy(VBase)
VB[ind2]    = A2eq
VB[ind2c]   = A2eq
λe          = GetEffectiveEigenvalue(G,VB,ind1,ind1c)
@printf("λ1,2 Effective - λr: %.5f, λi: %.5f\n",real(λe),imag(λe))

VB          = copy(VBase)
VB[ind1]    = A1eq
VB[ind1c]   = A1eq
λe          = GetEffectiveEigenvalue(G,VB,ind2,ind2c)
@printf("λ3,4 Effective - λr: %.5f, λi: %.5f\n",real(λe),imag(λe))

#---------------------------------------------------------------------- 

# Perform Floquet Analysis
#---------------------------------------------------------------------- 

Period      = 2*π/Ω2eq
nsteps      = 5000
dt          = Period/nsteps
Ω           = zeros(ComplexF64,nmodes)

# Ω[ind1]     =  Ω1eq*im
# Ω[ind1c]    = -Ω1eq*im
# VB          = copy(VBase)
# VB[ind1]    = A1eq
# VB[ind1c]   = A1eq

Ω[ind2]     =  Ω2eq*im
Ω[ind2c]    = -Ω2eq*im
VB          = copy(VBase)
VB[ind2]    = A2eq
VB[ind2c]   = A2eq

# Stepper-Arnoldi
#-------------------------------------------------- 
ifadjoint         = false
ifoptimal         = false
ifverbose         = false
verbosestep       = 500
#nsteps            = 500
#dt                = 1.0e-4
StpInp            = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)

ifarnoldi         = true 
ifverbose         = false
ifeigshift        = false
vlen              = nmodes
nev               = nmodes
ekryl             = 0
lkryl             = nev + ekryl 
eigshift          = 0.0 + 1.0im
ngs               = 2
bsize             = 1
outer_iterations  = 2
tol               = 1.0e-12
ArnInp            = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,ifeigshift,vlen,nev,ekryl,lkryl,eigshift,ngs,bsize,outer_iterations,tol)

B                 = ones(Float64,nmodes)  # Weight
ArnDir            = FloquetStepper(G,B,VB,Ω,StpInp,ArnInp)

display(abs.(ArnDir.ritz))

@printf "Done."
















