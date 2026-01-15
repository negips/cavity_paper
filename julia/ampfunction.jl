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
include("read_amplitudes.jl")

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

# |z3|^2
#-------------------------------------------------- 
# 1st order
# |z3|^2
i1_1  = [ind1; ind2] .- 1
n1_1  = CenterManifold.GetIndexNumber(i1_1,NV)
v1_1  = real.(H34[2][ind1,n1_1])


# 2nd order
# ϵ*|z3|^2
i2_1  = [1; ind1; ind2] .- 1
n2_1  = CenterManifold.GetIndexNumber(i2_1,NV)
v2_1  = real.(H34[3][ind1,n2_1])

# σ'*|z3|^2
i2_2  = [2; ind1; ind2] .- 1
n2_2  = CenterManifold.GetIndexNumber(i2_2,NV)
v2_2  = real.(H34[3][ind1,n2_2])

# 3rd order
# ϵ^2*|z3|^2
i3_1  = [1; 1; ind1; ind2] .- 1
n3_1  = CenterManifold.GetIndexNumber(i3_1,NV)
v3_1  = real.(H34[4][ind1,n3_1])

# ϵ*σ'*|z3|^2
i3_2  = [1; 2; ind1; ind2] .- 1
n3_2  = CenterManifold.GetIndexNumber(i3_2,NV)
v3_2  = real.(H34[4][ind1,n3_2])

# σ'^2*|z3|^2
i3_3  = [2; 2; ind1; ind2] .- 1
n3_3  = CenterManifold.GetIndexNumber(i3_3,NV)
v3_3  = real.(H34[4][ind1,n3_3])

# LHS
#------------------------------  
# |z3|^2*|z3|^2
i3_4  = [ind1; ind2; ind1; ind2] .- 1
n3_4  = CenterManifold.GetIndexNumber(i3_4,NV)
v3_4  = real.(H34[4][ind1,n3_4])
L3_34 = -v3_4

# |z5|^2*|z3|^2
i3_5  = [ind3; ind4; ind1; ind2] .- 1
n3_5  = CenterManifold.GetIndexNumber(i3_5,NV)
v3_5  = real.(H34[4][ind1,n3_5])
L3_56 = -v3_5

amplitude_z3(z3,z5,ϵ,σ) = v1_1 + v2_1*ϵ + v2_2*σ + v3_1*ϵ^2 + v3_2*ϵ*σ + v3_3*σ^2 - (z3^2)*L3_34 - (z5^2)*L3_56
nullc_z3(ϵ,σ) = amplitude_z3(0.0,0.0,ϵ,σ)
nullc_z3e(ϵ)  = amplitude_z3(0.0,0.0,ϵ,σ0)


# |z5|^2
#-------------------------------------------------- 
# 1st order
# |z5|^2
i1_1  = [ind3; ind4] .- 1
n1_1  = CenterManifold.GetIndexNumber(i1_1,NV)
w1_1  = real.(H56[2][ind3,n1_1])

# 2nd order
# ϵ*|z5|^2
i2_1  = [1; ind3; ind4] .- 1
n2_1  = CenterManifold.GetIndexNumber(i2_1,NV)
w2_1  = real.(H56[3][ind3,n2_1])

# σ'*|z5|^2
i2_2  = [2; ind3; ind4] .- 1
n2_2  = CenterManifold.GetIndexNumber(i2_2,NV)
w2_2  = real.(H56[3][ind3,n2_2])

# 3rd order
# ϵ^2*|z5|^2
i3_1  = [1; 1; ind3; ind4] .- 1
n3_1  = CenterManifold.GetIndexNumber(i3_1,NV)
w3_1  = real.(H56[4][ind3,n3_1])

# ϵ*σ'*|z5|^2
i3_2  = [1; 2; ind3; ind4] .- 1
n3_2  = CenterManifold.GetIndexNumber(i3_2,NV)
w3_2  = real.(H56[4][ind3,n3_2])

# σ'^2*|z5|^2
i3_3  = [2; 2; ind3; ind4] .- 1
n3_3  = CenterManifold.GetIndexNumber(i3_3,NV)
w3_3  = real.(H56[4][ind3,n3_3])

# LHS
#------------------------------  
# |z5|^2*|z5|^2
i3_4  = [ind3; ind4; ind3; ind4] .- 1
n3_4  = CenterManifold.GetIndexNumber(i3_4,NV)
w3_4  = real.(H56[4][ind3,n3_4])
L5_56 = -w3_4

# |z3|^2*|z5|^2
i3_5  = [ind1; ind2; ind3; ind4] .- 1
n3_5  = CenterManifold.GetIndexNumber(i3_5,NV)
w3_5  = real.(H56[4][ind3,n3_5])
L5_34 = -w3_5

amplitude_z5(z3,z5,ϵ,σ) = w1_1 + w2_1*ϵ + w2_2*σ + w3_1*ϵ^2 + w3_2*ϵ*σ + w3_3*σ^2 - (z3^2)*L5_34 - (z5^2)*L5_56
nullc_z5(ϵ,σ) = amplitude_z5(0.0,0.0,ϵ,σ)
nullc_z5e(ϵ)   = amplitude_z5(0.0,0.0,ϵ,σ0)

# Fixed Points
#----------------------------------------------------------------------
function GetFixedPoint(ampfunc,ϵval,σval,ax)

  fp  = zeros(Float64,2)

  nullc_ϵ = ampfunc(0.0,0.0,ϵval,σval)
  if ax == 1
    x1 = 1.0
    x2 = 0.0
  else
    x1 = 0.0
    x2 = 1.0
  end

  L         = - ampfunc(x1,x2,ϵval,σval)
  nullc_ϵ   = ampfunc(0.0,0.0,ϵval,σval)
  r         = nullc_ϵ/L

  if (r<0 && abs(r)<1.0e-16)
    r = 0.0
  end
  fp[ax]    =   sqrt(r)

  return fp
end
#----------------------------------------------------------------------
function GetEquilibriumAmpOmega(G,V0,ind,indc)

  # Only works for Normal Form representation
  # For Third Order

  v       = copy(V0)
  v[ind]  = 1.0
  v[indc] = 1.0

  nv      = length(v)
  OrdM    = length(G)

  if OrdM>=3
    dv1     = CenterManifold.GetAsymptoticField(v,G[1:3])

    index3  = [ind; ind; indc] .- 1
    n3      = CenterManifold.GetIndexNumber(index3,nv)
    α3      = G[3][ind,n3] 

    dv_ind   = dv1[ind] - α3

    amp2_eq  = -real(dv_ind)/real(α3)
    if (abs(amp2_eq)<1.0e-16)
      amp_eq = 0.0
    else
      amp_eq = sqrt(amp2_eq)
    end
    
    omega_eq = imag(dv_ind) + imag(α3*amp2_eq)

  end

  return amp_eq,omega_eq
end
#---------------------------------------------------------------------- 
function GetEffectiveEigenvalue(G,V0,ind,indc)

  # Only works for Normal Form representation
  # For Third Order

  v       = copy(V0)
  v[ind]  = 1.0
  v[indc] = 1.0

  nv      = length(v)
  OrdM    = length(G)

  if OrdM>=3
    dv1     = CenterManifold.GetAsymptoticField(v,G[1:3])

    index3  = [ind; ind; indc] .- 1
    n3      = CenterManifold.GetIndexNumber(index3,nv)
    α3      = G[3][ind,n3] 

    λeff    = dv1[ind] - α3

    λeff    = λeff + α3*V0[ind]*V0[indc]
  end

  return λeff
end
#---------------------------------------------------------------------- 
function GetEdgeState(famp1,famp2,ϵval,σval)

  nullc1 = famp1(0.0,0.0,ϵval,σval)
  nullc2 = famp2(0.0,0.0,ϵval,σval)

  L1_1      = -famp1(1.0,0.0,0.0,0.0)
  L1_2      = -famp1(0.0,1.0,0.0,0.0)

  # println(L1_1, " ", L1_2)

  L2_1      = -famp2(1.0,0.0,0.0,0.0)
  L2_2      = -famp2(0.0,1.0,0.0,0.0)

  # println(L2_1, " ", L2_2)

  nullc1    = famp1(0.0,0.0,ϵval,σval)
  nullc2    = famp2(0.0,0.0,ϵval,σval)
     
  lhs1      = (L1_1 - L1_2*L2_1/L2_2)
  rhs1      = nullc1 - L1_2*nullc2/L2_2
  
  r         = rhs1/lhs1
  if (abs(r)<1.0e-16)
    r       = 0.0
  elseif r < -1.0e-16
    println("No Edge State")
    return []
  end

  x1        = sqrt(r)

  xx        = nullc2/L2_2 - L2_1/L2_2*x1^2
  if (abs(xx)<1.0e-16)
    xx       = 0.0
  elseif xx < -1.0e-16
    println("No Edge State")
    return []
  end

  x2  = sqrt(xx)

  return [x1;x2]
end
#---------------------------------------------------------------------- 

println("NullCline function Generated.")







