#println("Time-Stepping for the reduced Center-Manifold Approximation")

include("ReadAmplitudeFile.jl")
#include("../../Module_CenterManifold/CenterManifold.jl")
include("$(JULIAHOME)/center_manifold/Module_CenterManifold/CenterManifold.jl")

using PyPlot
using Peaks
using Statistics
using Printf
using Random

include("$JULIACOMMON/SetZero.jl")

lafs   = 16
lgfs   = 12
lw     = 3        # Linewidth
mew    = 3        # markeredgewidth
mks    = 10       # markersize

Rec    = 4.131331106555E+03
Reci   = 1.0/Rec

close("all")
FOL    = "cavity_normal_form_pert/"

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

# Mode numbering starts from 0.
# A1  - Modes 4/5
# A2  - Modes 2/3

# For |A1|²
#---------------------------------------- 
λ4     = Coeffs[5][5,1]
λ5     = Coeffs[6][6,1]
λ45    = real(λ4 + λ5)

# ν|A1|²
α4_04  = Coeffs[5][5,2] 
α5_05  = Coeffs[6][6,2]
α45_0  = real(α4_04 + α5_05)

# σ|A1|²
α4_14  = Coeffs[5][10,2] 
α5_15  = Coeffs[6][11,2] 
α45_1  = real(α4_14 + α5_15)

# ν²|A1|²
β4_004 = Coeffs[5][5,3] 
β5_005 = Coeffs[6][6,3] 
β45_00 = real(β4_004 + β5_005)

# νσ|A1|²
β4_014 = Coeffs[5][10,3] 
β5_015 = Coeffs[6][11,3] 
β45_01 = real(β4_014 + β5_015)

# σ²|A1|²
β4_114 = Coeffs[5][25,3] 
β5_115 = Coeffs[6][26,3] 
β45_11 = real(β4_114 + β5_115)

# |A2|² ⋅ |A1|²
β4_234 = Coeffs[5][42,3] 
β5_235 = Coeffs[6][43,3] 
β45_23 = real(β4_234 + β5_235)

# |A1|² ⋅ |A1|²
β4_445 = Coeffs[5][54,3] 
β5_455 = Coeffs[6][55,3] 
β45_45 = real(β4_445 + β5_455)
#---------------------------------------- 

# For |A2|²
#---------------------------------------- 
λ2     = Coeffs[3][3,1]
λ3     = Coeffs[4][4,1]
λ23    = real(λ2 + λ3)

# ν|A2|²
α2_02  = Coeffs[3][3,2] 
α3_03  = Coeffs[4][4,2] 
α23_0  = real(α2_02 + α3_03)

# σ|A2|²
α2_12  = Coeffs[3][8,2] 
α3_13  = Coeffs[4][9,2] 
α23_1  = real(α2_12 + α3_13)

# ν²|A2|²
β2_002 = Coeffs[3][3,3] 
β3_003 = Coeffs[4][4,3] 
β23_00 = real(β2_002 + β3_003)

# νσ|A2|²
β2_012 = Coeffs[3][8,3] 
β3_013 = Coeffs[4][9,3] 
β23_01 = real(β2_012 + β3_013)

# σ²|A1|²
β2_112 = Coeffs[3][23,3] 
β3_113 = Coeffs[4][24,3] 
β23_11 = real(β2_112 + β3_113)

# |A2|² ⋅ |A2|²
β2_223 = Coeffs[3][38,3] 
β3_233 = Coeffs[4][41,3] 
β23_23 = real(β2_223 + β3_233)

# |A1|² ⋅ |A2|²
β2_245 = Coeffs[3][45,3] 
β3_345 = Coeffs[4][51,3] 
β23_45 = real(β2_245 + β3_345)

#---------------------------------------- 

ν      = 0.05017445479314724 + 0.00
#ν      =  0.08
#σ      = -9.1854401740e-02
σ      = -0.092

Amin        = 0.0
Amax        = 0.012
nstep       = 80000
A1A2        = LinRange(Amin,Amax,nstep)

close("all")

cmap  = get_cmap("tab10")

h1    = figure(num=1,figsize=[8,8])
ax1   = h1.subplots()

# Nullcline for A1 - Modes 4/5
A2                = A1A2
A1_null           = (-(λ45 + α45_0*ν + α45_1*σ + β45_00*ν*ν + β45_01*ν*σ + β45_11*σ*σ) .- β45_23.*A2.^2)./β45_45
nz_ind            = A1_null .> 0.0
tmp               = A1_null[nz_ind]
l                 = length(tmp)
A1_null2          = zeros(Float64,l+1)
A2_2              = zeros(Float64,l+1)
A1_null2[1:l]     = A1_null[nz_ind]
A2_2[1:l]         = A2[nz_ind]
# Last point      
A1_last           = 0.0
A2_last           = -(λ45 + α45_0*ν + α45_1*σ + β45_00*ν*ν + β45_01*ν*σ + β45_11*σ*σ)/β45_23
A1_null2[l+1]     = A1_last
A2_2[l+1]         = sqrt(A2_last)
A1_null2          = sqrt.(A1_null2)
A1_zeros          = zeros(Float64,length(A2)) #.+ 1.0e-7

ax1.plot(A1_null2,A2_2,linestyle="-",linewidth=lw,color=cmap(0),label="Nullcline |A1|")
ax1.plot(A1_zeros,A2,  linestyle="-",linewidth=lw,color=cmap(0))


# Nullcline for A2 - Modes 2/3
A1                = A1A2
A2_null           = (-(λ23 + α23_0*ν + α23_1*σ + β23_00*ν*ν + β23_01*ν*σ + β23_11*σ*σ) .- β23_45.*A1.^2)./β23_23
nz_ind            = A2_null .> 0.0
tmp               = A2_null[nz_ind]
l                 = length(tmp)
A2_null2          = zeros(Float64,l+1)
A1_2              = zeros(Float64,l+1)
A2_null2[1:l]     = A2_null[nz_ind]
A1_2[1:l]         = A1[nz_ind]
# Last point      
A2_last           = 0.0
A1_last           = -(λ23 + α23_0*ν + α23_1*σ + β23_00*ν*ν + β23_01*ν*σ + β23_11*σ*σ)/β23_45
A2_null2[l+1]     = A2_last
A1_2[l+1]         = sqrt(A1_last)
A2_null2          = sqrt.(A2_null2)
A2_zeros          = zeros(Float64,length(A1)) #.+ 1.0e-7
 
 
ax1.plot(A1_2,A2_null2,linestyle="-",linewidth=lw,color=cmap(1),label="Nullcline |A2|")
ax1.plot(A1,  A2_zeros,linestyle="-",linewidth=lw,color=cmap(1))

ax1.set_xlim([Amin-5.0e-4,Amax])
ax1.set_ylim([Amin-5.0e-4,Amax])
ax1.set_xlabel(L"|A1|",fontsize=lafs)
ax1.set_ylabel(L"|A2|",fontsize=lafs)

legend(fontsize=lgfs)

Rei    = Reci*(1.0 - ν)
Re     = 1.0/Rei
@printf("ν: %g ; Re = %g ;\n", ν, Re)

# Fixed Points
# Mode A1 fixed pooint
FP1_A2       = 0.0
tmp2         = -(λ45 + α45_0*ν + α45_1*σ + β45_00*ν*ν + β45_01*ν*σ + β45_11*σ*σ)/β45_45
if tmp2 >= 0.0
  FP1_A1       = sqrt(tmp2)
  ax1.plot(FP1_A1,FP1_A2,linestyle="none",marker="o",fillstyle="full",markeredgewidth=mew,markersize=mks,color="black")
end  

# Mode A2 fixed point
FP2_A1       = 0.0
tmp2         = -(λ23 + α23_0*ν + α23_1*σ + β23_00*ν*ν + β23_01*ν*σ + β23_11*σ*σ)/β23_23
if tmp2 >= 0.0
  FP2_A2       = sqrt(tmp2)
  ax1.plot(FP2_A1,FP2_A2,linestyle="none",marker="o",fillstyle="full",markeredgewidth=mew,markersize=mks,color="black")
end  

# Saddle node

c2           =  (λ23 + α23_0*ν + α23_1*σ + β23_00*ν*ν + β23_01*ν*σ + β23_11*σ*σ)
c1           =  (λ45 + α45_0*ν + α45_1*σ + β45_00*ν*ν + β45_01*ν*σ + β45_11*σ*σ)
c11          =  c1 - c2*β45_23/β23_23
b11          =  β45_45 - β45_23*β23_45/β23_23
A1_sq        = -c11/b11 
A2_sq        = -(c1 + β45_45*A1_sq)/β45_23

if (A1_sq >=0.0) && (A2_sq >= 0.0)
  FP3_A1     = sqrt(A1_sq)
  FP3_A2     = sqrt(A2_sq)
  ax1.plot(FP3_A1,FP3_A2,linestyle="none",marker="o",fillstyle="none",markeredgewidth=mew,markersize=mks,color="black")
  @printf("Edge State: (%6g,%6g)\n",FP3_A1,FP3_A2)
else
  @printf "No Edge state.\n"
end  


@printf "Done."
















