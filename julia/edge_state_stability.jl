# Eigenvalue of the Quasi-Steady State
#----------------------------------------------------------------------

# z5 is my z_{1}
# z3 is my z_{3}

# Start of Quasi-Steady solutions
#--------------------------------------------------  
nullc_z35e(ϵ) = nullc_z5e(ϵ) - L5_34/L3_34*nullc_z3e(ϵ)
QSϵ1          = find_zero(nullc_z35e,0.1)

# End of Quasi-steady solutions
#--------------------------------------------------  
nullc_z35e(ϵ) = nullc_z3e(ϵ) - L3_56/L5_56*nullc_z5e(ϵ)
QSϵ2          = find_zero(nullc_z35e,0.1) 

ϵrange = LinRange(QSϵ1,QSϵ2,nstep)
QSz3   = zeros(Float64,nstep)
QSz5   = zeros(Float64,nstep)

# Edge-States
for i in 1:nstep
  δ         = ϵrange[i]
  lhs2      = L3_34 - L3_56*L5_34/L5_56
  rhs2      = nullc_z3e(δ) - L3_56/L5_56*nullc_z5e(δ)
  ratio     = rhs2/lhs2
  if ratio<0
    ratio = 0.0
  end
  QSz3[i]   = sqrt(ratio)
  tmp3      = (nullc_z5e(δ) - L5_34*(QSz3[i]^2))/L5_56
  if tmp3 <0
    tmp3 = 0.0
  end
  QSz5[i]   = sqrt(tmp3)
end

close("all")
h10   = figure(num=10,figsize=[7.5,6])
ax10  = h10.subplots()

ax10.plot(ϵrange,QSz5,linestyle="-",color=cmap(0))
ax10.plot(ϵrange,QSz3,linestyle="-",color=cmap(1))
ax10.set_xlabel(L"ϵ",fontsize=lafs)
ax10.set_ylabel(L"|z_{i}|",fontsize=lafs)


# Build Eigenvalue Matrices

jj     = 10
ϵtmp   = ϵrange[jj]
Lmat1  = zeros(Float64,2,2)
Lmat1[1,1]  = nullc_z5e(ϵtmp) - 2*L5_56*QSz5[jj]^2
Lmat1[1,2]  = -L5_34*QSz3[jj]^2
Lmat1[2,1]  = -L3_56*QSz5[jj]^2
Lmat1[2,2]  = nullc_z3e(ϵtmp) - 2*L3_34*QSz3[jj]^2

EGVMat  = zeros(Float64,2,2,nstep)
EIGMat  = zeros(Float64,2,nstep)
for j in 1:nstep
  δ2        = ϵrange[j]
  Lmat      = zeros(Float64,2,2)
  Lmat[1,1] = nullc_z5e(δ2) - 2*L5_56*QSz5[j]^2
  Lmat[1,2] = -L5_34*QSz3[j]^2
  Lmat[2,1] = -L3_56*QSz5[j]^2
  Lmat[2,2] = nullc_z3e(δ2) - 2*L3_34*QSz3[j]^2

  F               = eigen(Lmat)
  EIGMat[:,j]     = F.values
  EGVMat[:,:,j]   = F.vectors
end

h11   = figure(num=11,figsize=[7.5,6])
ax11  = h11.subplots()
ax11.plot(ϵrange,EIGMat[1,:],linestyle="-",color=cmap(0))
ax11.plot(ϵrange,EIGMat[2,:],linestyle="-",color=cmap(1))

h12   = figure(num=12,figsize=[7.5,6])
ax12  = h12.subplots()
ax12.semilogy(ϵrange,EGVMat[1,1,:]./EGVMat[2,1,:],linestyle="-",color=cmap(0))
ax12.semilogy(ϵrange,EIGMat[1,2,:]./EGVMat[2,2,:],linestyle="-",color=cmap(1))


#ax7.plot(z5_vals,z3_range,linestyle="-",linewidth=lw,color=cmap(1),label=L"Nullcline - |z_{3}|\ LCO")
#ax7.plot(z5_line,z3_zeros,linestyle="--",linewidth=lw,color=cmap(1))

#ax7.set_xlabel(L"|z_{1}|",fontsize=lafs)
#ax7.set_ylabel(L"|z_{3}|",fontsize=lafs)
#ax7.set_xlim(xlims)
#ax7.set_ylim(ylims)

#tt2 = @sprintf("Re=%.1f",Re5)
#ax7.set_title(tt2,fontsize=lafs)
#ax7.legend(fontsize=lgfs,ncols=1)

if (nullclinesave)
  svname = "bifurcation5.eps"
  # h3.savefig(pdir*svname)
end  

println("Done.")







