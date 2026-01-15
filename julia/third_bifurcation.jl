# Plotting
#----------------------------------------------------------------------

# Start of Quasi-Steady solutions
#--------------------------------------------------  
nullc_z35e(ϵ) = nullc_z5e(ϵ) - L5_34/L3_34*nullc_z3e(ϵ)
ϵ3          = find_zero(nullc_z35e,0.1)
Rei3        = Reci*(1.0 - ϵ3)
Re3         = 1.0/Rei3

# z5 is my z_{1}
# z3 is my z_{3}

z3_min      = 0.0
z3_max      = sqrt(nullc_z5e(ϵ3)/L5_34)

z3_range    = LinRange(z3_min,z3_max,nstep)
z5_vals     = zeros(Float64,nstep)
for i in 1:nstep
  if i != nstep
    z5_vals[i] = sqrt((nullc_z5e(ϵ3) - (z3_range[i]^2)*L5_34)/L5_56)
  else
    z5_vals[i] = 0.0
  end
end  

ax4.plot(z5_vals,z3_range,linestyle="-",linewidth=lw,color=cmap(0),label=L"Nullcline - |z_{1}|\ LCO")
#ax4.plot(z5_zeros,z3_line,linestyle="--",linewidth=lw,color=cmap(0),label=L"Nullcline - |z_{1}|\ Fixed\ Point")
ax4.plot(z5_zeros,z3_line,linestyle="--",linewidth=lw,color=cmap(0))

# Fixed points
z3_fp1,z5_fp1 = GetFixedPoint(amplitude_z5,ϵ3,σ0,2)
#z3_fp1      = 0.0
#z5_fp1      = sqrt((nullc_z5e(ϵ3) - (z3_fp1^2)*L5_34)/L5_56)

# Origin
z3_fp2      = 0.0
z5_fp2      = 0.0

ax4.plot(z5_fp1,z3_fp1,linestyle="none",marker="o",color="black",markersize=8)
ax4.plot(z5_fp2,z3_fp2,linestyle="none",marker="o",color="black",markersize=8,fillstyle="none",markeredgewidth=2)


# z_{3} null-cline
#-------------------------------------------------- 

z3_min      = 0.0
z3_max      = sqrt(nullc_z3e(ϵ3)/L3_34)
nstep       = 1000
# z3_zeros    = zeros(Float64,nstep)
# z5_zeros    = zeros(Float64,nstep)
# z3_line     = LinRange(0.0,0.015,nstep)
# z5_line     = LinRange(0.0,0.015,nstep)

z3_range    = LinRange(z3_min,z3_max,nstep)
z5_vals     = zeros(Float64,nstep)
for i in 1:nstep
  if i != nstep
    z5_vals[i] = sqrt((nullc_z3e(ϵ3) - (z3_range[i]^2)*L3_34)/L3_56)
  else
    z5_vals[i] = 0.0
  end
end  

ax4.plot(z5_vals,z3_range,linestyle="-",linewidth=lw,color=cmap(1),label=L"Nullcline - |z_{3}|\ LCO")
#ax4.plot(z5_line,z3_zeros,linestyle="--",linewidth=lw,color=cmap(1),label=L"Nullcline - |z_{3}|\ Fixed\ Point")
ax4.plot(z5_line,z3_zeros,linestyle="--",linewidth=lw,color=cmap(1))

# Fixed points
z3_fp3,z5_fp3 = GetFixedPoint(amplitude_z3,ϵ3,σ0,1)
#z5_fp3      = 0.0
#z3_fp3      = sqrt((nullc_z3e(ϵ3) - (z5_fp3^2)*L3_56)/L3_34)
ax4.plot(z5_fp3,z3_fp3,linestyle="none",marker="o",color="black",markersize=8,fillstyle="top",markeredgewidth=2)

ax4.set_xlabel(L"|z_{1}|",fontsize=lafs)
#ax4.set_ylabel(L"|z_{3}|",fontsize=lafs)
ax4.set_xlim(xlims)
#ax4.set_ylim(ylims)

tt2 = @sprintf("Re=%.1f",Re3)
ax4.set_title(tt2,fontsize=lafs)
#ax4.legend(fontsize=lgfs,ncols=1)

if (nullclinesave)
  svname = "bifurcation3.eps"
  h2.savefig(pdir*svname)
end  

println("Done.")







