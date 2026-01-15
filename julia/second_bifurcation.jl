# Plotting
#----------------------------------------------------------------------

# First - λ3,4 emergence
#--------------------------------------------------  
ϵ2          = find_zero(nullc_z3e,0.1)
Rei2        = Reci*(1.0 - ϵ2)
Re2         = 1.0/Rei2

# z5 is my z_{1}
# z3 is my z_{3}

z3_min      = 0.0
z3_max      = sqrt(nullc_z5e(ϵ2)/L5_34)

z3_range    = LinRange(z3_min,z3_max,nstep)
z5_vals     = zeros(Float64,nstep)
for i in 1:nstep
  z5_vals[i] = sqrt((nullc_z5e(ϵ2) - (z3_range[i]^2)*L5_34)/L5_56)
end  

ax2.plot(z5_vals,z3_range,linestyle="-",linewidth=lw,color=cmap(0),label=L"Nullcline - |z_{1}|")
ax2.plot(z5_zeros,z3_line,linestyle="--",linewidth=lw,color=cmap(0))
ax2.plot(z5_line,z3_zeros,linestyle="--",linewidth=lw,color=cmap(1),label=L"Nullcline - |z_{3}|")


# Fixed points
z3_fp1,z5_fp1 = GetFixedPoint(amplitude_z5,ϵ2,σ0,2)
#z3_fp1      = 0.0
#z5_fp1      = sqrt((nullc_z5e(ϵ2) - (z3_fp1^2)*L5_34)/L5_56)

# Origin
#z3_fp2,z5_fp2 = GetFixedPoint(amplitude_z3,ϵ2,σ0,1)
z3_fp2      = 0.0
z5_fp2      = 0.0

ax2.plot(z5_fp1,z3_fp1,linestyle="none",marker="o",color="black",markersize=8)
ax2.plot(z5_fp2,z3_fp2,linestyle="none",marker="o",color="black",markersize=8,markerfacecolor="none",markeredgewidth=2)
#ax2.legend(fontsize=lgfs)

ax2.set_xlabel(L"|z_{1}|",fontsize=lafs)
ax2.set_ylabel(L"|z_{3}|",fontsize=lafs)
ax2.set_xlim(xlims)
ax2.set_ylim(ylims)

tt2 = @sprintf("Re=%.1f",Re2)
ax2.set_title(tt2,fontsize=lafs)


if (nullclinesave)
  svname = "bifurcation2.eps"
  h2.savefig(pdir*svname)
end  

println("Done.")







