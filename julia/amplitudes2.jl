using Roots
using PyPlot

# Generate the null-cline functions
include("ampfunction2.jl")

cmap  = get_cmap("tab10")

#-------------------------------------------------- 

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

#include("second_bifurcation.jl")
#include("second_third_bifurcation.jl")
#include("third_bifurcation.jl")

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

# #include("third_fourth_bifurcation.jl")
# include("re4600_bifurcation.jl")
# include("fourth_bifurcation.jl")
# include("fourth_fifth_bifurcation.jl")
# #include("fifth_bifurcation.jl")

println("Done.")




