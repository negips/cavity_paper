# Plotting the Spectrum

using CSV
using DataFrames
using PyPlot

close("all")

lafs  = 16
cm    = get_cmap("tab10")

pdir  = "plots/"

fname = "eigenvalues.txt"

f     = CSV.File(fname;delim=" ",ignorerepeated=true,stripwhitespace=true,skipto=2,header=false)
df    = DataFrame(f)

nm    = 21
λr    = zeros(Float64,nm)
λi    = zeros(Float64,nm)

λr[1:nm]    = df[1:nm,4]
λi[1:nm]    = df[1:nm,5]
# λr[22]      = 0.0       # Parameter Mode
# λi[22]      = 0.0       # Parameter Mode

h1  = figure(1,figsize=[7.0,5.0],layout="constrained")
ax1 = h1.subplots()
ax1.plot(λi,λr,linestyle="none",marker="o",markersize=6)
ax1.set_ylabel(L"\mathfrak{R}(λ)",fontsize=lafs)
ax1.set_xlabel(L"\mathfrak{Im}(λ)",fontsize=lafs)

ax1.set_xlim([-12,12])
ax1.set_ylim([-0.15,0.02])
ax1.plot([-12; 12],[0; 0],linewidth=0.5,linestyle="--",color="gray")

i1 = argmin(abs.(λi .- 7.5))
i2 = argmin(abs.(λi .- (-7.5)))

i3 = argmin(abs.(λi .- 10.0))
i4 = argmin(abs.(λi .- (-10.0)))

# First pair
λr1 = λr[i1]
λi1 = λi[i1]
λr2 = λr[i2]
λi2 = λi[i2]

# Second pair
λr3 = λr[i3]
λi3 = λi[i3]
λr4 = λr[i4]
λi4 = λi[i4]

ax1.plot(λi3,λr3,linestyle="none",marker="o",markersize=20,color="gray",markerfacecolor="none")
ax1.plot(λi4,λr4,linestyle="none",marker="o",markersize=20,color="gray",markerfacecolor="none")

# Annotation
#props=Dict("facecolor"=>"black", "arrowstyle"=>"->", "connectionstyle"=>"arc3,rad=0.2")
props=Dict("facecolor"=>"black", "arrowstyle"=>"->","linewidth"=>2,"color"=>cm(1))

δr = 0.005
ax1.annotate("",xy=(λi3, 0),xytext=(λi3, λr3+δr),arrowprops=props)
ax1.annotate("",xy=(λi4, 0),xytext=(λi4, λr4+δr),arrowprops=props)

svname = "cavity_spectra.eps"
h1.savefig(pdir*svname)





