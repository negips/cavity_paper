# Include the function files
function FloquetStepper(G::Vector{Matrix{T}},B::AbstractVector{S},V0::AbstractVector{T},Ω::AbstractVector{ComplexF64},StpInp,ArnInp) where {S,T<:Number}

  lbc             = false
  rbc             = false

  dtype           = T

  V,H             = StepperArnoldi.ArnKrylovInit(StpInp,ArnInp;Dtype=T)
  vinit           = StepperArnoldi.ArnInitVector(ArnInp.vlen,lbc,rbc,Dtype=T) 

  nev             = ArnInp.nev
  eigshift        = ArnInp.eigshift
  tol             = ArnInp.tol
  ngs             = ArnInp.ngs

  ifarnoldi       = ArnInp.ifarnoldi
  verbose         = ArnInp.ifverbose
  ifeigshift      = ArnInp.ifeigshift
  ifadjoint       = StpInp.ifadjoint         
  nsteps          = StpInp.nsteps            # Stepper Phase steps
  dt              = StpInp.timestep          # Time step

  nkryl           = 0
  lkryl           = ArnInp.lkryl
  block           = ArnInp.bsize 
  h,θ,v           = StepperArnoldi.ArnUpd(V,block,B,vinit,nkryl,ngs)
  V[:,1]          = v
  nkryl           = 0

  Rhs             = similar(v)

  ifconv          = false
  t               = 0.0           # Time
  i               = 0             # Istep

  maxouter_it     = ArnInp.outer_iterations
  major_it        = 1

  v1              = zeros(T,ArnInp.vlen) 
  v2              = zeros(T,ArnInp.vlen) 
  v3              = zeros(T,ArnInp.vlen)
  Bi              = 1.0./B          # Inverse Mass (Vector)

  # Eigenvalue Shift
  if ifeigshift
    Ωshift = exp(eigshift*nsteps*dt)
    if (verbose)
      println("EigShift: $(eigshift), Ω= $Ω, |Ω| = $(abs(Ωshift))")
    end
  else
    Ωshift = 0.0
  end

  # Start iterations
  println("Starting Stepper/Arnoldi Iterations")
  while (~ifconv)

    t = 0.0       # Time
    for i in 1:nsteps
      t = t + dt;
      # RK4!(v,L,v1,v2,v3,dt)
      FloquetRK4!(v,G,V0,Ω,t,dt)
    end  
  
    # Expand Krylov space
    V,H,nkryl,β,major_it = StepperArnoldi.IRAM!(V,H,B,v,nkryl,lkryl,major_it,nev,Ωshift,ngs,verbose)
    v   = V[:,nkryl]
  
    if (major_it>maxouter_it)
      break
    end  
  
    if (β < tol)
      @printf "Stoping Iteration, β: %12e\n" β
      ifconv = true
    end  
  end      # while ... 


  # Plotting after convergence
  Hr  = H[1:nev,1:nev]
  F   = eigen(Hr)
  evs = F.values
  DT  = dt*nsteps

  ritz = evs
  
  λr   = log.(abs.(evs))/DT
  λi   = atan.(imag(evs),real.(evs))/DT
  
  λ  = λr .+ im*λi
  # Eigenvectors  
  eigvec = V[:,1:nev]*F.vectors

  # Structured output
  arnout = StepperArnoldi.ArnoldiOutput(λ,eigvec,ritz,ifconv,tol) 
  
  return arnout
end
#---------------------------------------------------------------------- 


