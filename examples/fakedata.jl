using PolarimetryCalibration, Plots, LaTeXStrings

σ = 0.01
PAsrc = 36.0
PArange = 0:0.2:180

stokes = [(müller_rec()*
           müller_sky(PAaz*π/180)*
           lin_pol_source(pa = PAsrc*π/180))[i] + σ*randn()
          for PAaz = PArange, i in 1:4]
plot(PArange, stokes[:,2], label=L"ideal", xlabel = "PA", ylabel = "Stokes Q")

stokes = [(müller_rec(ΔG = 0.1)*
           müller_sky(PAaz*π/180)*
           lin_pol_source(pa = PAsrc*π/180))[i] + σ*randn()
          for PAaz = PArange, i in 1:4]
plot!(PArange, stokes[:,2], label=L"\Delta G = 0.1")

stokes = [(müller_rec(α = 0.5, ψ = 0.9)*
           müller_sky(PAaz*π/180)*
           lin_pol_source(pa = PAsrc*π/180))[i] + σ*randn()
          for PAaz = PArange, i in 1:4]
plot!(PArange, stokes[:,2], label=L"\alpha = 0.5, \psi = 0.9")
