using PolarimetryCalibration, Plots, LaTeXStrings

σ = 0.01
PAsrc = 36.0
PArange = 0:0.2:180
#const (I,Q,U,V)=(1,2,3,4)

stokes = [(müller_rec()*
           müller_sky(PAaz*π/180)*
           lin_pol_source(pa = PAsrc*π/180))[i] + σ*randn()
          for PAaz = PArange, i in 1:4]

plot(PArange, stokes[:,1],
     label="I",
     xlabel = "PA",
     ylabel = "Jy",
     title = L"$\Delta G = 0.1, \psi = 52$ $^\ocirc, PAsrc = 36$ $^\ocirc, p = 10\%, S = 1 Jy, \sigma = 10 mJy$")
#     title = L"$\Delta G = 0.1$, $\psi = 52$ $^\ocirc$, PAsrc = 36$ $^\ocirc$, p = 10\%, S = 1 Jy, $\sigma =$ 10 mJy")

plot!(PArange, stokes[:,2], label="V")
plot!(PArange, stokes[:,3], label="U")
plot!(PArange, -stokes[:,4], label="Q")

