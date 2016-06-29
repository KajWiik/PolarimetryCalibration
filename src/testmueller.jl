using PolarimetryCalibration, Plots, LaTeXStrings

PArange = -60:1.0:60
model = Float64[(müller_rec(ϵ = -0.002,
                            ϕ = 0.003,
                            ΔG = -0.002,
                            ψ = 0.11,
#                    α = fit_param[7])* # Circular feed
                            α = deg2rad(45.0))'* # Circular feed
                 müller_sky(deg2rad(PAaz))'*
                 lin_pol_source(flux = 10.11,
                        p = 0.1,
                        pa = deg2rad(-64.0)))[i] 
         for PAaz in PArange, i in 1:4]

# plot(plot(PArange, model[:,1], label = "I"), scatter!(pa, I, markersize = 3),
#      plot(PArange, model[:,2],label = "V"), scatter!(pa, V, markersize = 3),
#      plot(PArange, model[:,3],label = "U"), scatter!(pa, U, markersize = 3),
#      plot(PArange, -model[:,4],label = "Q"), scatter!(pa, Q, markersize = 3))

p1 = plot(PArange, model[:,1],title = "I")
p2 = plot(PArange, model[:,2],title = "V")
p3 = plot(PArange, model[:,3],title = "U")
p4 = plot(PArange, -model[:,4],title = "Q")
plot(p1, p2, p3, p4, legend = :none)

