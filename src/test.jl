using Plots

pa = readdlm("pa")
pa = reshape(pa', length(pa))
i = readdlm("i")
i = reshape(i', length(i))
q = readdlm("q")
q = reshape(q', length(q))
u = readdlm("u")
u = reshape(u', length(u))
v = readdlm("v")
v = reshape(v', length(v))

plot(pa,i)
