using JLD

i = readdlm("i")
q = readdlm("q")
u = readdlm("u")
v = readdlm("v")
pa = readdlm("pa")

save(ARGS[1]*".jld", "i", i, "q", q, "u", u, "v", v, "pa", pa)
