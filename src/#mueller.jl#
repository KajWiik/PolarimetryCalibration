type Linear end
type Circular end
type General end
type Permuted end

function müller_rec(; ϵ::Float64 = 0.0, ϕ::Float64 = 0.0, ΔG::Float64 = 0.0, α::Float64 = 45.0*π/180, ψ::Float64 = 0.0)
    return [1.0         (-2ϵ*sin(ϕ)*sin(2α) + ΔG/2*cos(2α)) 2ϵ*cos(ϕ) (2ϵ*sin(ϕ)*cos(2α) + ΔG/2*sin(2α))
            ΔG/2                     cos(2α)                   0.0                  sin(2α)
            2ϵ*cos(ϕ + ψ)          sin(2α)*sin(ψ)             cos(ψ)            -cos(2α)*sin(ψ)
            2ϵ*sin(ϕ + ψ)         -sin(2α)*cos(ψ)             sin(ψ)             cos(2α)*cos(ψ)]
end
function müller_rec(::Permuted; ϵ::Float64 = 0.0, ϕ::Float64 = 0.0, ΔG::Float64 = 0.0, α::Float64 = 45.0*π/180, ψ::Float64 = 0.0)
    return [1.0         (-2ϵ*sin(ϕ)*sin(2α) + ΔG/2*cos(2α)) 2ϵ*cos(ϕ) (2ϵ*sin(ϕ)*cos(2α) + ΔG/2*sin(2α))
            -2ϵ*sin(ϕ + ψ)         sin(2α)*cos(ψ)             -sin(ψ)             -cos(2α)*cos(ψ)
            2ϵ*cos(ϕ + ψ)          sin(2α)*sin(ψ)             cos(ψ)            -cos(2α)*sin(ψ)
            ΔG/2                     cos(2α)                   0.0                  sin(2α)]

end

function müller_rec(::Linear; ϵ::Float64 = 0.0, ϕ::Float64 = 0.0, ΔG::Float64 = 0.0, δα::Float64 = 0.0, ψ::Float64 = 0.0, s = +)
    pm = s(1.0)
    mp = -pm
    return [1.0           pm*ΔG/2       2ϵ*cos(ϕ) 2ϵ*sin(ϕ)
            ΔG/2          pm            0.0       pm*2δα
            2ϵ*cos(ϕ + ψ) pm*2δα*sin(ψ) cos(ψ)    mp*sin(ψ)
            2ϵ*sin(ϕ + ψ) mp*2δα*cos(ψ) sin(ψ)    pm*cos(ψ)]
end

function müller_rec(::Circular; ϵ::Float64 = 0.0, ϕ::Float64 = 0.0, ΔG::Float64 = 0.0, δα::Float64 = 0.0, ψ::Float64 = 0.0, s = +)
    pm = s(1.0)
    mp = -pm
    return [1.0           mp*2ϵ*sin(ϕ)  2ϵ*cos(ϕ)  pm*ΔG/2
            ΔG/2          mp*2δα          0.0        pm
            2ϵ*cos(ϕ + ψ) pm*sin(ψ)     cos(ψ)   pm*2δα*sin(ψ)
            2ϵ*sin(ϕ + ψ) mp*cos(ψ)     sin(ψ)   mp*2δα*cos(ψ)]
end

function müller_sky(PAaz::Float64)
    cpa = cos(2*PAaz)
    spa = sin(2*PAaz)
    return [1.0   0.0   0.0  0.0
            0.0   cpa   spa  0.0
            0.0  -spa   cpa  0.0
            0.0   0.0   0.0  1.0]
end

function lin_pol_source(;flux::Float64 = 1.0, p::Float64 = 0.1, pa::Float64 = 0.0)
    return Float64[flux
                   p*flux*cos(2*pa)
                   p*flux*sin(2*pa)
                   0.0]
end
