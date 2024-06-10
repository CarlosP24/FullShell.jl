using FullShell

base_model = (;
    R = 70,
    w = 0,
    d = 0,
    Δ0 = 0.23,
    ξd = 70,
    a0 = 5,
    preα = 0,
    g = 0,
    α = 0,
    μ = 0,
    τΓ = 1,
    Vexponent = 2,
    Vmin = 0,
    Vmax = 0,
    ishollow = true,
)

model = (;  base_model...,
    R = 30,
    w = 30,
    d = 10,
    Vmin = -30,
    g = 10,
    τΓ = 40,
    μ = 2,
    preα = 46.66,
    ishollow = false,
)

hSM, hSC, params = build_cyl(; model..., )
